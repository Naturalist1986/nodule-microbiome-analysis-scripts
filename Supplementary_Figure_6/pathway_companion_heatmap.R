#!/usr/bin/env Rscript
# ============================================================================
# PATHWAY/MODULE COMPANION HEATMAP
# ============================================================================
# Uses Boruta KOs as seeds to find the rest of their pathway/module members
# in normalized_kegg_results.tsv, tests whether companion KOs also correlate
# with the same plant trait, and draws a heatmap clustered by pathway/module.
# ============================================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

cat("================================================================================\n")
cat("PATHWAY / MODULE COMPANION HEATMAP\n")
cat("================================================================================\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

SUMMARY_FILE       <- "Boruta_KO_heatmap_with_abundance_summary.csv"
GENE_NAMES_FILE    <- "ko_gene_names.csv"
KO_ABUNDANCE_FILE  <- "normalized_kegg_results.tsv"
MEASUREMENTS_FILE  <- "Guy_measurements.xlsx"
KO_BINS_FILE       <- "ko_presence_absence_table.xlsx"    # bin presence/absence
KO_MODULE_CACHE    <- "ko_pathway_module_cache.rds"      # KO -> pathways/modules
PW_KO_CACHE        <- "pathway_to_kos_cache.rds"          # pathway/module -> all KOs
KEGG_DELAY         <- 0.25   # seconds between API calls

# Overview pathways to skip (too large/generic to be informative)
SKIP_PATHWAYS <- c(
  "map01100",   # Metabolic pathways          (4970 KOs)
  "map01120",   # Microbial metabolism         (1437 KOs)
  "map01110",   # Biosynthesis of secondary metabolites (2307 KOs)
  "map01250",   # Biosynthesis of nucleotide sugars
  "map01232",   # Nucleotide sugar biosynthesis
  "map01053",   # Biosynthesis of siderophore group nonribosomal peptides
  "map01310",   # Biosynthesis of terpenoids and polyketides  (sometimes added)
  "map03320"    # PPAR signaling (not relevant for bacteria)
)

# Modules to exclude from the analysis
SKIP_MODULES <- c(
  "M00897"    # Thiamine biosynthesis (thiN/K00949 and companions removed)
)

# ============================================================================
# 1.  LOAD BORUTA KO SUMMARY
# ============================================================================

cat("STEP 1: Loading Boruta KO summary...\n")

boruta_df  <- read_csv(SUMMARY_FILE,    show_col_types = FALSE)
gene_names <- read_csv(GENE_NAMES_FILE, show_col_types = FALSE)

boruta_kos <- boruta_df$KO
cat(sprintf("  %d Boruta KOs\n", length(boruta_kos)))

# ============================================================================
# 1b. IDENTIFY BORUTA KOs PRESENT IN AT LEAST ONE GENOMIC BIN
#     Mirrors the logic in Boruta_KO_Heatmap_with_abundance.R (Phase 11)
# ============================================================================

cat("\nSTEP 1b: Filtering seed KOs to those present in at least one genomic bin...\n")

ko_bins_raw <- read_excel(KO_BINS_FILE)

# Transpose if KOs are rows (same handling as Boruta_KO_Heatmap_with_abundance.R)
if ("KO" %in% colnames(ko_bins_raw)) {
  ko_col_vals      <- ko_bins_raw$KO
  ko_bins_matrix   <- as.matrix(ko_bins_raw[, -1])
  rownames(ko_bins_matrix) <- ko_col_vals
  ko_bins_t        <- t(ko_bins_matrix)
  ko_bins          <- as.data.frame(ko_bins_t)
  ko_bins$bin      <- rownames(ko_bins_t)
  ko_bins          <- ko_bins[, c("bin", setdiff(names(ko_bins), "bin"))]
} else {
  ko_bins <- as.data.frame(ko_bins_raw)
}

# All KOs (any KO ID) present (value == 1) in at least one bin
bin_long <- ko_bins %>%
  pivot_longer(cols = -bin, names_to = "KO", values_to = "present") %>%
  filter(present == 1)

all_kos_in_any_bin <- unique(bin_long$KO)   # used to filter ALL heatmap rows

# Subset: Boruta seeds present in at least one bin (used to filter seed_long)
kos_in_any_bin <- intersect(boruta_kos, all_kos_in_any_bin)

cat(sprintf("  %d / %d Boruta seed KOs detected in at least one genomic bin\n",
            length(kos_in_any_bin), length(boruta_kos)))
cat(sprintf("  Excluded seeds (not in any bin): %s\n",
            paste(setdiff(boruta_kos, kos_in_any_bin), collapse = ", ")))

# Map KO -> label (gene name or KO id) — starts with the 67 Boruta KOs
ko_label_map <- setNames(
  ifelse(!is.na(gene_names$Gene_Name) & gene_names$Gene_Name != "",
         gene_names$Gene_Name, gene_names$KO),
  gene_names$KO
)

# Fetch gene names for companion KOs not already covered.
# Called lazily after row_df is built; uses a persistent cache file.
ALL_KO_NAMES_CACHE <- "all_ko_gene_names_cache.rds"
if (file.exists(ALL_KO_NAMES_CACHE)) {
  all_ko_name_cache <- readRDS(ALL_KO_NAMES_CACHE)
} else {
  all_ko_name_cache <- list()
}

fetch_ko_short_name <- function(ko_id) {
  if (!is.null(all_ko_name_cache[[ko_id]])) return(all_ko_name_cache[[ko_id]])
  nm <- tryCatch({
    url_n <- paste0("https://rest.kegg.jp/list/", ko_id)
    con   <- url(url_n); on.exit(close(con), add = TRUE)
    txt   <- readLines(con, warn = FALSE)
    if (length(txt) == 0) return(ko_id)
    # Format: "ko:K00001\tgene1, gene2; full description"
    parts <- strsplit(txt[1], "\t")[[1]]
    if (length(parts) < 2) return(ko_id)
    # Extract gene name field (everything before first ";")
    gene_field <- trimws(sub(";.*", "", parts[2]))
    # Split by comma — KEGG separates multiple gene names with ", "
    candidates <- trimws(strsplit(gene_field, ",")[[1]])
    # Keep genuine short names: <=12 chars, no spaces, not an EC/bracket entry
    valid <- candidates[
      nchar(candidates) <= 12 &
      !grepl(" ", candidates) &
      !grepl("^\\[", candidates) &
      nchar(candidates) > 0
    ]
    if (length(valid) > 0) paste(valid, collapse = "/") else ko_id
  }, error = function(e) ko_id)
  all_ko_name_cache[[ko_id]] <<- nm
  Sys.sleep(KEGG_DELAY)
  nm
}

enrich_ko_labels <- function(kos) {
  # Re-try KOs whose cached value is just their own KO ID (previous failed lookup)
  # — they may now resolve correctly with the multi-name logic
  stale <- names(all_ko_name_cache)[
    unlist(all_ko_name_cache) == names(all_ko_name_cache)
  ]
  all_ko_name_cache[stale] <<- NULL

  missing <- setdiff(kos, names(ko_label_map))
  to_fetch <- union(missing, intersect(kos, stale))
  if (length(to_fetch) == 0) return(invisible(NULL))
  cat(sprintf("  Fetching gene names for %d KOs...\n", length(to_fetch)))
  for (ko in to_fetch) {
    nm <- fetch_ko_short_name(ko)
    ko_label_map[ko] <<- nm
  }
  saveRDS(all_ko_name_cache, ALL_KO_NAMES_CACHE)
  cat("  Gene name cache saved.\n")
}

get_label <- function(ko) {
  lbl <- ko_label_map[ko]
  ifelse(is.na(lbl), ko, lbl)
}

# Map KO -> primary trait
ko_trait_map <- setNames(boruta_df$Traits, boruta_df$KO)

# ============================================================================
# 2.  LOAD KEGG MEMBERSHIP CACHE  (KO → pathways + modules)
# ============================================================================

cat("\nSTEP 2: Loading KO-pathway/module cache...\n")

if (!file.exists(KO_MODULE_CACHE)) {
  stop(paste("Cache file not found:", KO_MODULE_CACHE,
             "\nRun map_KOs_to_pathways_modules.R first."))
}
ko_links <- readRDS(KO_MODULE_CACHE)

# Build long table: seed_KO  |  type (pathway/module)  |  id
seed_long <- bind_rows(lapply(names(ko_links), function(ko) {
  l <- ko_links[[ko]]
  # pathways (skip overview)
  pw  <- setdiff(l$pathways, SKIP_PATHWAYS)
  mod <- setdiff(l$modules,  SKIP_MODULES)
  rows <- NULL
  if (length(pw)  > 0) rows <- rbind(rows, data.frame(KO=ko, type="pathway", id=pw,  stringsAsFactors=FALSE))
  if (length(mod) > 0) rows <- rbind(rows, data.frame(KO=ko, type="module",  id=mod, stringsAsFactors=FALSE))
  rows
}))

cat(sprintf("  Seed entries after filtering overview pathways: %d\n", nrow(seed_long)))

# Remove seed KOs that were not detected in any genomic bin
seed_long <- seed_long %>% filter(KO %in% kos_in_any_bin)
cat(sprintf("  Seed entries after removing bin-absent seeds:    %d\n", nrow(seed_long)))

# ============================================================================
# 3.  LOAD ABUNDANCE DATA AND BUILD MERGED SAMPLE MATRIX
# ============================================================================

cat("\nSTEP 3: Loading abundance data...\n")

ko_data      <- read_tsv(KO_ABUNDANCE_FILE, show_col_types = FALSE)
measurements <- read_excel(MEASUREMENTS_FILE)
measurements$sample_id <- paste(measurements$Treatment, measurements$Replicate, sep = "_")

all_ko_ids    <- ko_data$kegg_number
all_sample_cols <- colnames(ko_data)[-1]

# Extract treatment + replicate from column name
extract_sample_info <- function(nm) {
  clean <- gsub("_DKDN.*", "", nm)
  parts <- strsplit(clean, "_")[[1]]
  if (length(parts) >= 2)
    return(c(parts[1], parts[2], paste(parts[1], parts[2], sep = "_")))
  c(NA, NA, NA)
}

sample_info <- data.frame(
  original_name = all_sample_cols,
  t(sapply(all_sample_cols, extract_sample_info)),
  stringsAsFactors = FALSE
)
colnames(sample_info) <- c("original_name", "treatment", "replicate", "sample_id")

# Build full abundance matrix and merge technical replicates
ko_matrix_all <- as.matrix(ko_data[, -1])
rownames(ko_matrix_all) <- all_ko_ids

unique_samples <- unique(sample_info$sample_id)
merged_matrix  <- matrix(NA, nrow = nrow(ko_matrix_all), ncol = length(unique_samples))
rownames(merged_matrix) <- rownames(ko_matrix_all)
colnames(merged_matrix) <- unique_samples

for (sid in unique_samples) {
  cols <- sample_info$original_name[sample_info$sample_id == sid]
  if (length(cols) == 1) {
    merged_matrix[, sid] <- ko_matrix_all[, cols]
  } else {
    merged_matrix[, sid] <- rowMeans(ko_matrix_all[, cols], na.rm = TRUE)
  }
}

# Align with measurements
common_samples  <- intersect(colnames(merged_matrix), measurements$sample_id)
merged_matrix   <- merged_matrix[, common_samples, drop = FALSE]
traits_filtered <- measurements %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

cat(sprintf("  Abundance matrix: %d KOs x %d samples\n",
            nrow(merged_matrix), ncol(merged_matrix)))

# ============================================================================
# 4.  QUERY KEGG: ALL KOs IN EACH PATHWAY/MODULE  (with cache)
# ============================================================================

cat("\nSTEP 4: Fetching all KOs for each pathway/module from KEGG...\n")

fetch_kos_for_pathway <- function(path_id) {
  # Works for both map-pathways and M-modules
  type_str <- if (grepl("^M", path_id)) "module" else "pathway"
  query_id <- if (grepl("^M", path_id)) paste0("md:", path_id) else paste0("path:", path_id)
  url_str  <- paste0("https://rest.kegg.jp/link/ko/", path_id)
  tryCatch({
    con <- url(url_str)
    on.exit(close(con), add = TRUE)
    txt <- readLines(con, warn = FALSE)
    # Format: "path:map00xxx\tko:K00001"
    ko_entries <- grep("ko:K", txt, value = TRUE)
    kos <- sub(".*ko:", "", ko_entries)
    kos
  }, error = function(e) character(0))
}

unique_ids <- unique(seed_long$id)

if (file.exists(PW_KO_CACHE)) {
  cat("  Loading cached pathway->KO mappings...\n")
  pw_ko_cache <- readRDS(PW_KO_CACHE)
} else {
  pw_ko_cache <- list()
}

ids_to_fetch <- setdiff(unique_ids, names(pw_ko_cache))
if (length(ids_to_fetch) > 0) {
  cat(sprintf("  Querying KEGG for %d new pathway/module IDs...\n", length(ids_to_fetch)))
  for (i in seq_along(ids_to_fetch)) {
    pid <- ids_to_fetch[i]
    cat(sprintf("  [%d/%d] %s\n", i, length(ids_to_fetch), pid))
    pw_ko_cache[[pid]] <- fetch_kos_for_pathway(pid)
    Sys.sleep(KEGG_DELAY)
  }
  saveRDS(pw_ko_cache, PW_KO_CACHE)
  cat("  Cache saved.\n")
} else {
  cat("  All pathways/modules already cached.\n")
}

# ============================================================================
# 5.  BUILD KO × PATHWAY/MODULE TABLE  — filter to KOs in abundance data
# ============================================================================

cat("\nSTEP 5: Building companion KO table...\n")

# Build long table: for each (seed_KO, pathway/module), list ALL member KOs
# that are present in the abundance dataset
rows_list <- lapply(seq_len(nrow(seed_long)), function(i) {
  seed_ko <- seed_long$KO[i]
  pw_id   <- seed_long$id[i]
  pw_type <- seed_long$type[i]
  all_member_kos <- pw_ko_cache[[pw_id]]
  if (length(all_member_kos) == 0) return(NULL)
  present_kos <- intersect(all_member_kos, rownames(merged_matrix))
  if (length(present_kos) < 2) return(NULL)   # need seed + >=1 companion
  data.frame(
    seed_KO    = seed_ko,
    seed_trait = unname(ko_trait_map[seed_ko]),
    pw_id      = pw_id,
    pw_type    = pw_type,
    KO         = present_kos,
    total_pw_kos = length(all_member_kos),
    present_pw_kos = length(present_kos),
    stringsAsFactors = FALSE
  )
})

companion_long <- bind_rows(rows_list)

if (nrow(companion_long) == 0) {
  cat("  No pathways/modules have >=2 KOs present in the data. Exiting.\n")
  quit(save = "no", status = 0)
}

cat(sprintf("  Pathway/module groups with >=2 KOs present: %d\n",
            n_distinct(companion_long$pw_id)))
cat(sprintf("  Total KO × pathway/module entries: %d\n", nrow(companion_long)))

# Flag whether each KO is a Boruta seed KO or a companion
companion_long <- companion_long %>%
  mutate(is_boruta = KO %in% boruta_kos)

# ============================================================================
# 6.  FETCH PATHWAY/MODULE NAMES  (reuse existing name cache if present)
# ============================================================================

cat("\nSTEP 6: Fetching pathway/module names...\n")

NAME_CACHE_FILE <- "pathway_module_names_cache.rds"
if (file.exists(NAME_CACHE_FILE)) {
  name_cache <- readRDS(NAME_CACHE_FILE)
} else {
  name_cache <- list()
}

get_name <- function(id) {
  if (!is.null(name_cache[[id]])) return(name_cache[[id]])
  url_n <- paste0("https://rest.kegg.jp/list/", id)
  nm <- tryCatch({
    con <- url(url_n); on.exit(close(con), add = TRUE)
    txt <- readLines(con, warn = FALSE)
    if (length(txt) > 0) {
      parts <- strsplit(txt[1], "\t")[[1]]
      if (length(parts) >= 2) parts[2] else id
    } else id
  }, error = function(e) id)
  name_cache[[id]] <<- nm
  Sys.sleep(KEGG_DELAY)
  nm
}

ids_need_name <- setdiff(unique(companion_long$pw_id), names(name_cache))
if (length(ids_need_name) > 0) {
  cat(sprintf("  Fetching %d new names from KEGG...\n", length(ids_need_name)))
  for (id in ids_need_name) get_name(id)
  saveRDS(name_cache, NAME_CACHE_FILE)
}

pw_name_df <- data.frame(
  pw_id   = unique(companion_long$pw_id),
  pw_name = sapply(unique(companion_long$pw_id), function(id) {
    nm <- name_cache[[id]]; if (is.null(nm)) id else nm
  }),
  stringsAsFactors = FALSE
)

# Trim long names for display
pw_name_df <- pw_name_df %>%
  mutate(pw_label = ifelse(nchar(pw_name) > 45,
                           paste0(substr(pw_name, 1, 42), "..."),
                           pw_name),
         pw_label = paste0(pw_id, ": ", pw_label))

companion_long <- companion_long %>%
  left_join(pw_name_df, by = "pw_id")

# ============================================================================
# 7.  COMPUTE SPEARMAN CORRELATIONS WITH THE SEEDED TRAIT
# ============================================================================

cat("\nSTEP 7: Computing Spearman correlations...\n")

trait_cols <- c("percent_N", "Plant_Biomass", "Fixation_per_Nodule", "NMF")

# For every unique KO present, compute r and p for ALL traits
ko_cor_df <- data.frame(KO = intersect(unique(companion_long$KO), rownames(merged_matrix)))

for (trait in trait_cols) {
  trait_vals <- traits_filtered[[trait]]
  ko_cor_df[[paste0("r_", trait)]] <- sapply(ko_cor_df$KO, function(ko) {
    x <- merged_matrix[ko, ]
    valid <- !is.na(x) & !is.na(trait_vals)
    if (sum(valid) >= 3) cor(x[valid], trait_vals[valid], method = "spearman") else NA_real_
  })
  ko_cor_df[[paste0("p_", trait)]] <- sapply(ko_cor_df$KO, function(ko) {
    x <- merged_matrix[ko, ]
    valid <- !is.na(x) & !is.na(trait_vals)
    if (sum(valid) >= 3) cor.test(x[valid], trait_vals[valid], method = "spearman", exact = FALSE)$p.value else NA_real_
  })
}

# For each (KO, pw_id) entry, pull out the r/p for the seeded trait
companion_long <- companion_long %>%
  left_join(ko_cor_df, by = "KO") %>%
  mutate(
    r_seeded = case_when(
      seed_trait == "percent_N"          ~ r_percent_N,
      seed_trait == "Plant_Biomass"      ~ r_Plant_Biomass,
      seed_trait == "Fixation_per_Nodule"~ r_Fixation_per_Nodule,
      seed_trait == "NMF"                ~ r_NMF,
      TRUE ~ NA_real_
    ),
    p_seeded = case_when(
      seed_trait == "percent_N"          ~ p_percent_N,
      seed_trait == "Plant_Biomass"      ~ p_Plant_Biomass,
      seed_trait == "Fixation_per_Nodule"~ p_Fixation_per_Nodule,
      seed_trait == "NMF"                ~ p_NMF,
      TRUE ~ NA_real_
    )
  )

# Significance stars
sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}
companion_long$sig_seeded <- sig_stars(companion_long$p_seeded)

# ============================================================================
# 8.  SELECT ONE pw_id per KO for row grouping (prefer module > pathway)
# ============================================================================

cat("\nSTEP 8: Selecting representative pathway/module per KO...\n")

# For KOs that belong to multiple seed groups, pick one:
# Priority: module > pathway; within same type, pick first by pw_id (alphabetical)
row_df <- companion_long %>%
  group_by(KO) %>%
  arrange(KO,
          desc(pw_type == "module"),    # modules first
          pw_id) %>%                    # then alphabetically by id
  slice(1) %>%
  ungroup() %>%
  arrange(seed_trait, pw_type, pw_id, desc(is_boruta), desc(abs(r_seeded)))

# For each selected group, count companions vs boruta seeds
group_stats <- companion_long %>%
  group_by(pw_id, pw_label, pw_type) %>%
  summarise(
    n_total     = n_distinct(KO),
    n_boruta    = sum(is_boruta & !duplicated(KO)),
    n_companion = n_total - n_boruta,
    trait       = paste(sort(unique(seed_trait)), collapse = "/"),
    .groups     = "drop"
  ) %>%
  filter(n_companion >= 1,    # must have at least 1 non-Boruta companion
         n_total     <= 15)   # skip large/generic pathways with >15 present KOs

# Keep only rows for groups that have companions
row_df <- row_df %>%
  filter(pw_id %in% group_stats$pw_id)

# ── Bin-presence filter (second pass) ────────────────────────────────────────
# Remove ANY row (Boruta seed or companion) whose KO was not detected in at
# least one genomic bin.  This catches:
#   (a) bin-absent Boruta seeds appearing as companions via another seed's group
#   (b) companion KOs (is_boruta=FALSE) that are present in the abundance data
#       but absent from all MAG bins (e.g. thiDE / K14153)
n_before <- nrow(row_df)
row_df <- row_df %>%
  filter(KO %in% all_kos_in_any_bin)
cat(sprintf("  Rows removed (not present in any genomic bin): %d\n",
            n_before - nrow(row_df)))

cat(sprintf("  Groups retained (>=1 companion KO): %d\n",  n_distinct(row_df$pw_id)))
cat(sprintf("  Total KOs in heatmap:               %d\n",  nrow(row_df)))

if (nrow(row_df) == 0) {
  cat("  No groups with companion KOs. Try relaxing SKIP_PATHWAYS. Exiting.\n")
  quit(save = "no", status = 0)
}

# Print summary
cat("\nPathway/Module groups in heatmap:\n")
cat(strrep("-", 80), "\n")
for (i in seq_len(nrow(group_stats))) {
  g <- group_stats[i, ]
  if (!g$pw_id %in% row_df$pw_id) next
  cat(sprintf("  [%s] %s\n    Trait: %s | %d Boruta KOs + %d companions (%d total)\n",
              g$pw_id, g$pw_label, g$trait,
              g$n_boruta, g$n_companion, g$n_total))
}

# ============================================================================
# 9.  SHARED SETUP  (same for both heatmaps)
# ============================================================================

cat("\nSTEP 9: Preparing shared data...\n")

# Fetch gene names for all companion KOs in the heatmap
enrich_ko_labels(unique(row_df$KO))

treatment_names <- c("ces"  = "CES", "RH"   = "RH",  "carR" = "CAR-G",
                     "carK" = "CAR-C", "hok" = "HUK", "mtz"  = "MTZ")
treatment_order <- c("CAR-C", "CAR-G", "CES", "HUK", "MTZ", "RH")

sample_treat <- traits_filtered %>%
  mutate(treat_label = treatment_names[Treatment]) %>%
  arrange(match(treat_label, treatment_order), sample_id)
ordered_samples <- sample_treat$sample_id

# Column labels: just the replicate number (treatment shown as facet title)
col_rep_labels <- sub("^[^_]+_", "", ordered_samples)

# Per-sample treatment factor used for column faceting
col_treat_vec <- factor(
  setNames(treatment_names[sample_treat$Treatment],
           sample_treat$sample_id)[ordered_samples],
  levels = treatment_order
)

# Colour scale: light grey for exactly 0 (absent), then blue→white→red for
# any detectable abundance.  The 0→0.001 step is narrower than one cell so
# it acts as a discrete absent/present boundary in the legend.
col_fun <- colorRamp2(
  c(0,        0.001,    sqrt(0.1), sqrt(0.25), sqrt(0.5), sqrt(1.0), sqrt(1.5), sqrt(2.0)),
  c("#d0d0d0","#2166AC","#4393C3", "#92C5DE",  "#f7f7f7", "#D6604D", "#B2182B", "#67001F")
)

# Fixed cell dimensions (mm) — equal for square cells
CELL_H_MM     <- 2.5
CELL_W_MM     <- 2.5
# Overhead: group strip + row names + Trait strip + Corr bar + facet gaps + right legends + margins
OVERHEAD_W_MM <- 200   # right-side legends + padding
OVERHEAD_H_MM <- 45    # column titles + sample names (no bottom legend block)

base_palette <- c(
  "#E41A1C","#377EB8","#4DAF4A","#FF7F00","#984EA3","#A65628",
  "#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3",
  "#A6D854","#FFD92F","#E5C494","#1B9E77","#D95F02","#7570B3",
  "#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#B3E2CD",
  "#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC"
)

trait_col <- c("percent_N"           = "#377EB8",
               "Plant_Biomass"       = "#4DAF4A",
               "Fixation_per_Nodule" = "#FF7F00",
               "NMF"                 = "#984EA3")

# Horizontal barplot annotation for row annotations
# bars run left (negative r) to right (positive r) within each row band
make_horiz_barplot_anno <- function(values, colors_vec, stars, xlim = c(-1, 1)) {
  function(index, k = NULL, n = NULL) {
    vals    <- values[index]
    strs    <- stars[index]
    cols    <- colors_vec[index]
    n_items <- length(vals)
    pushViewport(viewport(xscale = xlim, yscale = c(0.5, n_items + 0.5), clip = TRUE))
    # Vertical zero line
    grid.lines(x = unit(0, "native"), y = c(0, 1),
               gp = gpar(col = "grey50", lwd = 0.6))
    for (i in seq_along(vals)) {
      v     <- ifelse(is.na(vals[i]), 0, vals[i])
      v     <- max(xlim[1], min(xlim[2], v))
      col_i <- if (is.na(vals[i])) "grey90" else cols[i]
      just  <- if (v >= 0) c("left", "centre") else c("right", "centre")
      yi    <- n_items + 1 - i   # flip so top row = top of viewport
      grid.rect(x     = unit(0, "native"),
                y     = unit(yi, "native"),
                width = unit(abs(v), "native"),
                height = unit(0.7, "native"),
                just  = just,
                gp    = gpar(fill = col_i, col = NA))
      if (nchar(strs[i]) > 0) {
        sx <- if (v >= 0) unit(min(v + 0.05, xlim[2]), "native")
              else        unit(max(v - 0.05, xlim[1]), "native")
        grid.text(strs[i], x = sx, y = unit(yi, "native"),
                  gp   = gpar(fontsize = 5, fontface = "bold", col = "black"),
                  just = c(if (v >= 0) "left" else "right", "centre"))
      }
    }
    for (tx in c(-1, 0, 1))
      grid.lines(x = unit(tx, "native"), y = c(0, 0.04),
                 gp = gpar(col = "grey40", lwd = 0.5))
    popViewport()
  }
}

save_heatmap <- function(file, width, height, draw_fn) {
  pdf(paste0(file, ".pdf"), width = width, height = height)
  draw_fn(); dev.off()
  cat(sprintf("  Saved: %s.pdf\n", file))
  dpi <- 300
  png(paste0(file, ".png"), width = round(width * dpi),
      height = round(height * dpi), res = dpi)
  draw_fn(); dev.off()
  cat(sprintf("  Saved: %s.png\n", file))
  svg(paste0(file, ".svg"), width = width, height = height)
  draw_fn(); dev.off()
  cat(sprintf("  Saved: %s.svg\n", file))
}

# ============================================================================
# 10-11.  BUILD AND SAVE ONE HEATMAP PER TYPE  (module / pathway)
# ============================================================================

build_heatmap <- function(type, outfile) {

  sub_df <- row_df %>% filter(pw_type == type)
  if (nrow(sub_df) == 0) {
    cat(sprintf("  No groups of type '%s'. Skipping.\n", type)); return(invisible(NULL))
  }
  cat(sprintf("\n--- %s heatmap: %d KOs in %d groups ---\n",
              toupper(type), nrow(sub_df), n_distinct(sub_df$pw_id)))

  # Rows ordered by module/pathway group, seeds first within each group
  kos <- sub_df$KO

  # Abundance matrix: per-genome-normalized values (hits/genome), no z-scoring
  mat <- merged_matrix[kos, ordered_samples, drop = FALSE]
  # Sqrt-transform to emphasise differences near 0; clamp raw values first
  mat <- sqrt(pmin(mat, 2))

  # Row labels: short gene name (bold = Boruta seed; plain = companion)
  rlabels      <- setNames(sapply(kos, get_label), kos)
  row_fontface <- ifelse(sub_df$is_boruta, "bold", "plain")

  # Group color strip (left annotation)
  gids    <- unique(sub_df$pw_id)
  gcols   <- setNames(base_palette[seq_along(gids)], gids)
  glabels <- setNames(sapply(gids, function(id) {
    nm  <- name_cache[[id]]
    lbl <- if (is.null(nm) || nm == id) id else nm
    # keep full name — legend width will auto-expand to fit
    paste0(id, ": ", lbl)
  }), gids)
  row_glabel <- glabels[sub_df$pw_id]

  left_ha <- rowAnnotation(
    Group = row_glabel,
    col   = list(Group = setNames(gcols, glabels)),
    annotation_name_side = "top",
    annotation_name_gp   = gpar(fontsize = 6),
    simple_anno_size     = unit(3, "mm"),
    show_legend          = TRUE,
    annotation_legend_param = list(
      Group = list(
        title     = if (type == "module") "Module" else "Pathway",
        title_gp  = gpar(fontsize = 7, fontface = "bold"),
        labels_gp = gpar(fontsize = 5),
        nrow      = length(gids)
      )
    )
  )

  # Right annotation: Trait color strip + horizontal Corr barplot
  sr     <- sub_df$r_seeded
  sstars <- sub_df$sig_seeded
  scol   <- ifelse(!is.na(sr) & sr > 0, "#E41A1C",
            ifelse(!is.na(sr) & sr < 0, "#377EB8", "grey70"))

  row_ha <- rowAnnotation(
    Trait = sub_df$seed_trait,
    Corr  = AnnotationFunction(
      fun   = make_horiz_barplot_anno(sr, scol, sstars),
      which = "row",
      width = unit(1, "cm")
    ),
    col  = list(Trait = trait_col),
    annotation_name_side = "top",
    annotation_name_gp   = gpar(fontsize = 6),
    annotation_width     = list(
      Trait = unit(0.5, "mm"),
      Corr  = unit(2, "cm")
    )
  )

  n_rows    <- nrow(mat)
  n_cols    <- ncol(mat)
  body_w_mm <- n_cols * CELL_W_MM
  body_h_mm <- n_rows * CELL_H_MM
  PAGE_W    <- (body_w_mm + OVERHEAD_W_MM) / 25.4
  PAGE_H    <- (body_h_mm + OVERHEAD_H_MM) / 25.4
  cat(sprintf("  %d rows x %d cols  |  page %.1f x %.1f in\n",
              n_rows, n_cols, PAGE_W, PAGE_H))

  ht <- Heatmap(
    mat,
    name                   = "Abundance\n(hits/genome)",
    col                    = col_fun,
    # No clustering — rows arranged by module membership
    cluster_rows           = FALSE,
    show_row_dend          = FALSE,
    # Treatment facets
    column_split           = col_treat_vec,
    cluster_column_slices  = FALSE,
    cluster_columns        = FALSE,
    column_title_gp        = gpar(fontsize = 8, fontface = "bold"),
    column_title_side      = "top",
    column_gap             = unit(2, "mm"),
    # Row names on LEFT so layout is: [group] [gene name] [heatmap] [Trait][Corr]
    show_row_names         = TRUE,
    row_names_side         = "left",
    row_names_gp           = gpar(fontsize = 5, fontface = row_fontface),
    row_labels             = rlabels[kos],
    # Sample names below each facet: replicate number only, larger font
    show_column_names      = TRUE,
    column_labels          = col_rep_labels,
    column_names_gp        = gpar(fontsize = 7),
    column_names_rot       = 0,
    top_annotation         = NULL,   # treatment shown as facet titles
    left_annotation        = left_ha,
    right_annotation       = row_ha,
    width                  = unit(body_w_mm, "mm"),
    height                 = unit(body_h_mm, "mm"),
    heatmap_legend_param   = list(
      title_gp         = gpar(fontsize = 7, fontface = "bold"),
      labels_gp        = gpar(fontsize = 6),
      legend_direction = "vertical",
      title_position   = "topleft",
      at               = c(0, sqrt(0.1), sqrt(0.25), sqrt(0.5), sqrt(1.0), sqrt(1.5), sqrt(2.0)),
      labels           = c("0", "0.1", "0.25", "0.5", "1.0", "1.5", "2.0")
    ),
    border = FALSE
  )

  save_heatmap(outfile, PAGE_W, PAGE_H, function() {
    draw(ht,
         merge_legend            = FALSE,
         heatmap_legend_side     = "right",
         annotation_legend_side  = "right",
         padding                 = unit(c(4, 10, 15, 30), "mm"))
  })
}

cat("\nSTEP 10-11: Building and saving heatmaps...\n")
build_heatmap("module",  "pathway_companion_heatmap_modules")
build_heatmap("pathway", "pathway_companion_heatmap_pathways")

# ============================================================================
# 12.  SAVE COMPANION KO TABLE
# ============================================================================

cat("\nSTEP 12: Saving companion KO table...\n")

output_table <- companion_long %>%
  filter(pw_id %in% group_stats$pw_id) %>%
  mutate(gene_label = sapply(KO, get_label)) %>%
  select(pw_id, pw_name, pw_type, seed_KO, seed_trait,
         KO, gene_label, is_boruta,
         r_seeded, p_seeded, sig_seeded,
         r_percent_N, p_percent_N,
         r_Plant_Biomass, p_Plant_Biomass,
         r_Fixation_per_Nodule, p_Fixation_per_Nodule,
         r_NMF, p_NMF,
         total_pw_kos, present_pw_kos) %>%
  arrange(seed_trait, pw_type, pw_id, desc(is_boruta), desc(abs(r_seeded)))

write_csv(output_table, "pathway_companion_kos.csv")
cat("  Saved: pathway_companion_kos.csv\n")

# ============================================================================
# PRINT SUMMARY TABLE
# ============================================================================

cat("\n\nSUMMARY: KOs per group with significant correlation (p<0.05) to seeded trait\n")
cat(strrep("=", 80), "\n")

summary_out <- companion_long %>%
  filter(pw_id %in% group_stats$pw_id) %>%
  group_by(pw_id, pw_name, pw_type) %>%
  summarise(
    n_total     = n_distinct(KO),
    n_boruta    = sum(is_boruta & !duplicated(KO)),
    n_sig       = sum(!duplicated(KO) & !is.na(p_seeded) & p_seeded < 0.05),
    trait       = paste(sort(unique(seed_trait)), collapse = "/"),
    KOs_sig     = paste(KO[!is.na(p_seeded) & p_seeded < 0.05 & !duplicated(KO)],
                        collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sig), desc(n_total))

for (i in seq_len(nrow(summary_out))) {
  r <- summary_out[i, ]
  cat(sprintf("[%s] %s\n  Type: %s | Trait: %s | Total present: %d | Boruta seeds: %d | Sig (p<0.05): %d\n",
              r$pw_id, r$pw_name, r$pw_type, r$trait,
              r$n_total, r$n_boruta, r$n_sig))
  if (nchar(r$KOs_sig) > 0)
    cat(sprintf("  Sig KOs: %s\n", r$KOs_sig))
  cat("\n")
}

# ============================================================================
# BORUTA GENE MODULE SUMMARY
# For each original Boruta gene: list every module/pathway, total KEGG KOs,
# KOs present in the abundance data, and whether it was retained in the heatmap.
# ============================================================================

cat("\n\nSUMMARY: Modules/Pathways for each original Boruta gene (KO counts)\n")
cat(strrep("=", 80), "\n")

for (ko in boruta_kos) {
  gene_lbl  <- get_label(ko)
  in_bins   <- ko %in% kos_in_any_bin
  ko_mods   <- seed_long %>% filter(KO == ko)

  if (!in_bins) {
    cat(sprintf("\n%s (%s): [EXCLUDED — not present in any genomic bin]\n",
                gene_lbl, ko))
    next
  }
  if (nrow(ko_mods) == 0) {
    cat(sprintf("\n%s (%s): [no modules/pathways after filtering]\n",
                gene_lbl, ko))
    next
  }

  cat(sprintf("\n%s (%s):\n", gene_lbl, ko))
  for (j in seq_len(nrow(ko_mods))) {
    pid         <- ko_mods$id[j]
    ptype       <- ko_mods$type[j]
    nm          <- name_cache[[pid]]; if (is.null(nm)) nm <- pid
    total_kos   <- length(pw_ko_cache[[pid]])
    present_kos <- length(intersect(pw_ko_cache[[pid]], rownames(merged_matrix)))
    in_heatmap  <- pid %in% group_stats$pw_id
    cat(sprintf("  [%s] %s (%s)\n    KEGG total KOs: %d | Present in data: %d | In heatmap: %s\n",
                pid, nm, ptype, total_kos, present_kos,
                if (in_heatmap) "YES" else "NO"))
  }
}

cat("\n")
cat("================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")
