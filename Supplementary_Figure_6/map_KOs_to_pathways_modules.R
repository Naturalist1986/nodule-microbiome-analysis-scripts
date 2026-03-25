#!/usr/bin/env Rscript
# ============================================================================
# MAP BORUTA KOs TO KEGG PATHWAYS AND MODULES
# ============================================================================
# Queries the KEGG REST API to find which pathways and modules each KO belongs
# to, then summarises which pathways/modules are enriched among the Boruta KOs.
# ============================================================================

library(tidyverse)
library(readxl)

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

cat("================================================================================\n")
cat("MAPPING BORUTA KOs TO KEGG PATHWAYS AND MODULES\n")
cat("================================================================================\n\n")

# ============================================================================
# CONFIG
# ============================================================================

SUMMARY_FILE   <- "Boruta_KO_heatmap_with_abundance_summary.csv"
GENE_NAME_FILE <- "ko_gene_names.csv"
CACHE_FILE     <- "ko_pathway_module_cache.rds"
KEGG_DELAY     <- 0.2   # seconds between API requests

# ============================================================================
# 1. LOAD KOs
# ============================================================================

cat("STEP 1: Loading Boruta KOs...\n")
summary_df <- read_csv(SUMMARY_FILE, show_col_types = FALSE)
gene_names <- read_csv(GENE_NAME_FILE, show_col_types = FALSE)

kos <- summary_df$KO
cat(sprintf("  %d KOs to map\n", length(kos)))

# ============================================================================
# 2. FETCH PATHWAY AND MODULE MEMBERSHIP FROM KEGG (with cache)
# ============================================================================

cat("\nSTEP 2: Fetching pathway/module membership from KEGG...\n")

fetch_ko_links <- function(ko_id) {
  result <- list(pathways = character(0), modules = character(0),
                 pathway_names = character(0), module_names = character(0))
  tryCatch({
    # --- Pathways ---
    url_pw <- paste0("https://rest.kegg.jp/link/pathway/", ko_id)
    con <- url(url_pw); on.exit(close(con), add = TRUE)
    txt <- readLines(con, warn = FALSE)
    # format: "ko:K00001\tpath:map00010"
    paths <- sub(".*path:", "", grep("path:map", txt, value = TRUE))
    result$pathways <- paths

    Sys.sleep(KEGG_DELAY)

    # --- Modules ---
    url_md <- paste0("https://rest.kegg.jp/link/module/", ko_id)
    con2 <- url(url_md); on.exit(close(con2), add = TRUE)
    txt2 <- readLines(con2, warn = FALSE)
    mods <- sub(".*md:", "", grep("md:", txt2, value = TRUE))
    result$modules <- mods

    Sys.sleep(KEGG_DELAY)
  }, error = function(e) {})
  result
}

fetch_name <- function(id, type = "pathway") {
  tryCatch({
    prefix <- if (type == "pathway") "map" else "md:"
    url_n <- paste0("https://rest.kegg.jp/list/", id)
    con <- url(url_n); on.exit(close(con), add = TRUE)
    txt <- readLines(con, warn = FALSE)
    if (length(txt) > 0) {
      parts <- strsplit(txt[1], "\t")[[1]]
      if (length(parts) >= 2) return(parts[2])
    }
    return(id)
  }, error = function(e) id)
}

# Use cached results if available
if (file.exists(CACHE_FILE)) {
  cat("  Loading cached results...\n")
  ko_links <- readRDS(CACHE_FILE)
} else {
  cat(sprintf("  Querying KEGG for %d KOs (pathways + modules)...\n", length(kos)))
  ko_links <- list()
  for (i in seq_along(kos)) {
    ko <- kos[i]
    cat(sprintf("  [%d/%d] %s\n", i, length(kos), ko))
    ko_links[[ko]] <- fetch_ko_links(ko)
  }
  saveRDS(ko_links, CACHE_FILE)
  cat("  Cached results saved.\n")
}

# ============================================================================
# 3. BUILD LONG-FORMAT TABLE: KO × PATHWAY and KO × MODULE
# ============================================================================

cat("\nSTEP 3: Building KO-pathway and KO-module tables...\n")

ko_pathway_long <- bind_rows(lapply(names(ko_links), function(ko) {
  paths <- ko_links[[ko]]$pathways
  if (length(paths) == 0) return(NULL)
  data.frame(KO = ko, pathway_id = paths, stringsAsFactors = FALSE)
}))

ko_module_long <- bind_rows(lapply(names(ko_links), function(ko) {
  mods <- ko_links[[ko]]$modules
  if (length(mods) == 0) return(NULL)
  data.frame(KO = ko, module_id = mods, stringsAsFactors = FALSE)
}))

cat(sprintf("  KO-pathway pairs: %d\n", nrow(ko_pathway_long)))
cat(sprintf("  KO-module pairs:  %d\n", nrow(ko_module_long)))

# ============================================================================
# 4. FETCH PATHWAY AND MODULE NAMES (with cache)
# ============================================================================

cat("\nSTEP 4: Fetching pathway and module names...\n")

NAME_CACHE_FILE <- "pathway_module_names_cache.rds"

if (file.exists(NAME_CACHE_FILE)) {
  name_cache <- readRDS(NAME_CACHE_FILE)
} else {
  name_cache <- list()
}

get_name_cached <- function(id) {
  if (!is.null(name_cache[[id]])) return(name_cache[[id]])
  nm <- fetch_name(id)
  name_cache[[id]] <<- nm
  Sys.sleep(KEGG_DELAY)
  nm
}

# Fetch pathway names
unique_paths <- unique(ko_pathway_long$pathway_id)
cat(sprintf("  Fetching %d pathway names...\n", length(unique_paths)))
pathway_names_df <- data.frame(
  pathway_id = unique_paths,
  pathway_name = sapply(unique_paths, get_name_cached),
  stringsAsFactors = FALSE
)

# Fetch module names
unique_mods <- unique(ko_module_long$module_id)
cat(sprintf("  Fetching %d module names...\n", length(unique_mods)))
module_names_df <- data.frame(
  module_id = unique_mods,
  module_name = sapply(unique_mods, get_name_cached),
  stringsAsFactors = FALSE
)

saveRDS(name_cache, NAME_CACHE_FILE)

# ============================================================================
# 5. JOIN NAMES AND GENE NAMES, SUMMARISE
# ============================================================================

cat("\nSTEP 5: Summarising pathways and modules...\n")

# Add gene names and Boruta metadata
summary_slim <- summary_df %>%
  dplyr::select(KO, Functional_Category, Traits, Mean_Importance) %>%
  left_join(gene_names, by = "KO") %>%
  mutate(Label = ifelse(!is.na(Gene_Name) & Gene_Name != "", Gene_Name, KO))

# --- Pathway summary ---
pathway_summary <- ko_pathway_long %>%
  left_join(pathway_names_df, by = "pathway_id") %>%
  left_join(summary_slim, by = "KO") %>%
  group_by(pathway_id, pathway_name) %>%
  summarise(
    n_KOs      = n_distinct(KO),
    KOs        = paste(sort(unique(KO)),   collapse = "; "),
    Gene_Names = paste(sort(unique(Label)), collapse = "; "),
    Traits     = paste(sort(unique(Traits)), collapse = "; "),
    Mean_Importance = round(mean(Mean_Importance, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_KOs))

# --- Module summary ---
module_summary <- ko_module_long %>%
  left_join(module_names_df, by = "module_id") %>%
  left_join(summary_slim, by = "KO") %>%
  group_by(module_id, module_name) %>%
  summarise(
    n_KOs      = n_distinct(KO),
    KOs        = paste(sort(unique(KO)),   collapse = "; "),
    Gene_Names = paste(sort(unique(Label)), collapse = "; "),
    Traits     = paste(sort(unique(Traits)), collapse = "; "),
    Mean_Importance = round(mean(Mean_Importance, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_KOs))

# --- KO-level table with all pathways and modules ---
ko_full <- summary_slim %>%
  left_join(
    ko_pathway_long %>%
      left_join(pathway_names_df, by = "pathway_id") %>%
      group_by(KO) %>%
      summarise(Pathways = paste(pathway_name, collapse = " | "), .groups = "drop"),
    by = "KO"
  ) %>%
  left_join(
    ko_module_long %>%
      left_join(module_names_df, by = "module_id") %>%
      group_by(KO) %>%
      summarise(Modules = paste(module_name, collapse = " | "), .groups = "drop"),
    by = "KO"
  ) %>%
  mutate(
    Pathways = replace_na(Pathways, "—"),
    Modules  = replace_na(Modules, "—")
  )

# ============================================================================
# 6. PRINT SUMMARY
# ============================================================================

cat("\n\nTOP PATHWAYS (by number of Boruta KOs):\n")
cat(strrep("=", 80), "\n")
for (i in 1:min(20, nrow(pathway_summary))) {
  row <- pathway_summary[i, ]
  cat(sprintf("[%d KOs] %s (%s)\n         Genes: %s\n",
              row$n_KOs, row$pathway_name, row$pathway_id, row$Gene_Names))
}

cat("\n\nTOP MODULES (by number of Boruta KOs):\n")
cat(strrep("=", 80), "\n")
for (i in 1:min(20, nrow(module_summary))) {
  row <- module_summary[i, ]
  cat(sprintf("[%d KOs] %s (%s)\n         Genes: %s\n",
              row$n_KOs, row$module_name, row$module_id, row$Gene_Names))
}

# ============================================================================
# 7. SAVE OUTPUTS
# ============================================================================

cat("\nSTEP 6: Saving outputs...\n")

write_csv(pathway_summary, "boruta_KOs_pathway_summary.csv")
write_csv(module_summary,  "boruta_KOs_module_summary.csv")
write_csv(ko_full,         "boruta_KOs_with_pathways_and_modules.csv")

cat("  Saved: boruta_KOs_pathway_summary.csv\n")
cat("  Saved: boruta_KOs_module_summary.csv\n")
cat("  Saved: boruta_KOs_with_pathways_and_modules.csv\n")

cat("\n================================================================================\n")
cat("DONE\n")
cat("================================================================================\n")
