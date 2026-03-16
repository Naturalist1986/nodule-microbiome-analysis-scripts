#!/usr/bin/env Rscript
# ============================================================================
# BORUTA KO HEATMAP VISUALIZATION WITH KEGG ANNOTATION + KO ABUNDANCES
# ============================================================================
# Creates ONE heatmap of all Boruta-selected KO features:
# - Symbiont vs Endophyte presence/absence heatmap with KO abundance annotation on top
#
# KO Abundance annotation shows raw abundance (hits/genome) with:
# - Blue-White-Red color scale
# - Scale saturates at 1 (values >1 shown as max red)
#
# Includes KEGG pathway annotations and functional categories
# ============================================================================

library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)
if (requireNamespace("svglite",  quietly = TRUE)) library(svglite)
if (requireNamespace("KEGGREST", quietly = TRUE)) library(KEGGREST)
library(readxl)

# Fix select() conflict with MASS package
select <- dplyr::select

cat("================================================================================\n")
cat("BORUTA KO HEATMAP GENERATION WITH KEGG ANNOTATION + KO ABUNDANCES\n")
cat("================================================================================\n\n")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Data directory (relative to this script)
DATA_DIR <- "../data"

# Input files
BORUTA_RESULTS_FILE <- file.path(DATA_DIR, "Boruta_KO_important_KOs.csv")
KO_FILE <- file.path(DATA_DIR, "normalized_kegg_results.tsv")
MEASUREMENTS_FILE <- file.path(DATA_DIR, "Guy_measurements.xlsx")
KO_BINS_FILE <- file.path(DATA_DIR, "ko_presence_absence_table.xlsx")
ABUNDANCE_FILE <- file.path(DATA_DIR, "combined_abundance_long_renormalized.tsv")

# Output prefix
OUTPUT_PREFIX <- "Boruta_KO_heatmap_with_abundance"

# MAG-KO attribution directory
ATTRIBUTION_DIR <- DATA_DIR

# KEGG API rate limiting (requests per second)
KEGG_DELAY <- 0.15  # 150ms between requests

# ============================================================================
# PHASE 1: LOAD BORUTA RESULTS
# ============================================================================

cat("STEP 1: Loading Boruta results...\n")

boruta_kos <- read_csv(BORUTA_RESULTS_FILE, show_col_types = FALSE)
unique_kos <- unique(boruta_kos$KO)

cat(sprintf("  Loaded %d important KOs across %d traits\n",
            length(unique_kos), length(unique(boruta_kos$Trait))))

# Load short gene names (ko_gene_names.csv, same file used by pathway_companion_heatmap.R)
GENE_NAMES_FILE <- "ko_gene_names.csv"
if (file.exists(GENE_NAMES_FILE)) {
  gene_names_df <- read_csv(GENE_NAMES_FILE, show_col_types = FALSE)
  ko_short_names <- setNames(
    ifelse(!is.na(gene_names_df$Gene_Name) & gene_names_df$Gene_Name != "",
           gene_names_df$Gene_Name, gene_names_df$KO),
    gene_names_df$KO
  )
  cat(sprintf("  Loaded %d short gene names from %s\n", nrow(gene_names_df), GENE_NAMES_FILE))
} else {
  ko_short_names <- setNames(character(0), character(0))
  cat("  ko_gene_names.csv not found — will fall back to KO IDs\n")
}

# ============================================================================
# PHASE 2: FETCH KEGG ANNOTATIONS
# ============================================================================

ANNOTATIONS_CACHE <- paste0(OUTPUT_PREFIX, "_annotations.csv")

cat("\nSTEP 2: Fetching KEGG annotations (this may take a while)...\n")

# Use cached annotations if available (avoids slow API calls)
if (file.exists(ANNOTATIONS_CACHE)) {
  cat(sprintf("  Found cached annotations: %s — skipping API calls.\n", ANNOTATIONS_CACHE))
  kegg_annotations <- read_csv(ANNOTATIONS_CACHE, show_col_types = FALSE)
  cat(sprintf("  Loaded %d annotations from cache.\n", nrow(kegg_annotations)))
} else {

fetch_kegg_annotation <- function(ko_id) {
  result <- tryCatch({
    kegg_info <- keggGet(ko_id)[[1]]

    # Get definition/name
    definition <- if (!is.null(kegg_info$DEFINITION)) {
      kegg_info$DEFINITION
    } else if (!is.null(kegg_info$NAME)) {
      paste(kegg_info$NAME, collapse = "; ")
    } else {
      NA_character_
    }

    # Get pathways
    pathways <- if (!is.null(kegg_info$PATHWAY)) {
      paste(names(kegg_info$PATHWAY), collapse = "; ")
    } else {
      NA_character_
    }

    # Get pathway descriptions for categorization
    pathway_desc <- if (!is.null(kegg_info$PATHWAY)) {
      paste(kegg_info$PATHWAY, collapse = "; ")
    } else {
      NA_character_
    }

    # Get BRITE hierarchy for functional category
    brite <- if (!is.null(kegg_info$BRITE)) {
      # Extract top-level category
      brite_text <- paste(names(kegg_info$BRITE), collapse = "; ")
      brite_text
    } else {
      NA_character_
    }

    return(list(
      definition = definition,
      pathways = pathways,
      pathway_desc = pathway_desc,
      brite = brite
    ))
  }, error = function(e) {
    return(list(
      definition = NA_character_,
      pathways = NA_character_,
      pathway_desc = NA_character_,
      brite = NA_character_
    ))
  })

  Sys.sleep(KEGG_DELAY)
  return(result)
}

# Fetch annotations for all unique KOs
cat(sprintf("  Fetching annotations for %d KOs...\n", length(unique_kos)))

annotations_list <- list()
pb <- txtProgressBar(min = 0, max = length(unique_kos), style = 3)

for (i in seq_along(unique_kos)) {
  ko <- unique_kos[i]
  annotations_list[[ko]] <- fetch_kegg_annotation(ko)
  setTxtProgressBar(pb, i)
}
close(pb)

# Convert to dataframe
kegg_annotations <- data.frame(
  KO = unique_kos,
  Definition = sapply(annotations_list, function(x) x$definition),
  Pathways = sapply(annotations_list, function(x) x$pathways),
  Pathway_Desc = sapply(annotations_list, function(x) x$pathway_desc),
  BRITE = sapply(annotations_list, function(x) x$brite),
  stringsAsFactors = FALSE
)

# ============================================================================
# PHASE 3: ASSIGN FUNCTIONAL CATEGORIES
# ============================================================================

cat("\nSTEP 3: Assigning functional categories...\n")

# Define functional category mapping based on pathway descriptions
assign_functional_category <- function(pathway_desc, brite, definition) {
  # Convert to lowercase for matching
  text <- tolower(paste(pathway_desc, brite, definition, sep = " "))

  # Check patterns and assign categories
  if (grepl("nitrogen|nitr|ammonia|glutam|urea|nitrate|nitrite|fixation", text)) {
    return("Nitrogen metabolism")
  } else if (grepl("carbon|glycol|glucon|pentose|citrate|tca|krebs|pyruvate", text)) {
    return("Carbon metabolism")
  } else if (grepl("amino acid|aminoacid|protein|peptide|biosynthesis of amino", text)) {
    return("Amino acid metabolism")
  } else if (grepl("oxidative|antioxid|glutathione|redox|peroxid|catalase|superoxide", text)) {
    return("Redox & stress response")
  } else if (grepl("sulfur|sulfate|cysteine|methionine|sulph", text)) {
    return("Sulfur metabolism")
  } else if (grepl("lipid|fatty|phospholipid|membrane", text)) {
    return("Lipid metabolism")
  } else if (grepl("vitamin|cofactor|coenzyme|riboflavin|folate|biotin|thiamin|pqq", text)) {
    return("Cofactor biosynthesis")
  } else if (grepl("transport|abc|permease|pump|export|import|secretion", text)) {
    return("Transport systems")
  } else if (grepl("signal|regul|kinase|phosphat|response|two-component", text)) {
    return("Signaling & regulation")
  } else if (grepl("cell wall|peptidoglycan|polysaccharide|glycan|lps", text)) {
    return("Cell envelope")
  } else if (grepl("dna|rna|replica|transcript|translat|ribosom", text)) {
    return("Genetic information")
  } else if (grepl("energy|atp|electron|respir|oxidase|reductase|nadh", text)) {
    return("Energy metabolism")
  } else if (grepl("secondary|polyketide|terpene|alkaloid|antibiotic", text)) {
    return("Secondary metabolism")
  } else {
    return("Other/Unknown")
  }
}

kegg_annotations$Functional_Category <- mapply(
  assign_functional_category,
  kegg_annotations$Pathway_Desc,
  kegg_annotations$BRITE,
  kegg_annotations$Definition
)

cat("  Functional categories assigned:\n")
print(table(kegg_annotations$Functional_Category))

# Save annotations
write.csv(kegg_annotations, ANNOTATIONS_CACHE, row.names = FALSE)

} # end else (KEGG API fetch)

# ============================================================================
# PHASE 4: LOAD AND PROCESS ABUNDANCE DATA
# ============================================================================

cat("\nSTEP 4: Loading and processing abundance data...\n")

# Load KO abundance data
ko_data <- read_tsv(KO_FILE, show_col_types = FALSE)

# Load measurements for sample matching
measurements <- read_excel(MEASUREMENTS_FILE)

# Extract sample info
sample_cols <- colnames(ko_data)[-1]

extract_sample_info <- function(sample_name) {
  clean_name <- gsub("_DKDN.*", "", sample_name)
  parts <- strsplit(clean_name, "_")[[1]]
  if (length(parts) >= 2) {
    treatment <- parts[1]
    replicate <- parts[2]
    return(c(treatment, replicate, paste(treatment, replicate, sep = "_")))
  }
  return(c(NA, NA, NA))
}

sample_info <- data.frame(
  original_name = sample_cols,
  t(sapply(sample_cols, extract_sample_info)),
  stringsAsFactors = FALSE
)
colnames(sample_info) <- c("original_name", "treatment", "replicate", "sample_id")

# Create KO matrix
ko_matrix <- as.matrix(ko_data[, -1])
rownames(ko_matrix) <- ko_data$kegg_number

# Merge technical replicates
unique_samples <- unique(sample_info$sample_id)
merged_matrix <- matrix(NA, nrow = nrow(ko_matrix), ncol = length(unique_samples))
rownames(merged_matrix) <- rownames(ko_matrix)
colnames(merged_matrix) <- unique_samples

for (sample_id in unique_samples) {
  cols <- sample_info$original_name[sample_info$sample_id == sample_id]
  if (length(cols) == 1) {
    merged_matrix[, sample_id] <- ko_matrix[, cols]
  } else {
    merged_matrix[, sample_id] <- rowMeans(ko_matrix[, cols], na.rm = TRUE)
  }
}

# Match with measurements
measurements$sample_id <- paste(measurements$Treatment, measurements$Replicate, sep = "_")
common_samples <- intersect(colnames(merged_matrix), measurements$sample_id)

ko_matrix_filtered <- merged_matrix[unique_kos[unique_kos %in% rownames(merged_matrix)], common_samples]
traits_filtered <- measurements %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

cat(sprintf("  Processed %d KOs x %d samples\n", nrow(ko_matrix_filtered), ncol(ko_matrix_filtered)))

# ============================================================================
# PHASE 5: CALCULATE TREATMENT MEANS
# ============================================================================

cat("\nSTEP 5: Calculating treatment-level statistics...\n")

# Treatment mapping
treatment_names <- c(
  "ces" = "CES",
  "RH" = "RH",
  "carR" = "CAR-G",
  "carK" = "CAR-C",
  "hok" = "HUK",
  "mtz" = "MTZ"
)
treatment_order <- c("CAR-C", "CAR-G", "CES", "HUK", "MTZ", "RH")

# Calculate mean abundance per treatment
treatment_means <- data.frame(KO = rownames(ko_matrix_filtered))

for (treat in names(treatment_names)) {
  treat_samples <- traits_filtered$sample_id[traits_filtered$Treatment == treat]
  if (length(treat_samples) > 0) {
    treat_cols <- intersect(treat_samples, colnames(ko_matrix_filtered))
    if (length(treat_cols) > 0) {
      treatment_means[[treatment_names[treat]]] <- rowMeans(ko_matrix_filtered[, treat_cols, drop = FALSE], na.rm = TRUE)
    }
  }
}

# ============================================================================
# PHASE 6: CALCULATE CORRELATIONS WITH TRAITS
# ============================================================================

cat("\nSTEP 6: Calculating correlations with traits...\n")

trait_cols <- c("percent_N", "Plant_Biomass", "Fixation_per_Nodule", "NMF")

correlations <- data.frame(KO = rownames(ko_matrix_filtered))

for (trait in trait_cols) {
  trait_vals <- traits_filtered[[trait]]
  cors <- apply(ko_matrix_filtered, 1, function(x) {
    valid <- !is.na(x) & !is.na(trait_vals)
    if (sum(valid) >= 3) cor(x[valid], trait_vals[valid]) else NA
  })
  pvals <- apply(ko_matrix_filtered, 1, function(x) {
    valid <- !is.na(x) & !is.na(trait_vals)
    if (sum(valid) >= 3) cor.test(x[valid], trait_vals[valid])$p.value else NA
  })
  correlations[[paste0("cor_", trait)]] <- cors
  correlations[[paste0("p_",   trait)]] <- pvals
}

# ============================================================================
# PHASE 6b: CALCULATE REGRESSION SLOPES (EFFECT SIZES) WITH TRAITS
# ============================================================================

cat("\nSTEP 6b: Calculating regression slopes (effect sizes)...\n")

slopes <- data.frame(KO = rownames(ko_matrix_filtered))

for (trait in trait_cols) {
  trait_vals <- traits_filtered[[trait]]
  sl <- apply(ko_matrix_filtered, 1, function(x) {
    valid <- !is.na(x) & !is.na(trait_vals)
    if (sum(valid) >= 3) coef(lm(trait_vals[valid] ~ x[valid]))[2] else NA
  })
  slopes[[paste0("slope_", trait)]] <- sl
}

# ============================================================================
# PHASE 7: MERGE ALL DATA FOR HEATMAP
# ============================================================================

cat("\nSTEP 7: Preparing heatmap data...\n")

# Merge Boruta results with annotations
boruta_with_annot <- boruta_kos %>%
  left_join(kegg_annotations, by = "KO") %>%
  left_join(correlations, by = "KO") %>%
  left_join(slopes, by = "KO")

# Create summary - one row per KO with all associated traits
ko_summary <- boruta_with_annot %>%
  group_by(KO, Definition, Functional_Category) %>%
  summarise(
    Traits = paste(unique(Trait), collapse = "; "),
    n_traits = n_distinct(Trait),
    Mean_Importance = mean(MeanImp, na.rm = TRUE),
    cor_percent_N = first(cor_percent_N),
    cor_Plant_Biomass = first(cor_Plant_Biomass),
    cor_Fixation_per_Nodule = first(cor_Fixation_per_Nodule),
    cor_NMF = first(cor_NMF),
    slope_percent_N = first(slope_percent_N),
    slope_Plant_Biomass = first(slope_Plant_Biomass),
    slope_Fixation_per_Nodule = first(slope_Fixation_per_Nodule),
    slope_NMF = first(slope_NMF),
    .groups = "drop"
  ) %>%
  arrange(Functional_Category, desc(n_traits), desc(Mean_Importance))

# Create short definition for labeling
ko_summary$Short_Def <- sapply(ko_summary$Definition, function(x) {
  if (is.na(x)) return("Unknown")
  # Truncate at first semicolon or 50 chars
  x_clean <- strsplit(x, ";")[[1]][1]
  if (nchar(x_clean) > 50) {
    paste0(substr(x_clean, 1, 47), "...")
  } else {
    x_clean
  }
})

# ============================================================================
# PHASE 8: BUILD ABUNDANCE HEATMAP MATRIX
# ============================================================================

cat("\nSTEP 8: Building abundance heatmap matrix...\n")

# Get KOs that are in both summary and treatment means
kos_for_heatmap <- intersect(ko_summary$KO, treatment_means$KO)
ko_summary <- ko_summary %>% filter(KO %in% kos_for_heatmap)

# Order KOs by functional category
ko_order <- ko_summary$KO

# Build abundance matrix (treatments x KOs)
abundance_mat <- treatment_means %>%
  filter(KO %in% ko_order) %>%
  column_to_rownames("KO")

# Reorder to match ko_order
abundance_mat <- abundance_mat[ko_order, treatment_order[treatment_order %in% colnames(abundance_mat)]]

# Transpose for heatmap (treatments as rows)
heatmap_mat <- t(as.matrix(abundance_mat))

# Z-score normalize across treatments for each KO
heatmap_mat_scaled <- scale(heatmap_mat)
# Handle any NaN from constant columns
heatmap_mat_scaled[is.nan(heatmap_mat_scaled)] <- 0

cat(sprintf("  Heatmap matrix: %d treatments x %d KOs\n", nrow(heatmap_mat_scaled), ncol(heatmap_mat_scaled)))

# ============================================================================
# PHASE 8b: LOAD MAG-KO ATTRIBUTION DATA
# ============================================================================

cat("\nSTEP 8b: Loading MAG-KO attribution data...\n")

ko_mag_matrix_full <- read_tsv(file.path(ATTRIBUTION_DIR, "ko_mag_sample_matrix.tsv"), show_col_types = FALSE)
mag_taxon_map      <- read_tsv(file.path(ATTRIBUTION_DIR, "ko_mag_covariation.tsv"),    show_col_types = FALSE) %>%
  select(KO, mag_id, taxon_group) %>%
  distinct()

# Extract treatment from replicate name (e.g. "carK_1" -> "carK")
ko_mag_matrix_full <- ko_mag_matrix_full %>%
  mutate(treatment_raw = sub("_[0-9]+$", "", replicate))

# Sum gene_read_equiv by KO x treatment x taxon_group
ko_treat_taxon <- ko_mag_matrix_full %>%
  left_join(mag_taxon_map, by = c("KO", "mag_id")) %>%
  group_by(KO, treatment_raw, taxon_group) %>%
  summarise(total_equiv = sum(gene_read_equiv, na.rm = TRUE), .groups = "drop")

ko_treat_total <- ko_treat_taxon %>%
  group_by(KO, treatment_raw) %>%
  summarise(total = sum(total_equiv), .groups = "drop")

# % Bradyrhizobium per KO x treatment (NA when KO has no MAG signal in that treatment)
ko_treat_brady <- ko_treat_taxon %>%
  filter(taxon_group == "Bradyrhizobium") %>%
  left_join(ko_treat_total, by = c("KO", "treatment_raw")) %>%
  mutate(pct_brady = ifelse(total > 0, total_equiv / total * 100, NA_real_)) %>%
  select(KO, treatment_raw, pct_brady)

# Build wide matrix: rows = treatment (display names), cols = KOs
attribution_wide <- expand_grid(
  KO            = unique_kos,
  treatment_raw = names(treatment_names)
) %>%
  left_join(ko_treat_brady, by = c("KO", "treatment_raw")) %>%
  mutate(treatment_full = treatment_names[treatment_raw])

attribution_mat_wide <- attribution_wide %>%
  select(KO, treatment_full, pct_brady) %>%
  pivot_wider(names_from = KO, values_from = pct_brady) %>%
  arrange(match(treatment_full, treatment_order)) %>%
  column_to_rownames("treatment_full") %>%
  as.matrix()

cat(sprintf("  Attribution matrix: %d treatments x %d KOs\n",
            nrow(attribution_mat_wide), ncol(attribution_mat_wide)))

# Z-score each KO (column) across treatments so treatment differences are amplified
# NA cells stay NA; KOs with zero variance become all-NA (no signal)
attribution_mat_wide_z <- apply(attribution_mat_wide, 2, function(x) {
  valid <- !is.na(x)
  if (sum(valid) < 2 || sd(x[valid]) == 0) return(rep(NA_real_, length(x)))
  out <- rep(NA_real_, length(x))
  out[valid] <- (x[valid] - mean(x[valid])) / sd(x[valid])
  out
})
rownames(attribution_mat_wide_z) <- rownames(attribution_mat_wide)

cat(sprintf("  Z-scored attribution matrix: range %.2f – %.2f\n",
            min(attribution_mat_wide_z, na.rm = TRUE),
            max(attribution_mat_wide_z, na.rm = TRUE)))

# ============================================================================
# HELPER FUNCTIONS: SIGNIFICANCE STARS + BARPLOT WITH ASTERISKS
# ============================================================================

# Convert p-value to significance stars
sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# Custom anno function: barplot + significance stars above/below each bar
make_sig_barplot_anno <- function(values, colors, stars, ylim = c(-1, 1)) {
  function(index, k = NULL, n = NULL) {
    n_bars <- length(index)
    vals   <- values[index]
    cols   <- colors[index]
    strs   <- stars[index]

    pushViewport(viewport(xscale = c(0.5, n_bars + 0.5), yscale = ylim))

    # Border and zero line
    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.5))
    grid.segments(0.5, 0, n_bars + 0.5, 0,
                  gp = gpar(col = "grey60", lwd = 0.5),
                  default.units = "native")

    for (i in seq_along(vals)) {
      v    <- if (is.na(vals[i])) 0 else max(ylim[1], min(ylim[2], vals[i]))
      fill <- if (is.na(cols[i])) "grey" else cols[i]

      if (v >= 0) {
        grid.rect(x = i, y = 0, width = 0.8, height = v,
                  just = c("centre", "bottom"),
                  gp = gpar(fill = fill, col = NA),
                  default.units = "native")
      } else {
        grid.rect(x = i, y = 0, width = 0.8, height = abs(v),
                  just = c("centre", "top"),
                  gp = gpar(fill = fill, col = NA),
                  default.units = "native")
      }

      # Always place asterisk at the very top of the annotation viewport
      if (nchar(strs[i]) > 0) {
        grid.text(strs[i],
                  x = unit(i, "native"),
                  y = unit(1, "npc") - unit(1, "mm"),
                  gp = gpar(fontsize = 8, fontface = "bold", col = "black"),
                  just = c("centre", "top"))
      }
    }

    # Y-axis
    grid.yaxis(gp = gpar(fontsize = 7))
    popViewport()
  }
}

# ============================================================================
# PHASE 9: CREATE ANNOTATIONS FOR ABUNDANCE HEATMAP
# ============================================================================

cat("\nSTEP 9: Creating annotations...\n")

# Column annotation (KO metadata)
col_annot_df <- ko_summary %>%
  filter(KO %in% colnames(heatmap_mat_scaled)) %>%
  arrange(match(KO, colnames(heatmap_mat_scaled)))

# Functional category colors
category_colors <- c(
  "Nitrogen metabolism" = "#1f77b4",
  "Carbon metabolism" = "#ff7f0e",
  "Amino acid metabolism" = "#2ca02c",
  "Redox & stress response" = "#d62728",
  "Sulfur metabolism" = "#9467bd",
  "Lipid metabolism" = "#8c564b",
  "Cofactor biosynthesis" = "#e377c2",
  "Transport systems" = "#7f7f7f",
  "Signaling & regulation" = "#bcbd22",
  "Cell envelope" = "#17becf",
  "Genetic information" = "#aec7e8",
  "Energy metabolism" = "#ffbb78",
  "Secondary metabolism" = "#98df8a",
  "Other/Unknown" = "#c7c7c7"
)

# Correlation barplot values - use the dominant trait's correlation
dominant_trait <- boruta_kos %>%
  group_by(KO) %>%
  slice_max(MeanImp, n = 1) %>%
  ungroup()

bar_values <- sapply(colnames(heatmap_mat_scaled), function(ko) {
  trait <- dominant_trait$Trait[dominant_trait$KO == ko][1]
  if (is.na(trait)) return(0)
  cor_col <- paste0("cor_", trait)
  val <- col_annot_df[[cor_col]][col_annot_df$KO == ko][1]
  if (is.na(val)) return(0)
  val
})

bar_colors <- ifelse(bar_values > 0, "#4CAF50", "#E91E63")
bar_colors[bar_values == 0] <- "grey"

# Compute p-values and significance stars for Phase 9 correlation bars
bar_pvalues <- sapply(colnames(heatmap_mat_scaled), function(ko) {
  trait <- dominant_trait$Trait[dominant_trait$KO == ko][1]
  if (is.na(trait)) return(NA_real_)
  p_col <- paste0("p_", trait)
  val <- col_annot_df[[p_col]][col_annot_df$KO == ko][1]
  if (length(val) == 0) NA_real_ else val
})
sig_labels_bar <- sig_stars(bar_pvalues)

ylim_bar <- c(-1, 1)

# ============================================================================
# PHASE 10: Z-SCORE HEATMAP REMOVED - SKIP TO PRESENCE/ABSENCE HEATMAP
# ============================================================================

cat("\nSTEP 10: Skipping Z-score abundance heatmap (removed)...\n")

# Color scale for MAG attribution:
#   0%  Brady (all NRE)   -> orange  (#d95f02)
#   50% Brady             -> white   (#f7f7f7)
#   100% Brady (symbiont) -> blue    (#1f78b4)
#   grey90 = no MAG signal (NA)
brady_colors <- colorRamp2(
  c(-2, 0, 2),
  c("#d95f02", "#f7f7f7", "#1f78b4")
)

# Legend for the attribution annotation
brady_legend <- Legend(
  col_fun    = brady_colors,
  title      = "% Bradyrhizobium\n(z-score per KO)",
  title_gp   = gpar(fontsize = 9, fontface = "bold"),
  labels_gp  = gpar(fontsize = 8),
  legend_height = unit(3, "cm"),
  direction  = "vertical"
)

# Save function for P/A heatmap
save_heatmap <- function(filename, dev_type, ht_obj, w, h) {
  if (dev_type == "svg" && !requireNamespace("svglite", quietly = TRUE)) {
    cat(sprintf("  Skipping %s (svglite not available)\n", filename))
    return(invisible(NULL))
  }

  cat(sprintf("  Saving %s...\n", filename))

  if (dev_type == "png") {
    png(filename, width = w, height = h, units = "in", res = 300)
  } else if (dev_type == "pdf") {
    pdf(filename, width = w, height = h)
  } else {
    svglite::svglite(filename, width = w, height = h)
  }

  draw(ht_obj,
       merge_legend = TRUE,
       padding = unit(c(15, 25, 5, 25), "mm"))

  dev.off()
}

# ============================================================================
# PHASE 11: CREATE SYMBIONT VS ENDOPHYTE PRESENCE/ABSENCE HEATMAP
# ============================================================================

cat("\nSTEP 11: Creating Symbiont vs Endophyte presence/absence heatmap...\n")

# Load KO presence in bins
cat("  Loading KO presence/absence data...\n")
ko_bins_raw <- read_excel(KO_BINS_FILE)

# Check data structure - Transpose if KOs are in rows
if ("KO" %in% colnames(ko_bins_raw)) {
  cat("  Transposing KO presence table...\n")
  ko_col <- ko_bins_raw$KO
  ko_bins_matrix <- as.matrix(ko_bins_raw[, -1])
  rownames(ko_bins_matrix) <- ko_col
  
  ko_bins_t <- t(ko_bins_matrix)
  ko_bins <- as.data.frame(ko_bins_t)
  ko_bins$bin <- rownames(ko_bins_t)
  rownames(ko_bins) <- NULL
  ko_bins <- ko_bins[, c("bin", setdiff(names(ko_bins), "bin"))]
} else {
  ko_bins <- as.data.frame(ko_bins_raw)
}

# Load bin classification from abundance file
cat("  Loading bin classification data...\n")
abundance_data <- read_tsv(ABUNDANCE_FILE, show_col_types = FALSE)

# Create bins dataframe
bins <- abundance_data %>%
  distinct(Bin_Prefixed, .keep_all = TRUE) %>%
  transmute(
    bin = Bin_Prefixed,
    treatment = Treatment,
    organism_type = tolower(Organism)
  )

# Get Boruta KOs present in bin data
boruta_kos_present <- intersect(unique_kos, colnames(ko_bins))

if (length(boruta_kos_present) > 0) {
  
  # Create long format
  ko_cols_to_use <- c("bin", boruta_kos_present)
  top_long <- ko_bins[, ko_cols_to_use] %>%
    pivot_longer(cols = -bin, names_to = "KO", values_to = "present") %>%
    filter(present == 1)
  
  # Add treatment and organism type
  top_long2 <- top_long %>%
    left_join(bins[, c("bin", "treatment", "organism_type")], by = "bin")
  
  # Summarize presence per KO x treatment x organism_type
  pres <- top_long2 %>%
    distinct(KO, treatment, organism_type, bin) %>%
    count(KO, treatment, organism_type, name = "n_bins")
  
  pres_wide <- pres %>%
    mutate(present_flag = 1L) %>%
    select(-n_bins) %>%
    pivot_wider(
      names_from = organism_type,
      values_from = present_flag,
      values_fill = 0L
    )
  
  # Create complete grid
  treat_short <- c("ces", "RH", "carR", "carK", "hok", "mtz")
  
  grid <- expand_grid(
    KO = boruta_kos_present,
    treatment = treat_short
  )
  
  if (!"endophyte" %in% colnames(pres_wide)) pres_wide$endophyte <- 0L
  if (!"symbiont" %in% colnames(pres_wide)) pres_wide$symbiont <- 0L
  
  pres_full <- grid %>%
    left_join(pres_wide, by = c("KO", "treatment")) %>%
    mutate(
      endophyte = replace_na(endophyte, 0L),
      symbiont = replace_na(symbiont, 0L),
      State = case_when(
        endophyte == 1L & symbiont == 1L ~ "Both",
        endophyte == 1L & symbiont == 0L ~ "Endophyte_only",
        endophyte == 0L & symbiont == 1L ~ "Symbiont_only",
        TRUE ~ "Absent"
      )
    )
  
  # Add metadata
  pres_full <- pres_full %>%
    mutate(Treatment_Full = treatment_names[treatment]) %>%
    left_join(kegg_annotations[, c("KO", "Functional_Category", "Definition")], by = "KO")
  
  # Add MeanImp (Importance) to aggregation
  ko_traits <- boruta_kos %>%
    group_by(KO) %>%
    slice_max(MeanImp, n = 1) %>%
    ungroup() %>%
    select(KO, Trait, MeanImp)
  
  pres_full <- pres_full %>%
    left_join(ko_traits, by = "KO")
  
  # Create trait labels
  trait_labels <- c(
    "Plant_Biomass" = "Plant\nbiomass (g)",
    "percent_N" = "Leaf N\n(%)",
    "NMF" = "Nodule mass\nfraction",
    "Fixation_per_Nodule" = "Fixation rate\nper nodule"
  )
  
  pres_full <- pres_full %>%
    mutate(
      Trait_label = trait_labels[as.character(Trait)],
      Trait_label = ifelse(is.na(Trait_label), "Multiple", Trait_label)
    )
  
  # Reorder KOs by TRAIT first, then IMPORTANCE
  ko_order_pa <- pres_full %>%
    distinct(KO, Trait, MeanImp, Functional_Category) %>%
    arrange(Trait, desc(MeanImp), Functional_Category) %>%
    pull(KO)
  
  pres_full <- pres_full %>%
    mutate(
      KO = factor(KO, levels = ko_order_pa),
      State = factor(State, levels = c("Absent", "Symbiont_only", "Endophyte_only", "Both")),
      Treatment_Full = factor(Treatment_Full, levels = treatment_order)
    )
  
  # State colors
  state_cols <- c("Absent" = "#f0f0f0", "Symbiont_only" = "#8DD3C7", 
                  "Endophyte_only" = "#FFFFB3", "Both" = "#BEBADA")
  
  # Build state matrix
  mat_state_df <- pres_full %>%
    select(Treatment_Full, KO, State) %>%
    arrange(Treatment_Full, KO) %>%
    pivot_wider(names_from = KO, values_from = State)
  
  state_mat <- as.matrix(mat_state_df[, -1])
  rownames(state_mat) <- mat_state_df$Treatment_Full
  
  # Create column annotations
  col_annot_pa <- pres_full %>%
    distinct(KO, Trait, Trait_label, Functional_Category, MeanImp) %>%
    arrange(match(KO, colnames(state_mat)))
  
  # ---------------------------------------------------------
  # TRACK 1: Boruta Importance Barplot
  # ---------------------------------------------------------
  bar_values_imp <- col_annot_pa$MeanImp
  bar_colors_imp <- "#7b3294" # Purple
  max_imp <- max(bar_values_imp, na.rm = TRUE)
  ylim_bar_imp <- c(0, max_imp * 1.1)
  
  # ---------------------------------------------------------
  # TRACK 2: Correlation Barplot (Restored)
  # ---------------------------------------------------------
  bar_values_corr <- sapply(colnames(state_mat), function(ko) {
    # Find the trait associated with this KO
    trait_raw <- col_annot_pa$Trait[col_annot_pa$KO == ko]
    if (is.na(trait_raw)) return(0)
    
    cor_col <- paste0("cor_", trait_raw)
    
    # Get the correlation value from the correlation dataframe calculated in Phase 6
    val <- correlations[[cor_col]][correlations$KO == ko]
    if (length(val) == 0 || is.na(val)) return(0)
    val
  })
  
  # Colors: Green for positive, Pink/Red for negative
  bar_colors_corr <- ifelse(bar_values_corr > 0, "#4CAF50", "#E91E63")
  bar_colors_corr[bar_values_corr == 0] <- "grey"

  # Compute p-values and significance stars for correlation bars
  bar_pvalues_corr <- sapply(colnames(state_mat), function(ko) {
    trait_raw <- col_annot_pa$Trait[col_annot_pa$KO == ko]
    if (length(trait_raw) == 0 || is.na(trait_raw)) return(NA_real_)
    p_col <- paste0("p_", trait_raw)
    val <- correlations[[p_col]][correlations$KO == ko]
    if (length(val) == 0 || is.na(val)) NA_real_ else val
  })
  sig_labels_corr <- sig_stars(bar_pvalues_corr)

  # Fix Y-limits for correlation (-1 to 1)
  ylim_bar_corr <- c(-1, 1)
  
  # ---------------------------------------------------------
  # ANNOTATION SETUP
  # ---------------------------------------------------------
  
  # Trait colors for annotation
  trait_colors <- c(
    "Plant\nbiomass (g)" = "#e41a1c",
    "Leaf N\n(%)" = "#377eb8",
    "Nodule mass\nfraction" = "#4daf4a",
    "Fixation rate\nper nodule" = "#984ea3",
    "Multiple" = "#ff7f00"
  )
  
  # MAG attribution sub-heatmap setup
  pa_ko_order <- colnames(state_mat)

  # Align attribution matrix rows/cols to the heatmap order
  attr_kos    <- pa_ko_order[pa_ko_order %in% colnames(attribution_mat_wide)]
  attr_treats <- treatment_order[treatment_order %in% rownames(attribution_mat_wide)]

  attribution_mat_pa <- matrix(
    NA_real_,
    nrow = length(attr_treats),
    ncol = length(pa_ko_order),
    dimnames = list(attr_treats, pa_ko_order)
  )
  attribution_mat_pa[attr_treats, attr_kos] <- attribution_mat_wide_z[attr_treats, attr_kos]

  ko_attribution_anno_pa <- function(index, k = NULL, n = NULL) {
    n = length(index)
    mat_subset <- attribution_mat_pa[, index, drop = FALSE]
    pushViewport(viewport(xscale = c(0.5, n + 0.5), yscale = c(0, 1)))
    n_treatments <- nrow(mat_subset)
    n_kos        <- ncol(mat_subset)
    for (i in seq_len(n_treatments)) {
      for (j in seq_len(n_kos)) {
        val      <- mat_subset[i, j]
        fill_col <- if (is.na(val)) "grey90" else brady_colors(val)
        grid.rect(
          x = j, y = (i - 0.5) / n_treatments,
          width = 1, height = 1 / n_treatments,
          gp = gpar(fill = fill_col, col = "white", lwd = 0.5),
          default.units = "native"
        )
      }
    }
    popViewport()
  }
  
  ha_top_pa <- HeatmapAnnotation(
    "Boruta\nImportance" = anno_barplot(
      bar_values_imp,
      gp = gpar(fill = bar_colors_imp, col = NA),
      height = unit(2.5, "cm"),
      ylim = ylim_bar_imp,
      axis_param = list(side = "left", gp = gpar(fontsize = 8))
    ),
    "Correlation\n(r)" = make_sig_barplot_anno(bar_values_corr, bar_colors_corr, sig_labels_corr, ylim_bar_corr),
    "Associated\nTrait" = col_annot_pa$Trait_label,
    "Functional\nCategory" = col_annot_pa$Functional_Category,
    col = list(
      "Associated\nTrait" = trait_colors,
      "Functional\nCategory" = category_colors
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    gap = unit(2, "mm"),
    simple_anno_size = unit(0.5, "cm"),
    annotation_height = unit(c(2.5, 2.5, 0.5, 0.5), "cm")
  )
  
  # Create Heatmap
  cell_side_pa <- unit(8, "mm")
  ht_width_pa <- cell_side_pa * ncol(state_mat)
  ht_height_pa <- cell_side_pa * nrow(state_mat)
  
  col_labels_pa <- sapply(colnames(state_mat), function(ko) {
    nm <- ko_short_names[ko]
    if (!is.na(nm) && nm != "") return(nm)
    # Fallback: use KO ID
    ko
  })
  
  ht_pa <- Heatmap(
    state_mat,
    name = "Gene status",
    col = state_cols,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = ht_width_pa,
    height = ht_height_pa,
    show_column_names = TRUE,
    column_labels = col_labels_pa,
    column_names_side = "bottom",
    column_names_max_height = unit(8, "cm"),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 60,
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    top_annotation = ha_top_pa,
    
    # Split by Trait
    column_split = col_annot_pa$Trait_label, 
    
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 1),
    row_names_side = "left",
    border = TRUE,
    column_gap = unit(3, "mm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  # Save P/A heatmap
  cat("\nSTEP 12: Saving presence/absence heatmap...\n")
  
  calc_w_pa <- convertWidth(ht_width_pa + unit(20, "cm"), "in", valueOnly = TRUE)
  calc_h_pa <- convertHeight(ht_height_pa + unit(22, "cm"), "in", valueOnly = TRUE)
  calc_w_pa <- max(calc_w_pa, 14)
  calc_h_pa <- max(calc_h_pa, 12)
  
  save_heatmap(paste0(OUTPUT_PREFIX, "_presence_absence.png"), "png", ht_pa, calc_w_pa, calc_h_pa)
  save_heatmap(paste0(OUTPUT_PREFIX, "_presence_absence.pdf"), "pdf", ht_pa, calc_w_pa, calc_h_pa)
  save_heatmap(paste0(OUTPUT_PREFIX, "_presence_absence.svg"), "svg", ht_pa, calc_w_pa, calc_h_pa)

  # -----------------------------------------------------------------------
  # EFFECT SIZE VERSION  (regression slope β instead of Pearson r)
  # -----------------------------------------------------------------------
  cat("\nSTEP 12b: Creating effect-size (slope) version of heatmap...\n")

  bar_values_slope <- sapply(colnames(state_mat), function(ko) {
    trait_raw <- col_annot_pa$Trait[col_annot_pa$KO == ko]
    if (length(trait_raw) == 0 || is.na(trait_raw)) return(0)
    slope_col <- paste0("slope_", trait_raw)
    val <- slopes[[slope_col]][slopes$KO == ko]
    if (length(val) == 0 || is.na(val)) return(0)
    val
  })

  bar_colors_slope <- ifelse(bar_values_slope > 0, "#4CAF50", "#E91E63")
  bar_colors_slope[bar_values_slope == 0] <- "grey"

  max_abs_slope <- max(abs(bar_values_slope[bar_values_slope != 0]), na.rm = TRUE)
  ylim_bar_slope <- c(-max_abs_slope * 1.15, max_abs_slope * 1.15)

  ha_top_slope <- HeatmapAnnotation(
    "Boruta\nImportance" = anno_barplot(
      bar_values_imp,
      gp = gpar(fill = bar_colors_imp, col = NA),
      height = unit(2.5, "cm"),
      ylim = ylim_bar_imp,
      axis_param = list(side = "left", gp = gpar(fontsize = 8))
    ),
    "Effect size\n(\u03b2)" = make_sig_barplot_anno(bar_values_slope, bar_colors_slope, sig_labels_corr, ylim_bar_slope),
    "Associated\nTrait" = col_annot_pa$Trait_label,
    "Functional\nCategory" = col_annot_pa$Functional_Category,
    col = list(
      "Associated\nTrait" = trait_colors,
      "Functional\nCategory" = category_colors
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    gap = unit(2, "mm"),
    simple_anno_size = unit(0.5, "cm"),
    annotation_height = unit(c(2.5, 2.5, 0.5, 0.5), "cm")
  )

  ht_slope <- Heatmap(
    state_mat,
    name = "Gene status",
    col = state_cols,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    width = ht_width_pa,
    height = ht_height_pa,
    show_column_names = TRUE,
    column_labels = col_labels_pa,
    column_names_side = "bottom",
    column_names_max_height = unit(8, "cm"),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 60,
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    top_annotation = ha_top_slope,
    column_split = col_annot_pa$Trait_label,
    column_title_gp = gpar(fontsize = 11, fontface = "bold"),
    rect_gp = gpar(col = "white", lwd = 1),
    row_names_side = "left",
    border = TRUE,
    column_gap = unit(3, "mm"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )

  save_heatmap(paste0(OUTPUT_PREFIX, "_presence_absence_effectsize.png"), "png", ht_slope, calc_w_pa, calc_h_pa)
  save_heatmap(paste0(OUTPUT_PREFIX, "_presence_absence_effectsize.pdf"), "pdf", ht_slope, calc_w_pa, calc_h_pa)
  save_heatmap(paste0(OUTPUT_PREFIX, "_presence_absence_effectsize.svg"), "svg", ht_slope, calc_w_pa, calc_h_pa)
  cat("  Saved effect-size version.\n")

} else {
  cat("  WARNING: No Boruta KOs found in bin presence data. Skipping P/A heatmap.\n")
}

# ============================================================================
# PHASE 13: EXPORT SUMMARY TABLES
# ============================================================================

cat("\nSTEP 13: Exporting summary tables...\n")

# Full summary with all info
full_summary <- ko_summary %>%
  left_join(kegg_annotations[, c("KO", "Pathways", "BRITE")], by = "KO") %>%
  arrange(Functional_Category, desc(n_traits), desc(Mean_Importance))

write.csv(full_summary, paste0(OUTPUT_PREFIX, "_summary.csv"), row.names = FALSE)

# Trait-specific summaries
for (trait in trait_cols) {
  trait_kos <- boruta_kos %>%
    filter(Trait == trait) %>%
    left_join(kegg_annotations, by = "KO") %>%
    left_join(correlations, by = "KO") %>%
    arrange(desc(MeanImp))

  write.csv(trait_kos, paste0(OUTPUT_PREFIX, "_", trait, ".csv"), row.names = FALSE)
}

# ============================================================================
# PHASE 14: LINEAR REGRESSION ANALYSIS
# ============================================================================

cat("\nSTEP 14: Running linear regressions for all KO-Trait pairs...\n")

# Build regression summary by iterating over all Boruta KO-Trait pairs
regression_results_list <- vector("list", nrow(boruta_kos))

for (i in 1:nrow(boruta_kos)) {
  ko <- boruta_kos$KO[i]
  trait <- boruta_kos$Trait[i]

  # Skip if KO not in filtered matrix
  if (!(ko %in% rownames(ko_matrix_filtered))) next

  # Get abundance values across samples
  abundance <- ko_matrix_filtered[ko, ]

  # Get trait values (aligned to common_samples via traits_filtered)
  trait_vals <- traits_filtered[[trait]]

  # Remove NA pairs
  valid <- !is.na(abundance) & !is.na(trait_vals)
  if (sum(valid) < 3) next

  x <- abundance[valid]
  y <- trait_vals[valid]

  # Fit linear model: trait ~ KO_abundance
  fit <- lm(y ~ x)
  fit_summary <- summary(fit)

  # Extract statistics
  r <- cor(x, y)
  r_sq <- fit_summary$r.squared
  adj_r_sq <- fit_summary$adj.r.squared
  p_val <- fit_summary$coefficients[2, 4]  # p-value for slope
  slope <- fit_summary$coefficients[2, 1]
  intercept <- fit_summary$coefficients[1, 1]
  n <- sum(valid)

  # Get annotation info
  def <- kegg_annotations$Definition[kegg_annotations$KO == ko]
  if (length(def) == 0) def <- NA_character_
  func_cat <- kegg_annotations$Functional_Category[kegg_annotations$KO == ko]
  if (length(func_cat) == 0) func_cat <- NA_character_

  regression_results_list[[i]] <- data.frame(
    KO = ko,
    Trait = trait,
    Definition = def,
    Functional_Category = func_cat,
    Correlation = r,
    R_squared = r_sq,
    Adj_R_squared = adj_r_sq,
    P_value = p_val,
    Slope = slope,
    Intercept = intercept,
    N_samples = n,
    stringsAsFactors = FALSE
  )
}

# Combine and sort
regression_results <- bind_rows(regression_results_list) %>%
  arrange(Trait, P_value)

# Export combined regression table
write.csv(regression_results,
          paste0(OUTPUT_PREFIX, "_regression_summary.csv"),
          row.names = FALSE)

cat(sprintf("  Computed %d regressions, saved to %s_regression_summary.csv\n",
            nrow(regression_results), OUTPUT_PREFIX))

# ============================================================================
# PHASE 15: SCATTER PLOT REGRESSIONS
# ============================================================================

cat("\nSTEP 15: Generating scatter plots for KO-Trait regressions...\n")

# Create output subdirectory
regression_plot_dir <- "regression_plots"
if (!dir.exists(regression_plot_dir)) {
  dir.create(regression_plot_dir)
}

# Trait display labels for axis labeling
trait_axis_labels <- c(
  "percent_N" = "Leaf N (%)",
  "Plant_Biomass" = "Plant Biomass (g)",
  "Fixation_per_Nodule" = "Fixation Rate per Nodule",
  "NMF" = "Nodule Mass Fraction"
)

# Treatment colors for scatter points
treatment_point_colors <- c(
  "CES" = "#e41a1c",
  "RH" = "#377eb8",
  "CAR-G" = "#4daf4a",
  "CAR-C" = "#984ea3",
  "HUK" = "#ff7f00",
  "MTZ" = "#a65628"
)

for (i in 1:nrow(regression_results)) {
  ko <- regression_results$KO[i]
  trait <- regression_results$Trait[i]

  if (!(ko %in% rownames(ko_matrix_filtered))) next

  # Get sample-level data
  abundance <- ko_matrix_filtered[ko, ]
  trait_vals <- traits_filtered[[trait]]
  treatments <- traits_filtered$Treatment

  valid <- !is.na(abundance) & !is.na(trait_vals)

  # Build plot dataframe
  plot_df <- data.frame(
    abundance = abundance[valid],
    trait_value = trait_vals[valid],
    treatment = treatments[valid],
    stringsAsFactors = FALSE
  ) %>%
    mutate(treatment_full = treatment_names[treatment])

  # Get statistics
  r_val <- regression_results$Correlation[i]
  r_sq <- regression_results$R_squared[i]
  p_val <- regression_results$P_value[i]

  # Format p-value for display
  p_label <- if (p_val < 0.001) {
    "p < 0.001"
  } else {
    sprintf("p = %.3f", p_val)
  }

  # Get short definition for title
  def <- regression_results$Definition[i]
  short_def <- if (is.na(def)) {
    ko
  } else {
    def_clean <- strsplit(def, ";")[[1]][1]
    if (nchar(def_clean) > 60) paste0(substr(def_clean, 1, 57), "...") else def_clean
  }

  # Build scatter plot with regression line
  p <- ggplot(plot_df, aes(x = abundance, y = trait_value)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                color = "grey30", fill = "grey80", alpha = 0.3) +
    geom_point(aes(color = treatment_full), size = 3, alpha = 0.8) +
    scale_color_manual(
      values = treatment_point_colors,
      name = "Treatment"
    ) +
    annotate("text",
             x = -Inf, y = Inf,
             label = sprintf("r = %.3f\nR\u00b2 = %.3f\n%s", r_val, r_sq, p_label),
             hjust = -0.1, vjust = 1.2,
             size = 3.5, fontface = "italic") +
    labs(
      title = sprintf("%s: %s", ko, short_def),
      subtitle = sprintf("Boruta-selected for: %s", trait_axis_labels[trait]),
      x = sprintf("%s abundance (hits/genome)", ko),
      y = trait_axis_labels[trait]
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "right"
    )

  # Save in both PNG and PDF formats
  filename_base <- sprintf("%s/%s_%s", regression_plot_dir, ko, trait)
  ggsave(paste0(filename_base, ".png"), plot = p,
         width = 7, height = 5, dpi = 300)
  ggsave(paste0(filename_base, ".pdf"), plot = p,
         width = 7, height = 5)
}

cat(sprintf("  Saved %d scatter plots to '%s/' directory\n",
            nrow(regression_results), regression_plot_dir))

cat("\n================================================================================\n")
cat("HEATMAP GENERATION COMPLETE\n")
cat("================================================================================\n")
cat(sprintf("Total KOs visualized: %d\n", ncol(heatmap_mat_scaled)))
cat(sprintf("Treatments: %d\n", nrow(heatmap_mat_scaled)))
cat(sprintf("Functional categories: %d\n", length(unique(col_annot_df$Functional_Category))))
cat("\nOutput files:\n")
cat(sprintf("  - %s_presence_absence.png/pdf/svg (symbiont/endophyte heatmap with KO abundance)\n", OUTPUT_PREFIX))
cat(sprintf("  - %s_summary.csv (full summary table)\n", OUTPUT_PREFIX))
cat(sprintf("  - %s_annotations.csv (KEGG annotations)\n", OUTPUT_PREFIX))
cat(sprintf("  - %s_<trait>.csv (trait-specific results)\n", OUTPUT_PREFIX))
cat(sprintf("  - %s_regression_summary.csv (all KO-trait regression statistics)\n", OUTPUT_PREFIX))
cat(sprintf("  - regression_plots/*.png/pdf (%d scatter plots with regression lines)\n", nrow(regression_results)))
cat("================================================================================\n")

