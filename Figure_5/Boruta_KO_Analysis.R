#!/usr/bin/env Rscript
# ============================================================================
# BORUTA FEATURE SELECTION FOR KEGG ORTHOLOG (KO) ANALYSIS
# ============================================================================
# Protocol based on metagenomic best practices:
# - Prevalence and abundance filtering
# - Zero handling with Bayesian-multiplicative replacement
# - Centered Log-Ratio (CLR) transformation
# - Boruta with ranger backend for computational efficiency
# ============================================================================

# Load required libraries
library(Boruta)
library(ranger)
library(readxl)
library(tidyverse)
library(zCompositions)  # For zero handling (cmultRepl)
library(compositions)   # For CLR transformation

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input/Output files
KO_FILE <- "normalized_kegg_results.tsv"
MEASUREMENTS_FILE <- "Guy_measurements.xlsx"
OUTPUT_PREFIX <- "Boruta_KO"

# Threading - set based on your system
N_THREADS <- 4  # Adjust for your Linux system

# Filtering thresholds
PREVALENCE_THRESHOLD <- 0.50    # Remove KOs present in < 10% of samples
MIN_ABUNDANCE_THRESHOLD <- 0.01 # Remove KOs with max proportion < 1%

# Boruta parameters
MAX_RUNS <- 500          # Increased for high-dimensional data
P_VALUE <- 0.01          # Strict significance threshold
NUM_TREES <- 500         # Trees per random forest iteration

# Random seed for reproducibility
SEED <- 42

cat("\n")
cat("========================================================================\n")
cat("BORUTA FEATURE SELECTION FOR KEGG ORTHOLOGS\n")
cat("========================================================================\n")
cat("Threads:", N_THREADS, "\n")
cat("Max runs:", MAX_RUNS, "\n")
cat("P-value threshold:", P_VALUE, "\n")
cat("Prevalence threshold:", PREVALENCE_THRESHOLD, "\n")
cat("========================================================================\n\n")

# ============================================================================
# PHASE 1: DATA LOADING
# ============================================================================

cat("=== PHASE 1: Loading Data ===\n")

# Load KO abundance data
ko_data <- read_tsv(KO_FILE, show_col_types = FALSE)
cat("  Loaded KO data:", nrow(ko_data), "KOs x", ncol(ko_data) - 1, "samples\n")

# Load measurements
measurements <- read_excel(MEASUREMENTS_FILE)
cat("  Loaded measurements:", nrow(measurements), "samples\n")
cat("  Traits:", paste(colnames(measurements)[!colnames(measurements) %in% c("Treatment", "Replicate")], collapse = ", "), "\n")

# ============================================================================
# PHASE 2: SAMPLE MATCHING AND MERGING TECHNICAL REPLICATES
# ============================================================================

cat("\n=== PHASE 2: Sample Matching ===\n")

sample_cols <- colnames(ko_data)[-1]

# Extract sample info from column names (Treatment_Replicate_DKDN...)
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
rownames(sample_info) <- NULL

# Create KO matrix
ko_matrix <- as.matrix(ko_data[, -1])
rownames(ko_matrix) <- ko_data$kegg_number

# Merge technical replicates by averaging
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

cat("  Merged to", ncol(merged_matrix), "unique samples\n")

# Match with measurements
measurements$sample_id <- paste(measurements$Treatment, measurements$Replicate, sep = "_")
common_samples <- intersect(colnames(merged_matrix), measurements$sample_id)

ko_matrix_filtered <- merged_matrix[, common_samples]
traits_filtered <- measurements %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

cat("  Common samples:", length(common_samples), "\n")

# ============================================================================
# PHASE 3: PREPROCESSING (CRUCIAL FOR METAGENOMICS)
# ============================================================================

cat("\n=== PHASE 3: Preprocessing ===\n")

# Transpose for sample x feature format
datExpr <- t(ko_matrix_filtered)

cat("  Initial dimensions:", nrow(datExpr), "samples x", ncol(datExpr), "KOs\n")

# --- 3.1 Prevalence Filtering ---
# Remove KOs present in fewer than threshold% of samples
zero_prop <- colSums(datExpr == 0) / nrow(datExpr)
keep_prevalence <- zero_prop < (1 - PREVALENCE_THRESHOLD)
cat("  After prevalence filter (>", PREVALENCE_THRESHOLD * 100, "% presence):", sum(keep_prevalence), "KOs\n")

# --- 3.2 Abundance Filtering ---
# Remove KOs with extremely low maximum relative abundance
max_abundance <- apply(datExpr, 2, max)
keep_abundance <- max_abundance >= MIN_ABUNDANCE_THRESHOLD
cat("  After abundance filter (max >=", MIN_ABUNDANCE_THRESHOLD, "):", sum(keep_abundance), "KOs\n")

# Apply both filters
keep_kos <- keep_prevalence & keep_abundance
datExpr_filtered <- datExpr[, keep_kos]
cat("  After combined filtering:", ncol(datExpr_filtered), "KOs retained\n")

# --- 3.3 Additional filtering for high-zero columns ---
# cmultRepl will delete columns with >80% zeros, so we pre-filter to avoid dimension mismatch
cat("  Pre-filtering columns with >80% zeros (required for cmultRepl)...\n")
zero_prop_strict <- colSums(datExpr_filtered == 0) / nrow(datExpr_filtered)
keep_for_imputation <- zero_prop_strict < 0.80
cat("    Removing", sum(!keep_for_imputation), "columns with >80% zeros\n")
datExpr_filtered <- datExpr_filtered[, keep_for_imputation]
cat("    After strict zero filter:", ncol(datExpr_filtered), "KOs retained\n")

# --- 3.4 Zero Handling (Bayesian-multiplicative replacement) ---
cat("  Handling zeros with Bayesian-multiplicative replacement...\n")

# Check for zeros
n_zeros <- sum(datExpr_filtered == 0)
if (n_zeros > 0) {
  cat("    Found", n_zeros, "zeros (", round(100 * n_zeros / length(datExpr_filtered), 2), "%)\n")

  # Apply cmultRepl for zero replacement (preserves covariance structure)
  # Using GBM method (Geometric Bayesian Multiplicative)
  # Set z.delete = FALSE to prevent automatic column deletion
  datExpr_imputed <- tryCatch({
    result <- cmultRepl(datExpr_filtered, method = "GBM", output = "p-counts",
                        suppress.print = TRUE, z.warning = 0.8, z.delete = FALSE)
    # Ensure we keep column names
    if (is.null(colnames(result))) {
      colnames(result) <- colnames(datExpr_filtered)
    }
    result
  }, error = function(e) {
    cat("    Warning: cmultRepl failed:", e$message, "\n")
    cat("    Using simple multiplicative replacement...\n")
    # Fallback: replace zeros with small value (half of minimum non-zero)
    temp <- datExpr_filtered
    min_nonzero <- min(temp[temp > 0], na.rm = TRUE)
    temp[temp == 0] <- min_nonzero / 2
    temp
  })
} else {
  datExpr_imputed <- datExpr_filtered
}

# Verify dimensions match
cat("    After imputation:", ncol(datExpr_imputed), "KOs\n")

# --- 3.5 CLR Transformation ---
cat("  Applying Centered Log-Ratio (CLR) transformation...\n")

# Store column names before CLR (in case clr() changes them)
ko_names_before_clr <- colnames(datExpr_imputed)

# CLR transformation
datExpr_clr <- as.data.frame(clr(datExpr_imputed))

# Restore proper column and row names
colnames(datExpr_clr) <- ko_names_before_clr
rownames(datExpr_clr) <- rownames(datExpr_filtered)

cat("  Final preprocessed data:", nrow(datExpr_clr), "samples x", ncol(datExpr_clr), "KOs\n")

# ============================================================================
# PHASE 4: BORUTA ANALYSIS FOR EACH TRAIT
# ============================================================================

cat("\n=== PHASE 4: Running Boruta Analysis ===\n")
cat("  This may take a while with", ncol(datExpr_clr), "features...\n\n")

# Define traits to analyze
trait_names <- c("percent_N", "Plant_Biomass", "Fixation_per_Nodule", "NMF")

# Store results
all_boruta_results <- list()
all_important_kos <- data.frame()

set.seed(SEED)

for (trait in trait_names) {
  cat("--- Analyzing trait:", trait, "---\n")

  # Get trait values
  y <- traits_filtered[[trait]]

  # Check for missing values
  if (any(is.na(y))) {
    cat("  Warning: Removing", sum(is.na(y)), "samples with missing values\n")
    valid_idx <- !is.na(y)
    y_clean <- y[valid_idx]
    x_clean <- datExpr_clr[valid_idx, ]
  } else {
    y_clean <- y
    x_clean <- datExpr_clr
  }

  cat("  Running Boruta with", nrow(x_clean), "samples,", ncol(x_clean), "features\n")

  # Run Boruta with ranger backend
  boruta_result <- Boruta(
    x = x_clean,
    y = y_clean,
    doTrace = 2,
    maxRuns = MAX_RUNS,
    pValue = P_VALUE,
    getImp = getImpRfZ,  # Use ranger for speed
    num.threads = N_THREADS,
    num.trees = NUM_TREES
  )

  # Store raw result
  all_boruta_results[[trait]] <- boruta_result

  # Check for tentative features
  n_confirmed <- sum(boruta_result$finalDecision == "Confirmed")
  n_rejected <- sum(boruta_result$finalDecision == "Rejected")
  n_tentative <- sum(boruta_result$finalDecision == "Tentative")

  cat("  Initial results: Confirmed=", n_confirmed, ", Rejected=", n_rejected,
      ", Tentative=", n_tentative, "\n")

  # Resolve tentative features if needed
  if (n_tentative > 0) {
    cat("  Resolving", n_tentative, "tentative features with TentativeRoughFix...\n")
    boruta_final <- TentativeRoughFix(boruta_result)
  } else {
    boruta_final <- boruta_result
  }

  # Extract confirmed KOs
  confirmed_kos <- getSelectedAttributes(boruta_final, withTentative = FALSE)

  cat("  Final confirmed features:", length(confirmed_kos), "\n")

  if (length(confirmed_kos) > 0) {
    # Get importance scores from ImpHistory
    # ImpHistory is a matrix with runs as rows and features as columns
    imp_history <- boruta_final$ImpHistory

    # Calculate mean importance for each confirmed KO
    mean_importances <- sapply(confirmed_kos, function(ko) {
      if (ko %in% colnames(imp_history)) {
        mean(imp_history[, ko], na.rm = TRUE)
      } else {
        NA_real_
      }
    })

    importance_df <- data.frame(
      KO = confirmed_kos,
      Trait = trait,
      MeanImp = mean_importances,
      Decision = "Confirmed",
      stringsAsFactors = FALSE
    )

    all_important_kos <- rbind(all_important_kos, importance_df)
  }

  # Save individual trait results
  saveRDS(boruta_final, paste0(OUTPUT_PREFIX, "_", trait, "_boruta_result.rds"))

  # Create importance plot
  pdf(paste0(OUTPUT_PREFIX, "_", trait, "_importance_plot.pdf"), width = 14, height = 8)
  plot(boruta_final, xlab = "", xaxt = "n", main = paste("Boruta Feature Importance -", trait))

  # Add labels for confirmed features only (to avoid clutter)
  if (n_confirmed > 0 && n_confirmed <= 50) {
    confirmed_idx <- which(boruta_final$finalDecision == "Confirmed")
    axis(1, at = confirmed_idx, labels = names(boruta_final$finalDecision)[confirmed_idx],
         las = 2, cex.axis = 0.6)
  }
  dev.off()

  cat("\n")
}

# ============================================================================
# PHASE 5: COMBINE AND EXPORT RESULTS
# ============================================================================

cat("=== PHASE 5: Exporting Results ===\n")

# Sort by importance within each trait
all_important_kos <- all_important_kos %>%
  arrange(Trait, desc(MeanImp))

# Get unique KOs across all traits
unique_kos <- unique(all_important_kos$KO)
cat("  Total unique important KOs:", length(unique_kos), "\n")

# Create summary by trait
ko_summary <- all_important_kos %>%
  group_by(Trait) %>%
  summarise(
    n_KOs = n(),
    KOs = paste(KO, collapse = ", "),
    .groups = "drop"
  )

cat("\n  KOs per trait:\n")
print(ko_summary[, c("Trait", "n_KOs")])

# Save results
write.csv(all_important_kos, paste0(OUTPUT_PREFIX, "_important_KOs.csv"), row.names = FALSE)
write.csv(data.frame(KO = unique_kos), paste0(OUTPUT_PREFIX, "_unique_KOs.csv"), row.names = FALSE)

# Save workspace
save(all_boruta_results, all_important_kos, unique_kos, datExpr_clr, traits_filtered,
     file = paste0(OUTPUT_PREFIX, "_workspace.RData"))

cat("\n")
cat("========================================================================\n")
cat("BORUTA ANALYSIS COMPLETE\n")
cat("========================================================================\n")
cat("Total unique important KOs:", length(unique_kos), "\n")
cat("\nOutput files:\n")
cat("  - ", OUTPUT_PREFIX, "_important_KOs.csv (all important KOs with traits)\n", sep = "")
cat("  - ", OUTPUT_PREFIX, "_unique_KOs.csv (unique KO list)\n", sep = "")
cat("  - ", OUTPUT_PREFIX, "_<trait>_boruta_result.rds (individual results)\n", sep = "")
cat("  - ", OUTPUT_PREFIX, "_<trait>_importance_plot.pdf (importance plots)\n", sep = "")
cat("  - ", OUTPUT_PREFIX, "_workspace.RData (full workspace)\n", sep = "")
cat("========================================================================\n")
