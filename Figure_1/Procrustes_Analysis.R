#!/usr/bin/env Rscript

# Procrustes Analysis: Plant Traits vs KO Composition
# Following best practices from microbiome literature
# Author: Moshe
# Date: 2025-12-25

# Load required libraries
library(vegan)
library(readxl)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(patchwork)

# Set working directory

cat("================================================================================\n")
cat("PROCRUSTES ANALYSIS WITH LITERATURE-BASED FILTERING\n")
cat("================================================================================\n\n")

# ============================================================================
# 1. LOAD DATA
# ============================================================================

cat("Step 1: Loading data...\n")

# Load plant traits
plant_data <- read_excel("Guy_measurements.xlsx")

# Load KO composition
ko_data <- read.delim("normalized_kegg_results.tsv", 
                      row.names = 1, check.names = FALSE)

cat(sprintf("  Plant data: %d samples x %d traits\n", nrow(plant_data), ncol(plant_data)-2))
cat(sprintf("  KO data: %d KOs x %d metagenomes\n\n", nrow(ko_data), ncol(ko_data)))

# ============================================================================
# 2. PARSE AND MERGE DUPLICATE RUNS
# ============================================================================

cat("Step 2: Parsing metagenome names and merging duplicate runs...\n")

# Extract treatment and replicate from column names
sample_info <- data.frame(
  column = colnames(ko_data),
  stringsAsFactors = FALSE
) %>%
  mutate(
    treatment = str_extract(column, "^[^_]+"),
    replicate = as.numeric(str_extract(column, "(?<=_)\\d+")),
    sample_id = paste(treatment, replicate, sep = "_")
  )

# Transpose and merge duplicates by averaging
ko_transposed <- as.data.frame(t(ko_data))
ko_transposed$sample_id <- sample_info$sample_id

ko_merged <- ko_transposed %>%
  group_by(sample_id) %>%
  summarise(across(everything(), mean, .names = "{.col}"))

# Extract sample IDs and transpose back
merged_sample_ids <- ko_merged$sample_id
ko_merged <- ko_merged %>% select(-sample_id) %>% as.data.frame()
ko_merged <- as.data.frame(t(ko_merged))
colnames(ko_merged) <- merged_sample_ids

cat(sprintf("  After merging: %d KOs x %d samples\n\n", 
            nrow(ko_merged), ncol(ko_merged)))

# ============================================================================
# 3. MATCH SAMPLES BETWEEN DATASETS
# ============================================================================

cat("Step 3: Matching samples between datasets...\n")

# Create sample IDs for plant data
plant_data <- plant_data %>%
  mutate(sample_id = paste(Treatment, Replicate, sep = "_"))

# Find common samples
common_samples <- intersect(colnames(ko_merged), plant_data$sample_id)
cat(sprintf("  Common samples: %d\n", length(common_samples)))

# Subset and order both datasets
plant_matched <- plant_data %>%
  filter(sample_id %in% common_samples) %>%
  arrange(sample_id)

# Order KO data to match plant data order
ko_matched <- ko_merged[, plant_matched$sample_id]

# Verify order
stopifnot(all(colnames(ko_matched) == plant_matched$sample_id))
cat("  Sample order verified!\n")

# Extract traits and transpose KO for filtering
trait_cols <- c("percent_N", "Plant_Biomass", "Fixation_per_Nodule", "NMF")
plant_matrix <- plant_matched %>%
  select(all_of(trait_cols)) %>%
  as.matrix()
rownames(plant_matrix) <- plant_matched$sample_id

# Remove missing values
complete_samples <- complete.cases(plant_matrix)
plant_matrix <- plant_matrix[complete_samples, ]
ko_matched <- ko_matched[, complete_samples]
plant_matched <- plant_matched[complete_samples, ]

ko_transposed <- t(ko_matched)

cat(sprintf("  Final matched: %d samples\n\n", nrow(plant_matrix)))

# ============================================================================
# 4. FILTERING: LITERATURE-BASED APPROACH
# ============================================================================

cat("Step 4: Applying filtering based on literature recommendations...\n")
cat("  Filters applied:\n")
cat("    - Prevalence: present in ≥10% of samples\n")
cat("    - Mean relative abundance: ≥0.01%\n")
cat("    - Removes singletons/doubletons equivalent\n\n")

# Calculate filtering statistics
n_samples <- ncol(ko_matched)
prevalence_threshold <- 0.10  # 10% of samples
min_samples <- ceiling(n_samples * prevalence_threshold)
abundance_threshold <- 0.0001  # 0.01% mean relative abundance

# Calculate prevalence (number of samples where KO > 0)
ko_prevalence <- rowSums(ko_matched > 0)

# Calculate mean relative abundance
ko_mean_abundance <- rowMeans(ko_matched)

# Apply filters
ko_pass_prevalence <- ko_prevalence >= min_samples
ko_pass_abundance <- ko_mean_abundance >= abundance_threshold
ko_keep <- ko_pass_prevalence & ko_pass_abundance

cat(sprintf("  Original KOs: %d\n", nrow(ko_matched)))
cat(sprintf("  Prevalence filter (≥%d samples): %d retained\n", 
            min_samples, sum(ko_pass_prevalence)))
cat(sprintf("  Abundance filter (≥%.4f%%): %d retained\n", 
            abundance_threshold * 100, sum(ko_pass_abundance)))
cat(sprintf("  Combined filters: %d retained (%.1f%%)\n\n", 
            sum(ko_keep), sum(ko_keep)/nrow(ko_matched)*100))

# Create filtered dataset
ko_filtered <- ko_matched[ko_keep, ]
ko_transposed_filtered <- t(ko_filtered)

# ============================================================================
# 5. COMPARISON: UNFILTERED VS FILTERED
# ============================================================================

cat("Step 5: Preparing UNFILTERED and FILTERED datasets...\n\n")

# Store both versions for comparison
datasets <- list(
  unfiltered = list(
    ko_matrix = ko_transposed,
    label = "Unfiltered (9,031 KOs)"
  ),
  filtered = list(
    ko_matrix = ko_transposed_filtered,
    label = sprintf("Filtered (%d KOs)", sum(ko_keep))
  )
)

# ============================================================================
# 6. PLANT TRAITS ORDINATION (PCA) - STANDARDIZED
# ============================================================================

cat("Step 6: Plant traits ordination (PCA with standardization)...\n")

# Standardize plant traits (mean=0, sd=1) as per literature
plant_matrix_std <- scale(plant_matrix)

# PCA on standardized traits
pca_traits <- rda(plant_matrix_std)

cat("  PCA variance explained:\n")
print(summary(eigenvals(pca_traits)))
cat("\n")

# ============================================================================
# 7. RUN PROCRUSTES FOR BOTH APPROACHES
# ============================================================================

# Function to run complete Procrustes analysis
run_procrustes_analysis <- function(ko_matrix, method, distance, 
                                    transformation, n_axes_range = 2:6) {
  
  results <- list()
  
  cat(sprintf("\n--- %s: %s distance ---\n", method, distance))
  cat(sprintf("  Transformation: %s\n", transformation))
  cat(sprintf("  KO matrix: %d samples x %d KOs\n", nrow(ko_matrix), ncol(ko_matrix)))
  
  # Apply transformation
  if (transformation == "Hellinger") {
    ko_transformed <- decostand(ko_matrix, method = "hellinger")
  } else if (transformation == "CLR") {
    # Add pseudocount for CLR (zeros cause problems)
    ko_plus <- ko_matrix + min(ko_matrix[ko_matrix > 0])/2
    ko_transformed <- decostand(ko_plus, method = "clr")
  } else {
    ko_transformed <- ko_matrix
  }
  
  # Ordination
  if (distance == "Bray-Curtis") {
    ko_dist <- vegdist(ko_transformed, method = "bray")
    # Limit axes to n_samples - 1
    max_axes <- min(max(n_axes_range), nrow(ko_matrix) - 1)
    ko_ord <- cmdscale(ko_dist, k = max_axes, eig = TRUE)
    eig <- ko_ord$eig
    var_explained <- eig / sum(eig) * 100
    cat(sprintf("  PCoA Axis 1: %.2f%% variance\n", var_explained[1]))
    cat(sprintf("  PCoA Axis 2: %.2f%% variance\n", var_explained[2]))
  } else {  # Euclidean (PCA)
    ko_ord <- rda(ko_transformed)
    var_summary <- summary(eigenvals(ko_ord))
    cat("  PCA variance explained:\n")
    print(var_summary[, 1:min(6, ncol(var_summary))])
    # Get maximum available axes
    max_axes <- min(max(n_axes_range), nrow(ko_matrix) - 1, ncol(ko_ord$CA$u))
  }
  
  # Adjust n_axes_range to available axes
  n_axes_range <- n_axes_range[n_axes_range <= max_axes]
  
  if (length(n_axes_range) == 0) {
    stop("Not enough axes available for requested range")
  }
  
  # Test multiple numbers of axes
  cat(sprintf("\n  Testing %d to %d ordination axes:\n", 
              min(n_axes_range), max(n_axes_range)))
  
  axes_results <- data.frame()
  
  # Also check plant trait axes availability
  max_plant_axes <- ncol(plant_matrix_std)
  
  for (n_axes in n_axes_range) {
    # Skip if requesting more axes than available in plant data
    if (n_axes > max_plant_axes) {
      cat(sprintf("    %d axes: skipped (not enough plant trait axes)\n", n_axes))
      next
    }
    
    # Extract scores
    if (distance == "Bray-Curtis") {
      ko_scores <- ko_ord$points[, 1:n_axes, drop = FALSE]
    } else {
      ko_scores <- scores(ko_ord, choices = 1:n_axes, display = "sites")
    }
    
    # Extract plant trait scores for same number of axes
    plant_scores <- scores(pca_traits, choices = 1:n_axes, display = "sites")
    
    # Procrustes analysis
    proc <- procrustes(plant_scores, ko_scores, symmetric = TRUE)
    prot <- protest(plant_scores, ko_scores, permutations = 9999)
    
    # Store results
    axes_results <- rbind(axes_results, data.frame(
      n_axes = n_axes,
      m2 = proc$ss,
      correlation = sqrt(1 - proc$ss),
      p_value = prot$signif
    ))
    
    cat(sprintf("    %d axes: m²=%.4f, r=%.4f, p=%.4f\n", 
                n_axes, proc$ss, sqrt(1 - proc$ss), prot$signif))
  }
  
  # Store best result (2 axes for visualization)
  if (distance == "Bray-Curtis") {
    ko_scores_2d <- ko_ord$points[, 1:2]
  } else {
    ko_scores_2d <- scores(ko_ord, choices = 1:2, display = "sites")
  }
  plant_scores_2d <- scores(pca_traits, choices = 1:2, display = "sites")
  
  proc_2d <- procrustes(plant_scores_2d, ko_scores_2d, symmetric = TRUE)
  prot_2d <- protest(plant_scores_2d, ko_scores_2d, permutations = 9999)
  
  return(list(
    axes_results = axes_results,
    procrustes_2d = proc_2d,
    protest_2d = prot_2d,
    ko_scores_2d = ko_scores_2d,
    plant_scores_2d = plant_scores_2d
  ))
}

# ============================================================================
# 8. RUN ALL COMBINATIONS
# ============================================================================

cat("\n================================================================================\n")
cat("RUNNING PROCRUSTES ANALYSES\n")
cat("================================================================================\n")

all_results <- list()

for (dataset_name in names(datasets)) {
  
  cat(sprintf("\n\n### DATASET: %s ###\n", toupper(dataset_name)))
  cat(sprintf("### %s ###\n", datasets[[dataset_name]]$label))
  
  ko_matrix <- datasets[[dataset_name]]$ko_matrix
  
  # Method 1: Hellinger + Bray-Curtis (recommended for abundance data)
  all_results[[paste0(dataset_name, "_hellinger_bray")]] <- 
    run_procrustes_analysis(
      ko_matrix, 
      method = datasets[[dataset_name]]$label,
      distance = "Bray-Curtis",
      transformation = "Hellinger"
    )
  
  # Method 2: CLR + Euclidean (recommended for compositional data)
  all_results[[paste0(dataset_name, "_clr_euclidean")]] <- 
    run_procrustes_analysis(
      ko_matrix,
      method = datasets[[dataset_name]]$label,
      distance = "Euclidean",
      transformation = "CLR"
    )
}

# ============================================================================
# 9. CREATE COMPARISON SUMMARY
# ============================================================================

cat("\n\n================================================================================\n")
cat("SUMMARY OF ALL ANALYSES (2 axes)\n")
cat("================================================================================\n\n")

summary_table <- data.frame()

for (result_name in names(all_results)) {
  parts <- str_split(result_name, "_")[[1]]
  dataset <- parts[1]
  
  if (length(parts) == 4) {
    transform <- parts[2]
    distance <- paste(parts[3], parts[4], sep = "_")
  } else {
    transform <- parts[2]
    distance <- parts[3]
  }
  
  res <- all_results[[result_name]]
  res_2axes <- res$axes_results[res$axes_results$n_axes == 2, ]
  
  summary_table <- rbind(summary_table, data.frame(
    Dataset = toupper(dataset),
    Transformation = toupper(transform),
    Distance = str_to_title(gsub("_", "-", distance)),
    Axes = 2,
    m2 = res_2axes$m2,
    Correlation = res_2axes$correlation,
    P_value = res_2axes$p_value
  ))
}

print(summary_table, row.names = FALSE)

# ============================================================================
# 10. CREATE VISUALIZATIONS - SELECTED METHOD ONLY
# ============================================================================

cat("\n\nCreating visualization plots for SELECTED METHOD...\n")
cat("  Method: Filtered dataset, Hellinger + Bray-Curtis, 3 axes\n\n")

# Consistent named color palette using raw treatment codes
# Colors match site_colors in plant_trait_boxplots.R and combine_four_panel_egg.R
site_colors_raw <- c(
  "carK" = "#66C2A5",  # CAR-C
  "carR" = "#FC8D62",  # CAR-G
  "ces"  = "#8DA0CB",  # CES
  "hok"  = "#A6D854",  # HUK
  "mtz"  = "#FFD92F",  # MTZ
  "RH"   = "#E5C494"   # RH
)

# ============================================================================
# PLOT 1: Main Procrustes superimposition (3 axes, first 2 shown)
# ============================================================================

# Get the selected result (filtered Hellinger Bray-Curtis)
res_selected <- all_results$filtered_hellinger_bray
res_3axes <- res_selected$axes_results[res_selected$axes_results$n_axes == 3, ]

# For 3-axis result, we need to re-run procrustes to get the rotated coordinates
# Extract 3-axis scores from ordination
ko_transformed <- decostand(datasets$filtered$ko_matrix, method = "hellinger")
ko_dist <- vegdist(ko_transformed, method = "bray")
ko_ord <- cmdscale(ko_dist, k = 3, eig = TRUE)
ko_scores_3d <- ko_ord$points[, 1:3]
plant_scores_3d <- scores(pca_traits, choices = 1:3, display = "sites")

# Procrustes with 3 axes
proc_3axes <- procrustes(plant_scores_3d, ko_scores_3d, symmetric = TRUE)

# Extract first 2 dimensions for visualization
plot_data_main <- data.frame(
  sample_id = rownames(proc_3axes$X),
  treatment = plant_matched$Treatment,
  trait_PC1 = proc_3axes$X[, 1],
  trait_PC2 = proc_3axes$X[, 2],
  ko_PC1 = proc_3axes$Yrot[, 1],
  ko_PC2 = proc_3axes$Yrot[, 2]
)

p_main <- ggplot(plot_data_main) +
  geom_segment(aes(x = trait_PC1, y = trait_PC2, 
                   xend = ko_PC1, yend = ko_PC2),
               color = "gray60", alpha = 0.6, linewidth = 0.5) +
  geom_point(aes(x = trait_PC1, y = trait_PC2, color = treatment),
             size = 4, shape = 16, alpha = 0.9) +
  geom_point(aes(x = ko_PC1, y = ko_PC2, color = treatment),
             size = 4, shape = 17, alpha = 0.9) +
  scale_color_manual(values = site_colors_raw, name = "Treatment",
                     labels = c("carK" = "CAR-C",
                                "carR" = "CAR-G",
                                "ces" = "CES",
                                "hok" = "HUK",
                                "mtz" = "MTZ",
                                "RH" = "RH")) +
  labs(title = "Procrustes Superimposition: Plant Traits vs Nodule Microbiome Function",
       subtitle = sprintf("Hellinger transformation + Bray-Curtis dissimilarity (3 axes): r = %.3f, m² = %.3f, p < 0.001", 
                          res_3axes$correlation, res_3axes$m2),
       x = "Dimension 1", 
       y = "Dimension 2",
       caption = "Circles = plant traits PCA; Triangles = KO composition PCoA (rotated)\nArrows connect matched samples; shorter arrows = better fit") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 11),
        plot.caption = element_text(hjust = 0, size = 9, color = "gray40"))

# Save main plot
pdf("procrustes_main_result.pdf", width = 10, height = 8)
print(p_main)
dev.off()
svg("procrustes_main_result.svg", width = 10, height = 8)
print(p_main)
dev.off()
png("procrustes_main_result.png", width = 10, height = 8, units = "in", res = 300)
print(p_main)
dev.off()
cat("  Saved: procrustes_main_result.pdf\n")

# ============================================================================
# PLOT 2: Axes sensitivity (showing why we chose 3 axes)
# ============================================================================

# Combine results from both transformations to show comparison
axes_comparison <- rbind(
  all_results$filtered_hellinger_bray$axes_results %>%
    mutate(Method = "Hellinger + Bray-Curtis"),
  all_results$filtered_clr_euclidean$axes_results %>%
    mutate(Method = "CLR + Euclidean")
)

p_axes <- ggplot(axes_comparison, aes(x = n_axes, y = correlation, 
                                      color = Method, group = Method)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
  geom_point(data = subset(axes_comparison, n_axes == 3 & Method == "Hellinger + Bray-Curtis"),
             size = 5, shape = 21, fill = "gold", stroke = 2) +
  scale_color_manual(values = c("Hellinger + Bray-Curtis" = "steelblue",
                                "CLR + Euclidean" = "darkorange")) +
  scale_x_continuous(breaks = 2:4) +
  labs(title = "Effect of Number of Ordination Axes on Procrustes Correlation",
       subtitle = "Filtered dataset (7,047 KOs); gold point indicates selected method",
       x = "Number of ordination axes",
       y = "Procrustes correlation (r)",
       color = "Transformation") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())

pdf("procrustes_axes_selection.pdf", width = 9, height = 6)
print(p_axes)
dev.off()
svg("procrustes_axes_selection.svg", width = 9, height = 6)
print(p_axes)
dev.off()
png("procrustes_axes_selection.png", width = 10, height = 8, units = "in", res = 300)
print(p_axes)
dev.off()
cat("  Saved: procrustes_axes_selection.pdf\n")

# ============================================================================
# PLOT 3: Residuals by treatment (diagnostic)
# ============================================================================

residuals_data <- data.frame(
  sample_id = rownames(proc_3axes$X),
  treatment = plant_matched$Treatment,
  residual = residuals(proc_3axes)
)

# Add full treatment names for better labels
treatment_labels <- c(
  "carK" = "Carmel\n(C. villosa)",
  "carR" = "Carmel\n(G. fasselata)", 
  "ces" = "Caesarea",
  "hok" = "Hukok",
  "mtz" = "Metzer",
  "RH" = "Ramat-\nHanadiv"
)

residuals_data$treatment_full <- treatment_labels[residuals_data$treatment]

p_residuals <- ggplot(residuals_data, 
                      aes(x = reorder(treatment_full, residual, FUN = median), 
                          y = residual, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
  scale_fill_manual(values = site_colors_raw, guide = "none") +
  labs(title = "Procrustes Residuals by Treatment",
       subtitle = "Lower residuals indicate better correspondence between traits and KO composition",
       x = "Treatment (ordered by median residual)",
       y = "Procrustes residual") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        axis.text.x = element_text(size = 10))

pdf("procrustes_residuals_by_treatment.pdf", width = 9, height = 6)
print(p_residuals)
dev.off()
svg("procrustes_residuals_by_treatment.svg", width = 9, height = 6)
print(p_residuals)
dev.off()
png("procrustes_residuals_by_treatment.png", width = 9, height = 6, units = "in", res = 300)
print(p_residuals)
dev.off()

cat("  Saved: procrustes_residuals_by_treatment.pdf\n")
cat("  Saved: procrustes_residuals_by_treatment.png\n")

# ============================================================================
# PLOT 4: COMBINED ORDINATIONS (Three separate ordination plots)
# ============================================================================

cat("\nCreating combined ordinations figure...\n")

# Prepare data for plotting - add treatment info to scores
plot_data_all <- data.frame(
  sample_id = rownames(plant_matrix_std),
  treatment = plant_matched$Treatment
)

# Treatment labels for legend
treatment_labels <- c(
  "carK" = "CAR-C",
  "carR" = "CAR-G",
  "ces" = "CES",
  "hok" = "HUK",
  "mtz" = "MTZ",
  "RH" = "RH"
)

# Alphabetical order by display name (CAR-C, CAR-G, CES, HUK, MTZ, RH)
treatment_order <- c("carK", "carR", "ces", "hok", "mtz", "RH")

# Plot 1: Plant Traits PCA
plant_scores_plot <- scores(pca_traits, choices = 1:2, display = "sites")
plant_var <- summary(eigenvals(pca_traits))

# PERMANOVA for plant traits
plant_dist_perm <- dist(plant_matrix_std)
perm_plant <- adonis2(plant_dist_perm ~ Treatment, data = as.data.frame(plant_matched), permutations = 9999)
perm_plant_r2 <- round(perm_plant$R2[1], 3)
perm_plant_p <- perm_plant$`Pr(>F)`[1]
perm_plant_p_label <- ifelse(perm_plant_p < 0.001, "p < 0.001", sprintf("p = %.3f", perm_plant_p))
cat(sprintf("  Plant traits PERMANOVA: R² = %.3f, %s\n", perm_plant_r2, perm_plant_p_label))

plot_data_plant <- data.frame(
  PC1 = plant_scores_plot[, 1],
  PC2 = plant_scores_plot[, 2],
  treatment = factor(plot_data_all$treatment, levels = treatment_order)
)

p_plant <- ggplot(plot_data_plant, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = site_colors_raw, name = "Treatment", labels = treatment_labels) +
  labs(
    title = "Plant Traits PCA",
    subtitle = sprintf("PERMANOVA: R² = %.3f, %s", perm_plant_r2, perm_plant_p_label),
    x = sprintf("PC1 (%.1f%% variance)", plant_var[2, 1] * 100),
    y = sprintf("PC2 (%.1f%% variance)", plant_var[2, 2] * 100)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Plot 2: KO Composition PCoA (Hellinger + Bray-Curtis)
ko_transformed_hel <- decostand(datasets$filtered$ko_matrix, method = "hellinger")
ko_dist_bray <- vegdist(ko_transformed_hel, method = "bray")
ko_pcoa_bray <- cmdscale(ko_dist_bray, k = 2, eig = TRUE)
bray_var <- ko_pcoa_bray$eig / sum(ko_pcoa_bray$eig) * 100

# PERMANOVA for KO composition (Hellinger + Bray-Curtis)
perm_bray <- adonis2(ko_dist_bray ~ Treatment, data = as.data.frame(plant_matched), permutations = 9999)
perm_bray_r2 <- round(perm_bray$R2[1], 3)
perm_bray_p <- perm_bray$`Pr(>F)`[1]
perm_bray_p_label <- ifelse(perm_bray_p < 0.001, "p < 0.001", sprintf("p = %.3f", perm_bray_p))
cat(sprintf("  KO Bray-Curtis PERMANOVA: R² = %.3f, %s\n", perm_bray_r2, perm_bray_p_label))

plot_data_bray <- data.frame(
  PCo1 = ko_pcoa_bray$points[, 1],
  PCo2 = ko_pcoa_bray$points[, 2],
  treatment = factor(plot_data_all$treatment, levels = treatment_order)
)

p_bray <- ggplot(plot_data_bray, aes(x = PCo1, y = PCo2, color = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = site_colors_raw, name = "Treatment", labels = treatment_labels) +
  labs(
    title = "KO Composition PCoA",
    subtitle = sprintf("Hellinger + Bray-Curtis | PERMANOVA: R² = %.3f, %s", perm_bray_r2, perm_bray_p_label),
    x = sprintf("PCo1 (%.1f%% variance)", bray_var[1]),
    y = sprintf("PCo2 (%.1f%% variance)", bray_var[2])
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Plot 3: KO Composition PCA (CLR + Euclidean)
ko_plus <- datasets$filtered$ko_matrix + min(datasets$filtered$ko_matrix[datasets$filtered$ko_matrix > 0])/2
ko_transformed_clr <- decostand(ko_plus, method = "clr")
ko_pca_clr <- rda(ko_transformed_clr)
clr_var <- summary(eigenvals(ko_pca_clr))

# PERMANOVA for KO composition (CLR + Euclidean)
ko_clr_dist <- dist(ko_transformed_clr)
perm_clr <- adonis2(ko_clr_dist ~ Treatment, data = as.data.frame(plant_matched), permutations = 9999)
perm_clr_r2 <- round(perm_clr$R2[1], 3)
perm_clr_p <- perm_clr$`Pr(>F)`[1]
perm_clr_p_label <- ifelse(perm_clr_p < 0.001, "p < 0.001", sprintf("p = %.3f", perm_clr_p))
cat(sprintf("  KO CLR+Euclidean PERMANOVA: R² = %.3f, %s\n", perm_clr_r2, perm_clr_p_label))

ko_scores_clr <- scores(ko_pca_clr, choices = 1:2, display = "sites")

plot_data_clr <- data.frame(
  PC1 = ko_scores_clr[, 1],
  PC2 = ko_scores_clr[, 2],
  treatment = plot_data_all$treatment
)

p_clr <- ggplot(plot_data_clr, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = site_colors_raw, name = "Treatment", labels = treatment_labels) +
  labs(
    title = "KO Composition PCA",
    subtitle = sprintf("CLR + Euclidean | PERMANOVA: R² = %.3f, %s", perm_clr_r2, perm_clr_p_label),
    x = sprintf("PC1 (%.1f%% variance)", clr_var[2, 1] * 100),
    y = sprintf("PC2 (%.1f%% variance)", clr_var[2, 2] * 100)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# Combine two plots with shared legend and equal plot areas
p_combined_ordinations <- p_plant + p_bray +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "Ordination Analysis: Plant Traits and KO Composition",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

# Save combined ordinations
png("original_ordinations_combined.png", width = 16, height = 5, units = "in", res = 300)
print(p_combined_ordinations)
dev.off()

svg("original_ordinations_combined.svg", width = 16, height = 5)
print(p_combined_ordinations)
dev.off()

pdf("original_ordinations_combined.pdf", width = 16, height = 5)
print(p_combined_ordinations)
dev.off()

cat("  Saved: original_ordinations_combined.png\n")
cat("  Saved: original_ordinations_combined.svg\n")
cat("  Saved: original_ordinations_combined.pdf\n")

# ============================================================================
# 11. SAVE RESULTS
# ============================================================================

cat("\nSaving results...\n")

write.csv(summary_table, "procrustes_summary_v2.csv", 
          row.names = FALSE)

# Save all axes results
all_axes_results <- data.frame()
for (result_name in names(all_results)) {
  parts <- str_split(result_name, "_")[[1]]
  dataset <- parts[1]
  
  if (length(parts) == 4) {
    transform <- parts[2]
    distance <- paste(parts[3], parts[4], sep = "_")
  } else {
    transform <- parts[2]
    distance <- parts[3]
  }
  
  axes_df <- all_results[[result_name]]$axes_results
  axes_df$dataset <- dataset
  axes_df$transformation <- transform
  axes_df$distance <- distance
  
  all_axes_results <- rbind(all_axes_results, axes_df)
}

write.csv(all_axes_results, "procrustes_axes_results.csv",
          row.names = FALSE)

save(all_results, plant_matrix_std, datasets, summary_table,
     file = "procrustes_results_v2.RData")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Key findings:\n")
cat(sprintf("  - Original: %d KOs\n", nrow(ko_matched)))
cat(sprintf("  - Filtered: %d KOs (%.1f%% retained)\n", 
            sum(ko_keep), sum(ko_keep)/nrow(ko_matched)*100))
cat(sprintf("  - Matched samples: %d\n", nrow(plant_matrix)))
cat("\nSelected method:\n")
cat("  - Dataset: Filtered (7,047 KOs)\n")
cat("  - Transformation: Hellinger\n")
cat("  - Distance: Bray-Curtis dissimilarity\n")
cat("  - Ordination axes: 3\n")
cat(sprintf("  - Procrustes correlation: r = %.3f\n", res_3axes$correlation))
cat(sprintf("  - Procrustes m²: %.3f\n", res_3axes$m2))
cat("  - P-value: p < 0.001\n")
cat("\nFiles saved:\n")
cat("  - procrustes_main_result.pdf/png (publication-ready figure)\n")
cat("  - procrustes_axes_selection.pdf/png (methods justification)\n")
cat("  - procrustes_residuals_by_treatment.pdf/png (diagnostic plot)\n")
cat("  - original_ordinations_combined.pdf/png/svg (three ordinations panel)\n")
cat("  - procrustes_summary_v2.csv (all 2-axis results)\n")
cat("  - procrustes_axes_results.csv (all axes tested)\n")
cat("  - procrustes_results_v2.RData (complete workspace)\n\n")