#!/usr/bin/env Rscript

# Mantel Test Analysis with Integrated Figure using linkET
# Creates correlation heatmap + Mantel test network visualization
# Author: Moshe
# Date: 2025-01-06

# Load required libraries
library(vegan)
library(readxl)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

# Install linkET if not already installed
if (!require(linkET)) {
  if (!require(devtools)) install.packages("devtools")
  devtools::install_github("Hy4m/linkET")
  library(linkET)
}

cat("================================================================================\n")
cat("MANTEL TEST ANALYSIS WITH INTEGRATED FIGURE (linkET)\n")
cat("================================================================================\n\n")

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================

cat("Step 1: Loading data...\n")

plant_data <- read_excel("Guy_measurements_old.xlsx")
ko_data <- read.delim("normalized_kegg_results.tsv", 
                      row.names = 1, check.names = FALSE)

cat(sprintf("  Plant data: %d samples x %d traits\n", nrow(plant_data), ncol(plant_data)-2))
cat(sprintf("  KO data: %d KOs x %d metagenomes\n\n", nrow(ko_data), ncol(ko_data)))

# ============================================================================
# 2. PARSE AND MERGE DUPLICATE RUNS
# ============================================================================

cat("Step 2: Parsing and merging duplicate runs...\n")

sample_info <- data.frame(
  column = colnames(ko_data),
  stringsAsFactors = FALSE
) %>%
  mutate(
    treatment = str_extract(column, "^[^_]+"),
    replicate = as.numeric(str_extract(column, "(?<=_)\\d+")),
    sample_id = paste(treatment, replicate, sep = "_")
  )

ko_transposed <- as.data.frame(t(ko_data))
ko_transposed$sample_id <- sample_info$sample_id

ko_merged <- ko_transposed %>%
  group_by(sample_id) %>%
  summarise(across(everything(), mean, .names = "{.col}"))

merged_sample_ids <- ko_merged$sample_id
ko_merged <- ko_merged %>% select(-sample_id) %>% as.data.frame()
ko_merged <- as.data.frame(t(ko_merged))
colnames(ko_merged) <- merged_sample_ids

cat(sprintf("  After merging: %d KOs x %d samples\n\n", nrow(ko_merged), ncol(ko_merged)))

# ============================================================================
# 3. MATCH SAMPLES
# ============================================================================

cat("Step 3: Matching samples...\n")

plant_data <- plant_data %>%
  mutate(sample_id = paste(Treatment, Replicate, sep = "_"))

common_samples <- intersect(colnames(ko_merged), plant_data$sample_id)
cat(sprintf("  Common samples: %d\n", length(common_samples)))

plant_matched <- plant_data %>%
  filter(sample_id %in% common_samples) %>%
  arrange(sample_id)

ko_matched <- ko_merged[, plant_matched$sample_id]

# Select ALL trait columns (exclude Replicate and Treatment)
trait_cols <- c("percent_N", "Height", "Nodule_Biomass", "Plant_Biomass", 
                "Fixation_per_Nodule", "Fixation_per_Plant", "RMF", "NMF")

plant_matrix <- plant_matched %>%
  select(all_of(trait_cols)) %>%
  as.matrix()
rownames(plant_matrix) <- plant_matched$sample_id

complete_samples <- complete.cases(plant_matrix)
plant_matrix <- plant_matrix[complete_samples, ]
ko_matched <- ko_matched[, complete_samples]

ko_transposed <- t(ko_matched)

cat(sprintf("  Final matched: %d samples\n\n", nrow(plant_matrix)))

# ============================================================================
# 4. FILTERING
# ============================================================================

cat("Step 4: Filtering KOs...\n")

n_samples <- ncol(ko_matched)
prevalence_threshold <- 0.10
min_samples <- ceiling(n_samples * prevalence_threshold)
abundance_threshold <- 0.0001

ko_prevalence <- rowSums(ko_matched > 0)
ko_mean_abundance <- rowMeans(ko_matched)
ko_keep <- (ko_prevalence >= min_samples) & (ko_mean_abundance >= abundance_threshold)

ko_filtered <- ko_matched[ko_keep, ]
ko_transposed_filtered <- t(ko_filtered)

cat(sprintf("  Retained: %d KOs (%.1f%%)\n\n", sum(ko_keep), sum(ko_keep)/nrow(ko_matched)*100))

# ============================================================================
# 5. PREPARE DATA FOR linkET
# ============================================================================

cat("Step 5: Preparing data for linkET...\n")

# Clean trait names for display (all 8 traits)
trait_label_map <- c(
  "percent_N" = "Leaf N (%)",
  "Height" = "Plant Height (cm)",
  "Nodule_Biomass" = "Nodule Biomass (g)",
  "Plant_Biomass" = "Plant Biomass (g)",
  "Fixation_per_Nodule" = "N₂ Fixation per Nodule",
  "Fixation_per_Plant" = "N₂ Fixation per Plant",
  "RMF" = "Root Mass Fraction",
  "NMF" = "Nodule Mass Fraction"
)

# Create trait dataframe with clean names
trait_df <- as.data.frame(plant_matrix)
colnames(trait_df) <- trait_label_map[colnames(trait_df)]

# For linkET, we need the KO composition as a "species" matrix
# and traits as "environmental" variables
species_df <- as.data.frame(ko_transposed_filtered)

cat(sprintf("  Trait data: %d samples x %d traits\n", nrow(trait_df), ncol(trait_df)))
cat(sprintf("  KO data: %d samples x %d KOs\n\n", nrow(species_df), ncol(species_df)))

# ============================================================================
# 6. CALCULATE CORRELATIONS AND MANTEL TESTS
# ============================================================================

cat("Step 6: Calculating correlations and Mantel tests...\n")

# A. Pearson correlations between traits (for heatmap)
cor_res <- correlate(trait_df, method = "pearson")

# B. Mantel tests between traits and KO composition
# linkET's mantel_test expects:
#   - spec: the "species" data (KO composition)
#   - env: the "environmental" data (traits)
mantel_res <- mantel_test(species_df, trait_df,
                          spec_select = list(`Functional composition` = 1:ncol(species_df)),
                          mantel_fun = "mantel") %>%
  # Categorize p-values for coloring
  mutate(
    p_cat = cut(p, 
                breaks = c(-Inf, 0.01, 0.05, Inf),
                labels = c("< 0.01", "0.01 - 0.05", "≥ 0.05"),
                right = FALSE)
  ) %>%
  # Categorize |r| for line width
  mutate(
    r_cat = cut(abs(r), 
                breaks = c(-Inf, 0.1, 0.2, Inf),
                labels = c("< 0.1", "0.1 - 0.2", "≥ 0.2"),
                right = FALSE)
  ) %>%
  # Direction for line type
  mutate(
    sign_cat = ifelse(r > 0, "Positive", "Negative")
  )

cat("  Correlation matrix calculated\n")
cat("  Mantel tests completed\n\n")

# Print Mantel results
cat("=== MANTEL TEST RESULTS ===\n")
for (i in 1:nrow(mantel_res)) {
  sig_marker <- ifelse(mantel_res$p[i] < 0.01, "**",
                       ifelse(mantel_res$p[i] < 0.05, "*", ""))
  cat(sprintf("  %-30s r = %6.4f, p = %.4f %s\n",
              mantel_res$env[i],
              mantel_res$r[i],
              mantel_res$p[i],
              sig_marker))
}
cat("\n** p < 0.01, * p < 0.05\n\n")

# ============================================================================
# 7. CREATE INTEGRATED FIGURE WITH linkET
# ============================================================================

cat("Step 7: Creating integrated figure...\n")

# Define color palettes
heatmap_palette <- rev(brewer.pal(n = 11, name = "RdBu"))
mantel_p_colors <- c(
  "< 0.01" = "#D95F02",      # Orange
  "0.01 - 0.05" = "#7570B3", # Purple
  "≥ 0.05" = "grey70"        # Grey
)

# Create the integrated plot
p <- qcorrplot(cor_res, type = "lower", diag = FALSE) +
  
  # Add correlation heatmap
  geom_square() +
  scale_fill_gradientn(
    colours = heatmap_palette,
    limits = c(-1, 1),
    name = "Pearson r",
    guide = guide_colorbar(order = 1)
  ) +
  
  # Add Mantel test connections
  geom_couple(
    data = mantel_res,
    aes(colour = p_cat, size = r_cat, linetype = sign_cat),
    curvature = nice_curvature(),
    nudge_x = 0.1,
    label.size = 5
  ) +
  
  # Customize Mantel test aesthetics
  scale_colour_manual(
    values = mantel_p_colors,
    name = "Mantel p-value",
    guide = guide_legend(
      order = 2,
      override.aes = list(size = 2, linewidth = 1)
    )
  ) +
  
  scale_size_manual(
    values = c("< 0.1" = 0.5, "0.1 - 0.2" = 1.5, "≥ 0.2" = 2.5),
    name = "Mantel |r|",
    guide = guide_legend(order = 3)
  ) +
  
  scale_linetype_manual(
    values = c("Positive" = "solid", "Negative" = "dashed"),
    name = "Direction",
    guide = guide_legend(
      order = 4,
      override.aes = list(size = 1.5, linewidth = 1)
    )
  ) +
  
  # Theme adjustments (optimized for 8x8 matrix)
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.text.y.right = element_text(size = 24),
    legend.position = "right",
    legend.key = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 14),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  
  # Add title
  labs(
    title = "Plant Trait Correlations and Associations with Nodule Microbiome Function",
    subtitle = "Left: Pearson correlations between traits | Right: Mantel tests with KO composition"
  )

# Display
print(p)

cat("  Figure created successfully\n\n")

# ============================================================================
# 8. SAVE FIGURES
# ============================================================================

cat("Step 8: Saving figures...\n")

# Save in multiple formats (larger size for 8x8 matrix)
ggsave("mantel_linkET_figure.pdf", plot = p, width = 14, height = 10, dpi = 300)
ggsave("mantel_linkET_figure.png", plot = p, width = 14, height = 10, dpi = 300)
ggsave("mantel_linkET_figure.svg", plot = p, width = 14, height = 10)

cat("  Saved: mantel_linkET_figure (pdf/png/svg)\n\n")

# ============================================================================
# OPTIONAL: Create filtered version (excluding math-dependent correlations)
# ============================================================================

cat("Creating filtered version (optional)...\n")

# Define mathematically dependent pairs
dependent_pairs <- list(
  c("N₂ Fixation per Plant", "Nodule Biomass (g)"),
  c("N₂ Fixation per Plant", "N₂ Fixation per Nodule"),
  c("Nodule Mass Fraction", "Nodule Biomass (g)"),
  c("Nodule Mass Fraction", "Plant Biomass (g)"),
  c("Root Mass Fraction", "Plant Biomass (g)")
)

# Function to filter correlation matrix
filter_cor_matrix <- function(cor_res, dependent_pairs) {
  # Create a copy
  filtered <- cor_res
  
  # Set dependent pairs to NA
  for (pair in dependent_pairs) {
    var1 <- pair[1]
    var2 <- pair[2]
    
    # Find indices
    idx1 <- which(filtered$var1 == var1 & filtered$var2 == var2)
    idx2 <- which(filtered$var1 == var2 & filtered$var2 == var1)
    
    # Set to NA (will not be plotted)
    if (length(idx1) > 0) filtered$r[idx1] <- NA
    if (length(idx2) > 0) filtered$r[idx2] <- NA
  }
  
  return(filtered)
}

# Apply filter
cor_res_filtered <- filter_cor_matrix(cor_res, dependent_pairs)

# Create filtered plot
p_filtered <- qcorrplot(cor_res_filtered, type = "lower", diag = FALSE) +
  geom_square() +
  scale_fill_gradientn(
    colours = heatmap_palette,
    limits = c(-1, 1),
    name = "Pearson r",
    guide = guide_colorbar(order = 1)
  ) +
  geom_couple(
    data = mantel_res,
    aes(colour = p_cat, size = r_cat, linetype = sign_cat),
    curvature = nice_curvature(),
    nudge_x = 0.4,
    nudge_y = 0.1
  ) +
  scale_colour_manual(
    values = mantel_p_colors,
    name = "Mantel p-value",
    guide = guide_legend(
      order = 2,
      override.aes = list(size = 2, linewidth = 1)
    )
  ) +
  scale_size_manual(
    values = c("< 0.1" = 0.5, "0.1 - 0.2" = 1.5, "≥ 0.2" = 2.5),
    name = "Mantel |r|",
    guide = guide_legend(order = 3)
  ) +
  scale_linetype_manual(
    values = c("Positive" = "solid", "Negative" = "dashed"),
    name = "Direction",
    guide = guide_legend(
      order = 4,
      override.aes = list(size = 1.5, linewidth = 1)
    )
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 9, face = "bold"),
    legend.position = "right",
    legend.key = element_blank(),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(
    title = "Plant Trait Associations (Biologically Independent Correlations Only)",
    subtitle = "Mathematically dependent correlations removed"
  )

# Save filtered version
ggsave("mantel_linkET_figure_filtered.pdf", plot = p_filtered, width = 14, height = 10, dpi = 300)
ggsave("mantel_linkET_figure_filtered.png", plot = p_filtered, width = 14, height = 10, dpi = 300)

cat("  Saved: mantel_linkET_figure_filtered (pdf/png) - math-dependent correlations removed\n\n")

# ============================================================================
# 9. SAVE NUMERICAL RESULTS
# ============================================================================

cat("Step 9: Saving numerical results...\n")

# Extract correlation matrix
cor_matrix <- cor(trait_df, method = "pearson")
write.csv(cor_matrix, "trait_correlation_matrix.csv")

# Save Mantel results
mantel_results_export <- mantel_res %>%
  select(env, r, p, p_cat, r_cat, sign_cat) %>%
  rename(
    Trait = env,
    Mantel_r = r,
    P_value = p,
    Significance = p_cat,
    Correlation_strength = r_cat,
    Direction = sign_cat
  )

write.csv(mantel_results_export, "mantel_test_results.csv", row.names = FALSE)

cat("  Saved: trait_correlation_matrix.csv\n")
cat("  Saved: mantel_test_results.csv\n\n")

# ============================================================================
# 10. SUMMARY
# ============================================================================

cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Summary:\n")
cat(sprintf("  %d plant traits analyzed\n", ncol(trait_df)))
cat(sprintf("  %d KO features included\n", ncol(species_df)))
cat(sprintf("  %d samples used\n\n", nrow(trait_df)))

cat("Significant trait-microbiome associations:\n")
sig_results <- mantel_res %>% filter(p < 0.05)
if (nrow(sig_results) > 0) {
  for (i in 1:nrow(sig_results)) {
    cat(sprintf("  %s: r = %.3f, p = %.4f\n",
                sig_results$env[i],
                sig_results$r[i],
                sig_results$p[i]))
  }
} else {
  cat("  No significant associations found (p < 0.05)\n")
}

cat("\nFiles created:\n")
cat("  - mantel_linkET_figure (pdf/png/svg) - ALL traits\n")
cat("  - mantel_linkET_figure_filtered (pdf/png) - biologically independent only\n")
cat("  - trait_correlation_matrix.csv\n")
cat("  - mantel_test_results.csv\n\n")

cat("Interpretation:\n")
cat("  The integrated figure shows:\n")
cat("  - LEFT: Correlation structure among plant traits\n")
cat("  - RIGHT: Mantel test associations between each trait and KO composition\n")
cat("  - Line color indicates p-value significance\n")
cat("  - Line thickness indicates correlation strength (|r|)\n")
cat("  - Line style indicates direction (solid = positive, dashed = negative)\n\n")

cat("Note: Mantel tests assess whether trait dissimilarity patterns\n")
cat("      are correlated with microbiome functional dissimilarity.\n\n")

cat("IMPORTANT: Some trait pairs are mathematically dependent:\n")
cat("  - N₂ Fixation per Plant = f(Nodule Biomass, Fixation per Nodule)\n")
cat("  - NMF = Nodule Biomass / Plant Biomass\n")
cat("  - RMF involves Plant Biomass in its calculation\n")
cat("These will show strong correlations due to mathematical relationships,\n")
cat("not necessarily biological relationships. Consider this when interpreting.\n\n")