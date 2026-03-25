#!/usr/bin/env Rscript

# Genus-Level Taxonomy Barplot - Top 20 Taxa with Normalized Abundances
# Creates stacked barplot with good color contrast
# Author: Moshe
# Date: 2024

library(tidyverse)
library(RColorBrewer)
library(scales)
library(vegan)

# ===== Configuration =====
input_file <- "combined_by_sample_genus.tsv"
output_file <- "taxonomy_barplot_genus_top20.pdf"

# Number of top genera to show
n_top_genera <- 20

# Figure dimensions
fig_width <- 16
fig_height <- 9

# ===== Read and Process Data =====
cat("Reading genus data from", input_file, "...\n")

data <- read_tsv(input_file, show_col_types = FALSE)

# Convert to long format
data_long <- data %>%
  pivot_longer(
    cols = -taxonomy,
    names_to = "sample_full",
    values_to = "abundance"
  )

# Extract treatment and sample number
data_long <- data_long %>%
  mutate(
    treatment = str_extract(sample_full, "^[^_]+"),
    sample_num = str_extract(sample_full, "^[^_]+_[0-9]+"),
    sample = sample_num
  )

# Extract genus name
data_long <- data_long %>%
  mutate(
    genus = case_when(
      taxonomy == "unassigned" ~ "Unassigned",
      TRUE ~ str_extract(taxonomy, "g__[^;]+$") %>% 
        str_remove("g__")
    )
  )

# ===== NORMALIZE ABUNDANCES TO 100% PER SAMPLE =====
cat("Normalizing abundances to 100% per sample...\n")

# Check current sums
sample_sums <- data_long %>%
  group_by(sample) %>%
  summarise(total = sum(abundance, na.rm = TRUE))

cat("Abundance sums before normalization (should be ~100):\n")
print(summary(sample_sums$total))

# Normalize to 100%
data_long <- data_long %>%
  group_by(sample) %>%
  mutate(abundance = (abundance / sum(abundance, na.rm = TRUE)) * 100) %>%
  ungroup()

# Verify normalization
sample_sums_after <- data_long %>%
  group_by(sample) %>%
  summarise(total = sum(abundance, na.rm = TRUE))

cat("\nAbundance sums after normalization:\n")
print(summary(sample_sums_after$total))

# ===== Calculate mean abundance and select top genera =====
genus_means <- data_long %>%
  group_by(genus) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE)) %>%
  arrange(desc(mean_abundance))

cat("\nTotal genera detected:", nrow(genus_means), "\n")

# Select top N genera
top_genera <- genus_means %>% 
  head(n_top_genera) %>% 
  pull(genus)

cat("Showing top", n_top_genera, "genera, grouping rest as 'Other'\n")

# Group low-abundance genera into "Other"
plot_data <- data_long %>%
  mutate(
    genus_plot = if_else(genus %in% top_genera, genus, "Other (< top 20)")
  ) %>%
  group_by(treatment, sample, genus_plot) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Extract replicate number for x-axis labels
plot_data <- plot_data %>%
  mutate(replicate = str_extract(sample, "[0-9]+$"))

# Set factor levels for treatment order
treatment_order <- c("carK", "carR", "ces", "hok", "mtz", "RH")

# Treatment labels — abbreviated, matching Figure 1 and Figure 2 barplot
treatment_labels <- c(
  "carK" = "CAR-C",
  "carR" = "CAR-G",
  "ces"  = "CES",
  "hok"  = "HUK",
  "mtz"  = "MTZ",
  "RH"   = "RH"
)

# Unified treatment colors — matching make_figure.R (Procrustes/Figure 1)
treatment_colors <- c(
  "CAR-C" = "#66C2A5",
  "CAR-G" = "#FC8D62",
  "CES"   = "#8DA0CB",
  "HUK"   = "#A6D854",
  "MTZ"   = "#FFD92F",
  "RH"    = "#E5C494"
)

# Unified treatment shapes — all circles
treatment_shapes <- c(
  "CAR-C" = 16,
  "CAR-G" = 16,
  "CES"   = 16,
  "HUK"   = 16,
  "MTZ"   = 16,
  "RH"    = 16
)

plot_data$treatment <- factor(plot_data$treatment, 
                              levels = treatment_order,
                              labels = treatment_labels[treatment_order])

# Order genera by abundance (with "Other" last)
genus_order <- c(
  genus_means %>% 
    filter(genus %in% top_genera) %>% 
    pull(genus),
  "Other (< top 20)"
)
plot_data$genus_plot <- factor(plot_data$genus_plot, levels = genus_order)

# ===== Generate High-Contrast Color Palette =====
n_colors <- length(unique(plot_data$genus_plot)) - 1  # Exclude "Other"

cat("Generating", n_colors + 1, "high-contrast colors\n")

# Use a combination of distinct palettes
colors <- c(
  # Bradyrhizobium - coral/salmon
  "#FB8072",
  # Other major genera - varied distinct colors
  "#80B1D3",  # light blue
  "#FDB462",  # light orange
  "#B3DE69",  # olive/tan
  "#BEBADA",  # lavender
  "#8DD3C7",  # cyan
  "#FCCDE5",  # pink
  "#D9D9D9",  # gray
  "#BC80BD",  # purple
  "#CCEBC5",  # light green
  "#FFED6F",  # yellow
  "#FFFFCC",  # light yellow/cream
  # Additional colors
  "#E41A1C",  # red
  "#377EB8",  # blue
  "#4DAF4A",  # green
  "#984EA3",  # purple
  "#FF7F00",  # orange
  "#A65628",  # brown
  "#F781BF",  # pink
  "#999999"   # gray
)

# Ensure we have enough colors
if (length(colors) < n_colors) {
  additional <- hcl.colors(n_colors - length(colors), palette = "Set 2")
  colors <- c(colors, additional)
}

# Take only what we need and add gray for "Other"
colors <- c(colors[1:n_colors], "gray75")
names(colors) <- levels(plot_data$genus_plot)

# ===== Create Plot =====
cat("Creating plot...\n")

p <- ggplot(plot_data, aes(x = replicate, y = abundance, fill = genus_plot)) +
  geom_bar(stat = "identity", position = "stack", width = 0.9) +
  facet_grid(~ treatment, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colors, name = "Genus") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  labs(
    title = "Taxonomic Composition - All Samples by Treatment (Top 20 Genera)",
    x = "Replicate",
    y = "Relative Abundance (%)"
  ) +
  theme_classic() +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    panel.spacing = unit(0.5, "lines")
  ) +
  guides(fill = guide_legend(ncol = 1))

# ===== Save Plot =====
cat("Saving plot to", output_file, "...\n")

ggsave(
  filename = output_file,
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300
)

cat("Done! Plot saved to", output_file, "\n")

# ===== Print summary statistics =====
cat("\n=== Summary Statistics ===\n")
cat("Total samples:", length(unique(plot_data$sample)), "\n")
cat("Genera shown individually:", n_colors, "\n")
cat("Treatments:", paste(unique(plot_data$treatment), collapse = ", "), "\n")

cat("\n--- Top 20 Most Abundant Genera ---\n")
print(genus_means %>% head(20), n = 20)

# Calculate "Other" percentage
other_percent <- plot_data %>%
  filter(genus_plot == "Other (< top 20)") %>%
  summarise(mean_other = mean(abundance)) %>%
  pull(mean_other)

cat("\nMean abundance in 'Other' category:", round(other_percent, 2), "%\n")

# ===== Save PNG version =====
ggsave(
  filename = str_replace(output_file, ".pdf", ".png"),
  plot = p,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300
)

cat("Also saved PNG version\n")

# ===== ORDINATION ANALYSIS =====
cat("\n=== Performing Ordination Analysis ===\n")

# Prepare data for ordination: sample x genus matrix
cat("Preparing sample x genus matrix...\n")

ordination_data <- data_long %>%
  group_by(sample, genus) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = genus,
    values_from = abundance,
    values_fill = 0
  )

# Extract sample names and create metadata
sample_names <- ordination_data$sample
metadata <- data.frame(
  sample = sample_names,
  treatment = str_extract(sample_names, "^[^_]+")
) %>%
  mutate(treatment = factor(treatment, 
                            levels = treatment_order,
                            labels = treatment_labels[treatment_order]))

# Create abundance matrix (samples as rows, genera as columns)
abundance_matrix <- ordination_data %>%
  select(-sample) %>%
  as.matrix()

rownames(abundance_matrix) <- sample_names

cat("Matrix dimensions:", nrow(abundance_matrix), "samples x", 
    ncol(abundance_matrix), "genera\n")

# Calculate Bray-Curtis distance
cat("Calculating Bray-Curtis dissimilarity matrix...\n")
bray_dist <- vegdist(abundance_matrix, method = "bray")

# --- PCoA (Principal Coordinates Analysis) ---
cat("Performing PCoA...\n")
pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Extract variance explained
pcoa_var <- round(pcoa_result$eig / sum(pcoa_result$eig) * 100, 2)

# Create PCoA data frame
pcoa_df <- data.frame(
  sample = sample_names,
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2]
) %>%
  left_join(metadata, by = "sample")

# --- PERMANOVA (run before plotting so results can be annotated) ---
cat("Performing PERMANOVA...\n")
permanova_result <- adonis2(bray_dist ~ treatment, data = metadata, permutations = 999)
perm_r2 <- round(permanova_result$R2[1], 3)
perm_p  <- permanova_result$`Pr(>F)`[1]
perm_p_label <- ifelse(perm_p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", perm_p)))
permanova_label <- paste0("PERMANOVA\nR\u00b2 = ", perm_r2, ", ", perm_p_label)
cat("PERMANOVA R2:", perm_r2, " p:", perm_p, "\n")

# Create PCoA plot
cat("Creating PCoA plot...\n")

# Calculate convex hulls for each treatment
pcoa_hulls <- pcoa_df %>%
  group_by(treatment) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = treatment,
                                 shape = treatment)) +
  geom_polygon(data = pcoa_hulls, aes(fill = treatment), alpha = 0.15,
               show.legend = FALSE) +
  geom_point(size = 4, alpha = 0.8) +
  annotate("text", x = Inf, y = Inf, label = permanova_label,
           hjust = 1.05, vjust = 1.4, size = 3.5, fontface = "italic",
           colour = "grey30") +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_fill_manual(values  = treatment_colors, guide = "none") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  labs(
    title = "PCoA - Bacterial Community Composition (Bray-Curtis)",
    x = paste0("PC1 (", pcoa_var[1], "%)"),
    y = paste0("PC2 (", pcoa_var[2], "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Save PCoA plot
pcoa_file <- str_replace(output_file, ".pdf", "_PCoA.pdf")
ggsave(
  filename = pcoa_file,
  plot = pcoa_plot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)
cat("PCoA plot saved to:", pcoa_file, "\n")

# Save PNG version
ggsave(
  filename = str_replace(pcoa_file, ".pdf", ".png"),
  plot = pcoa_plot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

# Save SVG version
ggsave(
  filename = str_replace(pcoa_file, ".pdf", ".svg"),
  plot = pcoa_plot,
  width = 10,
  height = 7,
  units = "in"
)
cat("PCoA plot also saved as PNG and SVG\n")

# --- NMDS (Non-metric Multidimensional Scaling) ---
cat("Performing NMDS...\n")
set.seed(123)  # For reproducibility
nmds_result <- metaMDS(abundance_matrix, distance = "bray", k = 2, 
                       trymax = 100, trace = FALSE)

cat("NMDS stress:", round(nmds_result$stress, 3), "\n")
if (nmds_result$stress > 0.2) {
  cat("Warning: NMDS stress > 0.2, interpretation should be cautious\n")
}

# Create NMDS data frame
nmds_df <- data.frame(
  sample = sample_names,
  NMDS1 = nmds_result$points[, 1],
  NMDS2 = nmds_result$points[, 2]
) %>%
  left_join(metadata, by = "sample")

# Create NMDS plot
cat("Creating NMDS plot...\n")

# Calculate convex hulls for each treatment
nmds_hulls <- nmds_df %>%
  group_by(treatment) %>%
  slice(chull(NMDS1, NMDS2)) %>%
  ungroup()

nmds_plot <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = treatment,
                                 shape = treatment)) +
  geom_polygon(data = nmds_hulls, aes(fill = treatment), alpha = 0.15,
               show.legend = FALSE) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = treatment_colors, name = "Treatment") +
  scale_fill_manual(values  = treatment_colors, guide = "none") +
  scale_shape_manual(values = treatment_shapes, name = "Treatment") +
  labs(
    title = paste0("NMDS - Bacterial Community Composition (Stress = ", 
                   round(nmds_result$stress, 3), ")"),
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Save NMDS plot
nmds_file <- str_replace(output_file, ".pdf", "_NMDS.pdf")
ggsave(
  filename = nmds_file,
  plot = nmds_plot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)
cat("NMDS plot saved to:", nmds_file, "\n")

# Save PNG version
ggsave(
  filename = str_replace(nmds_file, ".pdf", ".png"),
  plot = nmds_plot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300
)

# Save SVG version
ggsave(
  filename = str_replace(nmds_file, ".pdf", ".svg"),
  plot = nmds_plot,
  width = 10,
  height = 7,
  units = "in"
)
cat("NMDS plot also saved as PNG and SVG\n")

cat("\nPERMANOVA Results:\n")
print(permanova_result)

# Save PERMANOVA results
permanova_file <- str_replace(output_file, ".pdf", "_PERMANOVA.txt")
sink(permanova_file)
cat("=== PERMANOVA Results ===\n")
cat("Testing for differences in bacterial community composition by treatment\n")
cat("Distance metric: Bray-Curtis\n")
cat("Permutations: 999\n\n")
print(permanova_result)
sink()
cat("\nPERMANOVA results saved to:", permanova_file, "\n")

# --- Beta-dispersion test ---
cat("\nTesting homogeneity of dispersion (betadisper)...\n")
betadisp_result <- betadisper(bray_dist, metadata$treatment)
betadisp_test <- permutest(betadisp_result, permutations = 999)
cat("\nBeta-dispersion test (PERMDISP):\n")
print(betadisp_test)

# Save beta-dispersion results
betadisp_file <- str_replace(output_file, ".pdf", "_betadispersion.txt")
sink(betadisp_file)
cat("=== Beta-dispersion Test (PERMDISP) ===\n")
cat("Testing homogeneity of multivariate dispersions\n")
cat("This tests whether variation within groups is similar\n")
cat("Significant result suggests different dispersion among treatments\n\n")
print(betadisp_test)
cat("\n\nAverage distances to centroids by treatment:\n")
print(betadisp_result$group.distances)
sink()
cat("Beta-dispersion results saved to:", betadisp_file, "\n")

# Save ordination coordinates
pcoa_coords_file <- str_replace(output_file, ".pdf", "_PCoA_coordinates.txt")
write_tsv(pcoa_df, pcoa_coords_file)
cat("PCoA coordinates saved to:", pcoa_coords_file, "\n")

nmds_coords_file <- str_replace(output_file, ".pdf", "_NMDS_coordinates.txt")
write_tsv(nmds_df, nmds_coords_file)
cat("NMDS coordinates saved to:", nmds_coords_file, "\n")

cat("\n=== Ordination Analysis Complete ===\n\n")

# ===== Copy PCoA SVG to Final_Figures/Figure_2 =====
FINAL_FIG2 <- "G:/My Drive/Moshe/Efrat_Guy_Project/Final_Figures/Figure_2"
dir.create(FINAL_FIG2, showWarnings = FALSE, recursive = TRUE)
pcoa_svg <- str_replace(pcoa_file, ".pdf", ".svg")
if (file.exists(pcoa_svg)) {
  file.copy(pcoa_svg,
            file.path(FINAL_FIG2, "Figure_2D_Genus_PCoA.svg"), overwrite = TRUE)
  cat("Copied PCoA SVG to Final_Figures/Figure_2/Figure_2D_Genus_PCoA.svg\n")
}
if (file.exists(pcoa_file)) {
  file.copy(pcoa_file,
            file.path(FINAL_FIG2, "Figure_2D_Genus_PCoA.pdf"), overwrite = TRUE)
  cat("Copied PCoA PDF to Final_Figures/Figure_2/Figure_2D_Genus_PCoA.pdf\n")
}

# ===== Create genus list file =====
genus_list_file <- str_replace(output_file, ".pdf", "_genus_list.txt")
write_tsv(
  genus_means %>% 
    mutate(
      rank = row_number(),
      shown_in_plot = if_else(rank <= n_top_genera, "Yes", "No (in Other)")
    ),
  genus_list_file
)
cat("Complete genus list saved to:", genus_list_file, "\n")

# ===== Save treatment summary =====
treatment_summary <- plot_data %>%
  filter(genus_plot != "Other (< top 20)") %>%
  group_by(treatment, genus_plot) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop") %>%
  arrange(treatment, desc(mean_abundance))

treatment_summary_file <- str_replace(output_file, ".pdf", "_treatment_summary.txt")
write_tsv(treatment_summary, treatment_summary_file)
cat("Treatment summary saved to:", treatment_summary_file, "\n")