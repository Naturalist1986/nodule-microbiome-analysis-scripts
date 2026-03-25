#!/usr/bin/env Rscript
# ============================================================================
# BIN ABUNDANCE-BASED ALPHA DIVERSITY BOXPLOTS BY TREATMENT
# ============================================================================
# One boxplot per diversity metric (Richness, Shannon, Simpson)
# Treatment colors and labels from Procrustes_Analysis.R
# ============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR <- "/mnt/c/Users/owner/My Drive (moshe.alon@mail.huji.ac.il)/Moshe/Efrat_Guy_Project"

BIN_ABUND_FILE <- file.path(BASE_DIR,
  "New_Binning_Pipeline_Coassembly/combined_abundance_long_renormalized.tsv")

OUTPUT_DIR <- file.path(BASE_DIR, "Boruta_New/bin_diversity_boruta")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTPUT_PREFIX <- file.path(OUTPUT_DIR, "bin_diversity_boxplots")

# ── Treatment colors and labels (from Procrustes_Analysis.R) ─────────────────
site_colors <- c(
  "carK" = "#66C2A5",   # CAR-C
  "carR" = "#FC8D62",   # CAR-G
  "ces"  = "#8DA0CB",   # CES
  "hok"  = "#A6D854",   # HUK
  "mtz"  = "#FFD92F",   # MTZ
  "RH"   = "#E5C494"    # RH
)

treatment_labels <- c(
  "carK" = "CAR-C",
  "carR" = "CAR-G",
  "ces"  = "CES",
  "hok"  = "HUK",
  "mtz"  = "MTZ",
  "RH"   = "RH"
)

# Ordered factor levels (controls first, then sites alphabetically)
treatment_order <- c("carK", "carR", "ces", "hok", "mtz", "RH")

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading bin abundance data...\n")
bin_abund <- read_tsv(BIN_ABUND_FILE, show_col_types = FALSE)
cat("  Rows:", nrow(bin_abund), "\n")
cat("  Columns:", paste(colnames(bin_abund), collapse = ", "), "\n\n")

# ============================================================================
# COMPUTE ALPHA DIVERSITY PER SAMPLE
# ============================================================================

cat("Computing alpha diversity metrics...\n")

bin_diversity <- bin_abund %>%
  group_by(Sample) %>%
  summarise(
    richness = sum(Relative_Abundance > 0),
    shannon  = {
      p <- Relative_Abundance[Relative_Abundance > 0] / 100
      -sum(p * log(p))
    },
    simpson  = {
      p <- Relative_Abundance[Relative_Abundance > 0] / 100
      1 - sum(p^2)
    },
    .groups = "drop"
  ) %>%
  rename(sample_id = Sample) %>%
  # Extract treatment from sample_id (e.g. "carK_1" -> "carK")
  mutate(
    treatment = sub("_\\d+$", "", sample_id),
    treatment = factor(treatment, levels = treatment_order)
  )

cat("  Samples computed:", nrow(bin_diversity), "\n")
cat("  Treatments found:", paste(sort(unique(bin_diversity$treatment)), collapse = ", "), "\n\n")
print(bin_diversity %>% arrange(treatment, sample_id))

# ============================================================================
# KRUSKAL-WALLIS + PAIRWISE WILCOXON (for significance annotation)
# ============================================================================

cat("\nKruskal-Wallis tests:\n")
kw_results <- list()
for (metric in c("richness", "shannon", "simpson")) {
  kw <- kruskal.test(reformulate("treatment", response = metric), data = bin_diversity)
  kw_results[[metric]] <- kw
  cat(sprintf("  %s: chi2 = %.3f, df = %d, p = %.4f\n",
              metric, kw$statistic, kw$parameter, kw$p.value))
}

# Helper: format p-value for annotation
fmt_p <- function(p) {
  if (p < 0.001) "p < 0.001"
  else if (p < 0.01) sprintf("p = %.3f", p)
  else sprintf("p = %.3f", p)
}

# ============================================================================
# PLOTTING FUNCTION
# ============================================================================

make_diversity_boxplot <- function(data, metric, y_label, title_label, kw_p) {

  # colours for treatments present in the data
  present_treatments <- levels(droplevels(data$treatment))
  col_vals  <- site_colors[present_treatments]
  lab_vals  <- treatment_labels[present_treatments]

  # y-position for KW annotation (5% above max)
  y_max  <- max(data[[metric]], na.rm = TRUE)
  y_min  <- min(data[[metric]], na.rm = TRUE)
  y_ann  <- y_max + (y_max - y_min) * 0.07

  p <- ggplot(data, aes(x = treatment, y = .data[[metric]],
                        fill = treatment, color = treatment)) +
    # Shaded box
    geom_boxplot(
      width         = 0.55,
      alpha         = 0.65,
      outlier.shape = NA,
      linewidth     = 0.6
    ) +
    # Individual points (jittered)
    geom_jitter(
      width  = 0.18,
      size   = 2.4,
      alpha  = 0.85,
      shape  = 21,
      fill   = "white",
      stroke = 1.1
    ) +
    # Kruskal-Wallis p annotation
    annotate(
      "text",
      x     = length(present_treatments) / 2 + 0.5,
      y     = y_ann,
      label = paste0("Kruskal-Wallis\n", fmt_p(kw_p)),
      size  = 3.4,
      color = "grey30",
      hjust = 0.5
    ) +
    scale_fill_manual(
      values = col_vals,
      labels = lab_vals,
      name   = "Treatment"
    ) +
    scale_color_manual(
      values = col_vals,
      labels = lab_vals,
      name   = "Treatment"
    ) +
    scale_x_discrete(labels = lab_vals) +
    coord_cartesian(ylim = c(y_min - (y_max - y_min) * 0.05,
                             y_ann + (y_max - y_min) * 0.08)) +
    labs(
      title = title_label,
      x     = NULL,
      y     = y_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x        = element_text(size = 11, face = "bold"),
      axis.text.y        = element_text(size = 10),
      axis.title.y       = element_text(size = 11),
      legend.position    = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank()
    )

  return(p)
}

# ============================================================================
# BUILD THREE INDIVIDUAL PLOTS
# ============================================================================

cat("\nBuilding plots...\n")

p_richness <- make_diversity_boxplot(
  bin_diversity, "richness",
  y_label     = "Bin Richness (count)",
  title_label = "Bin Richness",
  kw_p        = kw_results$richness$p.value
)

p_shannon <- make_diversity_boxplot(
  bin_diversity, "shannon",
  y_label     = "Shannon Index (H')",
  title_label = "Shannon Diversity",
  kw_p        = kw_results$shannon$p.value
)

p_simpson <- make_diversity_boxplot(
  bin_diversity, "simpson",
  y_label     = "Simpson Index (1 - D)",
  title_label = "Simpson Diversity",
  kw_p        = kw_results$simpson$p.value
)

# ============================================================================
# SAVE INDIVIDUAL PLOTS
# ============================================================================

for (info in list(
  list(plot = p_richness, name = "richness",  w = 6, h = 5),
  list(plot = p_shannon,  name = "shannon",   w = 6, h = 5),
  list(plot = p_simpson,  name = "simpson",   w = 6, h = 5)
)) {
  pdf(paste0(OUTPUT_PREFIX, "_", info$name, ".pdf"),
      width = info$w, height = info$h)
  print(info$plot)
  dev.off()

  png(paste0(OUTPUT_PREFIX, "_", info$name, ".png"),
      width = info$w, height = info$h,
      units = "in", res = 300)
  print(info$plot)
  dev.off()

  cat("  Saved:", info$name, "\n")
}

# ============================================================================
# COMBINED THREE-PANEL FIGURE
# ============================================================================

cat("\nBuilding combined figure...\n")

p_combined <- p_richness + p_shannon + p_simpson +
  plot_layout(ncol = 3) +
  plot_annotation(
    title    = "Bin Abundance-Based Alpha Diversity by Treatment",
    subtitle = paste0("n = ", nrow(bin_diversity), " samples | ",
                      "Kruskal-Wallis p shown on each panel"),
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

pdf(paste0(OUTPUT_PREFIX, "_combined.pdf"), width = 16, height = 5.5)
print(p_combined)
dev.off()

png(paste0(OUTPUT_PREFIX, "_combined.png"),
    width = 16, height = 5.5, units = "in", res = 300)
print(p_combined)
dev.off()

svg(paste0(OUTPUT_PREFIX, "_combined.svg"), width = 16, height = 5.5)
print(p_combined)
dev.off()

cat("  Saved: combined panel\n")

# ============================================================================
# DONE
# ============================================================================

cat("\n========================================================================\n")
cat("DONE\n")
cat("========================================================================\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Files produced:\n")
for (f in sort(list.files(OUTPUT_DIR, pattern = "bin_diversity_boxplots"))) {
  cat(" -", f, "\n")
}
cat("========================================================================\n")
