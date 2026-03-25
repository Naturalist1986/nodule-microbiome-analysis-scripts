#!/usr/bin/env Rscript

setwd("G:/My Drive/Moshe/Efrat_Guy_Project/Procrustes")

# =============================================================================
# make_figure.R — Master Figure Assembly Script
# Consolidates: Archive/plant_trait_boxplots.R, Mantel.R,
#               Procrustes_Analysis.R, combine_four_panel_egg.R
#
# 5-panel layout:
#   Row 1: A (inoculation bar chart) | B (plant trait boxplots)
#   Rows 2-3: C (Mantel figure) | D (Plant Traits PCA) | E (KO PCoA)
#
# Expected output files:
#   plant_traits_final_formatted.png
#   mantel_linkET_figure.png
#   original_ordinations_combined.png
#   combined_four_panel_egg.pdf
#   combined_four_panel_egg.png
# =============================================================================

# ── LIBRARIES ──────────────────────────────────────────────────────────────────
library(readxl)
library(ggplot2)
library(tidyverse)   # includes dplyr, tidyr, stringr, purrr
library(gridExtra)
library(grid)
library(png)
library(patchwork)
library(car)
library(multcomp)
library(emmeans)
library(multcompView)
library(vegan)
library(RColorBrewer)

# Install linkET if not already installed
if (!require(linkET, quietly = TRUE)) {
  if (!require(devtools, quietly = TRUE)) install.packages("devtools")
  devtools::install_github("Hy4m/linkET")
  library(linkET)
}

# ── SHARED CONSTANTS ───────────────────────────────────────────────────────────

# Site colors by display name (includes DGN / YATIR for Panel A)
site_colors <- c(
  "CAR-C" = "#66C2A5",
  "CAR-G" = "#FC8D62",
  "CES"   = "#8DA0CB",
  "DGN"   = "#E78AC3",
  "HUK"   = "#A6D854",
  "MTZ"   = "#FFD92F",
  "RH"    = "#E5C494",
  "YATIR" = "#B3B3B3"
)

# Site colors by raw treatment code (for ordination panels D & E)
site_colors_raw <- c(
  "carK" = "#66C2A5",   # CAR-C
  "carR" = "#FC8D62",   # CAR-G
  "ces"  = "#8DA0CB",   # CES
  "hok"  = "#A6D854",   # HUK
  "mtz"  = "#FFD92F",   # MTZ
  "RH"   = "#E5C494"    # RH
)

# Treatment shapes — all circles
site_shapes_raw <- c(
  "carK" = 16,
  "carR" = 16,
  "ces"  = 16,
  "hok"  = 16,
  "mtz"  = 16,
  "RH"   = 16
)

# Treatment display labels and ordering
treatment_labels <- c(
  "carK" = "CAR-C",
  "carR" = "CAR-G",
  "ces"  = "CES",
  "hok"  = "HUK",
  "mtz"  = "MTZ",
  "RH"   = "RH"
)
treatment_order <- c("carK", "carR", "ces", "hok", "mtz", "RH")

# Panel tag theme — places bold letter in top-left margin, outside the plot area
tag_theme <- theme(
  plot.tag          = element_text(face = "bold", size = 16),
  plot.tag.position = "topleft"
)

# ══ SECTION 1: LOAD SHARED DATA ═══════════════════════════════════════════════
cat("================================================================================\n")
cat("SECTION 1: LOADING SHARED DATA\n")
cat("================================================================================\n\n")

plant_data_shared <- read_excel("Guy_measurements.xlsx")
ko_data <- read.delim("normalized_kegg_results.tsv",
                      row.names = 1, check.names = FALSE)

cat(sprintf("  Plant data (Guy_measurements.xlsx): %d samples x %d traits\n",
            nrow(plant_data_shared), ncol(plant_data_shared) - 2))
cat(sprintf("  KO data (normalized_kegg_results.tsv): %d KOs x %d metagenomes\n\n",
            nrow(ko_data), ncol(ko_data)))

# ══ SECTION 2: PANEL B — Plant Trait Boxplots ═════════════════════════════════
cat("================================================================================\n")
cat("SECTION 2: PANEL B — Plant Trait Boxplots\n")
cat("================================================================================\n\n")

# Prepare data (add display Treatment and NMF %)
data_traits <- plant_data_shared %>%
  mutate(
    Treatment_Full = factor(Treatment,
                            levels = names(treatment_labels),
                            labels = treatment_labels),
    NMF_percent = NMF * 100
  )

# Sub-panel configuration for the four traits
traits_config <- list(
  list(col = "percent_N",           title = "Leaf N",
       ylab = "Leaf N (%)",                                             tag = "(a)"),
  list(col = "Plant_Biomass",       title = "Plant Biomass",
       ylab = "Plant Biomass (g)",                                      tag = "(b)"),
  list(col = "Fixation_per_Nodule", title = "Fixation rate",
       ylab = expression(paste(mu, "mol C"[2], "H"[4], " h"^-1, " g"^-1)), tag = "(c)"),
  list(col = "NMF_percent",         title = "NMF",
       ylab = "Nodule mass fraction (%)",                               tag = "(d)")
)

process_trait <- function(config, raw_data) {
  trait_col    <- config$col
  formula_obj  <- as.formula(paste(trait_col, "~ Treatment_Full"))
  model        <- aov(formula_obj, data = raw_data)
  emm          <- emmeans(model, ~ Treatment_Full)
  cld_obj      <- cld(emm, Letters = letters, adjust = "tukey", decreasing = TRUE)

  cld_df <- as.data.frame(cld_obj) %>%
    mutate(
      letters = trimws(.group),
      y_pos   = max(raw_data[[trait_col]], na.rm = TRUE) * 1.08
    )

  p <- ggplot(raw_data,
              aes(x = Treatment_Full, y = .data[[trait_col]],
                  fill = Treatment_Full)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7) +
    geom_point(position = position_jitter(width = 0.15),
               size = 1.2, alpha = 0.4) +
    # Tukey letters above boxes
    geom_text(data = cld_df,
              aes(x = Treatment_Full, y = y_pos, label = letters),
              size = 6, fontface = "bold", inherit.aes = FALSE, vjust = 0) +
    # Sub-panel letter (a/b/c/d) — stays as annotate within the subplot
    annotate("text", x = -Inf, y = Inf, label = config$tag,
             vjust = 1.5, hjust = -0.5, fontface = "bold", size = 8) +
    scale_fill_manual(values = site_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(title = config$title, y = config$ylab) +
    theme_bw() +
    theme(
      axis.title.x       = element_blank(),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(hjust = 0.5, face = "bold", size = 17),
      axis.title.y       = element_text(hjust = 0.5, face = "bold", size = 17),
      legend.position    = "none",
      legend.text        = element_text(hjust = 0.5, face = "bold", size = 17),
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 16),
      panel.grid.major.x = element_blank()
    )

  return(list(plot = p))
}

results_traits   <- lapply(traits_config, process_trait, raw_data = data_traits)
all_trait_plots  <- lapply(results_traits, `[[`, "plot")

# Legend extraction
legend_plot <- ggplot(data_traits,
                      aes(x = Treatment_Full, y = percent_N, fill = Treatment_Full)) +
  geom_boxplot() +
  scale_fill_manual(values = site_colors) +
  labs(fill = "Treatment Site") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title    = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.text     = element_text(size = 13)
  ) +
  guides(fill = guide_legend(ncol = 1))

get_leg <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

traits_grid <- grid.arrange(
  do.call(arrangeGrob, c(all_trait_plots, ncol = 4)),
  get_leg(legend_plot),
  ncol   = 2,
  widths = c(10, 2.5)
)

# Output directory for Figure 1 panels (vector files)
FINAL_FIG1 <- "G:/My Drive/Moshe/Efrat_Guy_Project/Final_Figures/Figure_1"
dir.create(FINAL_FIG1, showWarnings = FALSE, recursive = TRUE)

# Save standalone PNG
ggsave("plant_traits_final_formatted.png", traits_grid,
       width = 18, height = 5.5, dpi = 300)
cat("  Saved: plant_traits_final_formatted.png\n\n")

# Save Panel B as true vector to Final_Figures
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
library(svglite)
svglite(file.path(FINAL_FIG1, "Figure_1B_PlantTraits_Panel.svg"), width = 18, height = 5.5)
grid.draw(traits_grid)
dev.off()
pdf(file.path(FINAL_FIG1, "Figure_1B_PlantTraits_Panel.pdf"), width = 18, height = 5.5)
grid.draw(traits_grid)
dev.off()
cat("  Saved: Figure_1B_PlantTraits_Panel.svg/pdf\n\n")

# Use traits_grid directly (no PNG round-trip) for combined figure
panel_b_grob <- gridExtra::arrangeGrob(
  traits_grid,
  top = grid::textGrob("B", x = unit(0.01, "npc"), hjust = 0,
                       gp = grid::gpar(fontface = "bold", fontsize = 16))
)

# ══ SECTION 3: PANEL C — Mantel Test Figure ═══════════════════════════════════
cat("================================================================================\n")
cat("SECTION 3: PANEL C — Mantel Test Figure\n")
cat("================================================================================\n\n")

# Mantel uses the *older* plant data file (Guy_measurements_old.xlsx)
plant_data_old <- read_excel("Guy_measurements_old.xlsx")
cat(sprintf("  Loaded Guy_measurements_old.xlsx: %d rows\n", nrow(plant_data_old)))

# --- Merge duplicate metagenome runs (using shared ko_data) ---
sample_info_m <- data.frame(column = colnames(ko_data), stringsAsFactors = FALSE) %>%
  mutate(
    treatment = str_extract(column, "^[^_]+"),
    replicate = as.numeric(str_extract(column, "(?<=_)\\d+")),
    sample_id = paste(treatment, replicate, sep = "_")
  )

ko_tr_m <- as.data.frame(t(ko_data))
ko_tr_m$sample_id <- sample_info_m$sample_id

ko_merged_m <- ko_tr_m %>%
  group_by(sample_id) %>%
  summarise(across(everything(), mean, .names = "{.col}"))

ids_m       <- ko_merged_m$sample_id
ko_merged_m <- ko_merged_m %>% dplyr::select(-sample_id) %>% as.data.frame()
ko_merged_m <- as.data.frame(t(ko_merged_m))
colnames(ko_merged_m) <- ids_m

# --- Match samples ---
plant_data_old <- plant_data_old %>%
  mutate(sample_id = paste(Treatment, Replicate, sep = "_"))

common_m      <- intersect(colnames(ko_merged_m), plant_data_old$sample_id)
cat(sprintf("  Common samples: %d\n", length(common_m)))

plant_matched_m <- plant_data_old %>%
  filter(sample_id %in% common_m) %>%
  arrange(sample_id)
ko_matched_m <- ko_merged_m[, plant_matched_m$sample_id]

# --- Select traits and remove missing values ---
trait_cols_m <- c("percent_N", "Height", "Nodule_Biomass", "Plant_Biomass",
                  "Fixation_per_Nodule", "Fixation_per_Plant", "RMF", "NMF")

plant_matrix_m <- plant_matched_m %>%
  dplyr::select(all_of(trait_cols_m)) %>%
  as.matrix()
rownames(plant_matrix_m) <- plant_matched_m$sample_id

complete_m       <- complete.cases(plant_matrix_m)
plant_matrix_m   <- plant_matrix_m[complete_m, ]
ko_matched_m     <- ko_matched_m[, complete_m]

# --- Filter KOs ---
n_s_m         <- ncol(ko_matched_m)
min_s_m       <- ceiling(n_s_m * 0.10)
ko_prev_m     <- rowSums(ko_matched_m > 0)
ko_mean_m     <- rowMeans(ko_matched_m)
ko_keep_m     <- (ko_prev_m >= min_s_m) & (ko_mean_m >= 0.0001)
ko_filtered_m <- ko_matched_m[ko_keep_m, ]
species_df_m  <- as.data.frame(t(ko_filtered_m))   # samples × KOs

cat(sprintf("  KOs retained: %d (%.1f%%)\n\n",
            sum(ko_keep_m), sum(ko_keep_m) / nrow(ko_matched_m) * 100))

# --- Trait labels and data frame for linkET ---
trait_label_map_m <- c(
  "percent_N"           = "Leaf N (%)",
  "Height"              = "Plant Height (cm)",
  "Nodule_Biomass"      = "Nodule Biomass (g)",
  "Plant_Biomass"       = "Plant Biomass (g)",
  "Fixation_per_Nodule" = "N\u2082 Fixation per Nodule",
  "Fixation_per_Plant"  = "N\u2082 Fixation per Plant",
  "RMF"                 = "Root Mass Fraction",
  "NMF"                 = "Nodule Mass Fraction"
)

trait_df_m             <- as.data.frame(plant_matrix_m)
colnames(trait_df_m)   <- trait_label_map_m[colnames(trait_df_m)]

# --- Correlations and Mantel tests ---
cat("  Calculating correlations and Mantel tests...\n")
cor_res_m  <- correlate(trait_df_m, method = "pearson")

mantel_res_m <- mantel_test(
  species_df_m, trait_df_m,
  spec_select = list(`Functional composition` = 1:ncol(species_df_m)),
  mantel_fun  = "mantel"
) %>%
  mutate(
    p_cat    = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                   labels = c("< 0.01", "0.01 - 0.05", "\u2265 0.05"),
                   right  = FALSE),
    r_cat    = cut(abs(r), breaks = c(-Inf, 0.1, 0.2, Inf),
                   labels = c("< 0.1", "0.1 - 0.2", "\u2265 0.2"),
                   right  = FALSE),
    sign_cat = ifelse(r > 0, "Positive", "Negative")
  )

cat("  Mantel tests complete\n\n")

# --- Build Mantel figure ---
heatmap_pal_m  <- rev(brewer.pal(n = 11, name = "RdBu"))
mantel_p_col_m <- c("< 0.01"       = "#D95F02",
                    "0.01 - 0.05"  = "#7570B3",
                    "\u2265 0.05"  = "grey70")

p_mantel <- qcorrplot(cor_res_m, type = "lower", diag = FALSE) +
  geom_tile(aes(fill = r), colour = "white", linewidth = 0.5) +  # geom_tile = vector; geom_square uses rasterGrob
  scale_fill_gradientn(colours  = heatmap_pal_m,
                       limits   = c(-1, 1),
                       name     = "Pearson r",
                       guide    = guide_colorbar(order = 1)) +
  geom_couple(data      = mantel_res_m,
              aes(colour = p_cat, size = r_cat, linetype = sign_cat),
              curvature  = nice_curvature(),
              nudge_x    = 0.1,
              label.size = 5) +
  scale_colour_manual(values = mantel_p_col_m,
                      name   = "Mantel p-value",
                      guide  = guide_legend(order        = 2,
                                            override.aes = list(size = 2,
                                                                linewidth = 1))) +
  scale_size_manual(values = c("< 0.1" = 0.5, "0.1 - 0.2" = 1.5, "\u2265 0.2" = 2.5),
                    name   = "Mantel |r|",
                    guide  = guide_legend(order = 3)) +
  scale_linetype_manual(values = c("Positive" = "solid", "Negative" = "dashed"),
                        guide  = "none") +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 14),
    axis.text.y        = element_text(size = 14),
    axis.text.y.right  = element_text(size = 24),
    legend.position    = "right",
    legend.key         = element_blank(),
    legend.title       = element_text(size = 14, face = "bold"),
    legend.text        = element_text(size = 14),
    plot.title         = element_text(size = 15, face = "bold"),
    plot.subtitle      = element_text(size = 14),
    plot.margin        = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.background    = element_rect(fill = "white", colour = NA)
  ) +
  labs(
    title    = "Plant Trait Correlations and Associations with Nodule Microbiome Function",
    subtitle = "Left: Pearson correlations between traits | Right: Mantel tests with KO composition"
  )

# Save Mantel figure and numerical results
ggsave("mantel_linkET_figure.pdf", plot = p_mantel, width = 14, height = 10, dpi = 300)
ggsave("mantel_linkET_figure.png", plot = p_mantel, width = 14, height = 10, dpi = 300)
svglite("mantel_linkET_figure.svg", width = 14, height = 10)
print(p_mantel); dev.off()
svglite(file.path(FINAL_FIG1, "Figure_1C_Mantel_Panel.svg"), width = 14, height = 10)
print(p_mantel); dev.off()
ggsave(file.path(FINAL_FIG1, "Figure_1C_Mantel_Panel.pdf"), plot = p_mantel, width = 14, height = 10)
cat("  Saved: Figure_1C_Mantel_Panel.svg/pdf\n")

cor_matrix_m <- cor(trait_df_m, method = "pearson")
write.csv(cor_matrix_m, "trait_correlation_matrix.csv")

mantel_export_m <- mantel_res_m %>%
  dplyr::select(env, r, p, p_cat, r_cat, sign_cat) %>%
  rename(Trait               = env,
         Mantel_r            = r,
         P_value             = p,
         Significance        = p_cat,
         Correlation_strength = r_cat,
         Direction           = sign_cat)
write.csv(mantel_export_m, "mantel_test_results.csv", row.names = FALSE)

cat("  Saved: mantel_linkET_figure.png\n")
cat("  Saved: trait_correlation_matrix.csv\n")
cat("  Saved: mantel_test_results.csv\n\n")

# Use p_mantel directly as panel_c (no PNG round-trip)
panel_c <- p_mantel + labs(tag = "C") + tag_theme

# ══ SECTION 4: PANELS D & E — Procrustes + Ordinations ════════════════════════
cat("================================================================================\n")
cat("SECTION 4: PANELS D & E — Procrustes Analysis + Ordinations\n")
cat("================================================================================\n\n")

# --- Merge duplicate metagenome runs (using shared ko_data) ---
sample_info_p <- data.frame(column = colnames(ko_data), stringsAsFactors = FALSE) %>%
  mutate(
    treatment = str_extract(column, "^[^_]+"),
    replicate = as.numeric(str_extract(column, "(?<=_)\\d+")),
    sample_id = paste(treatment, replicate, sep = "_")
  )

ko_tr_p <- as.data.frame(t(ko_data))
ko_tr_p$sample_id <- sample_info_p$sample_id

ko_merged_p <- ko_tr_p %>%
  group_by(sample_id) %>%
  summarise(across(everything(), mean, .names = "{.col}"))

ids_p       <- ko_merged_p$sample_id
ko_merged_p <- ko_merged_p %>% dplyr::select(-sample_id) %>% as.data.frame()
ko_merged_p <- as.data.frame(t(ko_merged_p))
colnames(ko_merged_p) <- ids_p

# --- Match samples (uses Guy_measurements.xlsx via plant_data_shared) ---
plant_data_proc <- plant_data_shared %>%
  mutate(sample_id = paste(Treatment, Replicate, sep = "_"))

common_p      <- intersect(colnames(ko_merged_p), plant_data_proc$sample_id)
cat(sprintf("  Common samples: %d\n", length(common_p)))

plant_matched_p <- plant_data_proc %>%
  filter(sample_id %in% common_p) %>%
  arrange(sample_id)
ko_matched_p <- ko_merged_p[, plant_matched_p$sample_id]
stopifnot(all(colnames(ko_matched_p) == plant_matched_p$sample_id))
cat("  Sample order verified\n")

# --- Select 4 traits, remove missing values ---
trait_cols_p <- c("percent_N", "Plant_Biomass", "Fixation_per_Nodule", "NMF")

plant_matrix_p <- plant_matched_p %>%
  dplyr::select(all_of(trait_cols_p)) %>%
  as.matrix()
rownames(plant_matrix_p) <- plant_matched_p$sample_id

complete_p      <- complete.cases(plant_matrix_p)
plant_matrix_p  <- plant_matrix_p[complete_p, ]
ko_matched_p    <- ko_matched_p[, complete_p]
plant_matched_p <- plant_matched_p[complete_p, ]

# --- Filter KOs ---
n_s_p         <- ncol(ko_matched_p)
min_s_p       <- ceiling(n_s_p * 0.10)
ko_prev_p     <- rowSums(ko_matched_p > 0)
ko_mean_p     <- rowMeans(ko_matched_p)
ko_keep_p     <- (ko_prev_p >= min_s_p) & (ko_mean_p >= 0.0001)
ko_filtered_p <- ko_matched_p[ko_keep_p, ]
ko_tr_filt_p  <- t(ko_filtered_p)   # samples × KOs

cat(sprintf("  KOs retained: %d (%.1f%%)\n",
            sum(ko_keep_p), sum(ko_keep_p) / nrow(ko_matched_p) * 100))
cat(sprintf("  Final samples: %d\n\n", nrow(plant_matrix_p)))

# --- Plant traits PCA (standardised) ---
plant_std_p  <- scale(plant_matrix_p)
pca_traits_p <- rda(plant_std_p)
plant_var_p  <- summary(eigenvals(pca_traits_p))

# --- Procrustes analysis: 3 axes, Hellinger + Bray-Curtis ---
cat("  Running Procrustes analysis (3 axes, Hellinger + Bray-Curtis)...\n")
ko_hel_p      <- decostand(ko_tr_filt_p, method = "hellinger")
ko_bray_p     <- vegdist(ko_hel_p, method = "bray")
ko_ord3_p     <- cmdscale(ko_bray_p, k = 3, eig = TRUE)

plant_s3_p    <- scores(pca_traits_p, choices = 1:3, display = "sites")
ko_s3_p       <- ko_ord3_p$points[, 1:3]

proc3_p       <- procrustes(plant_s3_p, ko_s3_p, symmetric = TRUE)
prot3_p       <- protest(plant_s3_p, ko_s3_p, permutations = 9999)

r_3ax   <- sqrt(1 - proc3_p$ss)
m2_3ax  <- proc3_p$ss
p_3ax   <- prot3_p$signif
cat(sprintf("  r = %.3f, m\u00B2 = %.4f, p = %.4f\n\n", r_3ax, m2_3ax, p_3ax))

# --- Save Procrustes superimposition plot (standalone) ---
plot_data_proc <- data.frame(
  sample_id = rownames(proc3_p$X),
  treatment  = plant_matched_p$Treatment,
  trait_PC1  = proc3_p$X[, 1],
  trait_PC2  = proc3_p$X[, 2],
  ko_PC1     = proc3_p$Yrot[, 1],
  ko_PC2     = proc3_p$Yrot[, 2]
)

p_proc_main <- ggplot(plot_data_proc) +
  geom_segment(aes(x = trait_PC1, y = trait_PC2, xend = ko_PC1, yend = ko_PC2),
               color = "gray60", alpha = 0.6, linewidth = 0.5) +
  geom_point(aes(x = trait_PC1, y = trait_PC2, color = treatment),
             size = 4, shape = 16, alpha = 0.9) +
  geom_point(aes(x = ko_PC1, y = ko_PC2, color = treatment),
             size = 4, shape = 17, alpha = 0.9) +
  scale_color_manual(values = site_colors_raw, name = "Treatment",
                     labels = treatment_labels) +
  labs(
    title    = "Procrustes Superimposition: Plant Traits vs Nodule Microbiome Function",
    subtitle = sprintf("Hellinger + Bray-Curtis (3 axes): r = %.3f, m\u00B2 = %.3f, p < 0.001",
                       r_3ax, m2_3ax),
    x        = "Dimension 1",
    y        = "Dimension 2",
    caption  = paste0("Circles = plant traits PCA; Triangles = KO composition PCoA (rotated)\n",
                      "Arrows connect matched samples; shorter arrows = better fit")
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", size = 14),
    plot.subtitle   = element_text(size = 11),
    plot.caption    = element_text(hjust = 0, size = 9, color = "gray40")
  )

pdf("procrustes_main_result.pdf", width = 10, height = 8)
print(p_proc_main)
dev.off()
png("procrustes_main_result.png", width = 10, height = 8, units = "in", res = 300)
print(p_proc_main)
dev.off()
cat("  Saved: procrustes_main_result.pdf / .png\n\n")

# --- PANEL D: Plant Traits PCA (live ggplot object) ---
cat("  Building Panel D: Plant Traits PCA...\n")

plant_scores_p <- scores(pca_traits_p, choices = 1:2, display = "sites")

perm_plant_p   <- adonis2(dist(plant_std_p) ~ Treatment,
                           data = as.data.frame(plant_matched_p),
                           permutations = 9999)
perm_plant_r2  <- round(perm_plant_p$R2[1], 3)
perm_plant_pv  <- perm_plant_p$`Pr(>F)`[1]
perm_plant_lbl <- ifelse(perm_plant_pv < 0.001, "p < 0.001",
                         sprintf("p = %.3f", perm_plant_pv))
cat(sprintf("  Plant traits PERMANOVA: R\u00B2 = %.3f, %s\n",
            perm_plant_r2, perm_plant_lbl))

plot_data_plant <- data.frame(
  PC1       = plant_scores_p[, 1],
  PC2       = plant_scores_p[, 2],
  treatment = factor(plant_matched_p$Treatment, levels = treatment_order)
)

plant_hulls <- plot_data_plant %>%
  group_by(treatment) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup()

p_plant <- ggplot(plot_data_plant, aes(x = PC1, y = PC2, color = treatment,
                                        shape = treatment)) +
  geom_polygon(data = plant_hulls, aes(fill = treatment), alpha = 0.15,
               show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = site_colors_raw, name = "Treatment",
                     labels = treatment_labels) +
  scale_fill_manual(values  = site_colors_raw, guide = "none") +
  scale_shape_manual(values = site_shapes_raw, name = "Treatment",
                     labels = treatment_labels) +
  labs(
    title    = "Plant Traits PCA",
    subtitle = sprintf("PERMANOVA: R\u00B2 = %.3f, %s", perm_plant_r2, perm_plant_lbl),
    x        = sprintf("PC1 (%.1f%% variance)", plant_var_p[2, 1] * 100),
    y        = sprintf("PC2 (%.1f%% variance)", plant_var_p[2, 2] * 100),
    tag      = "D"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",          # legend removed — shared with Panel E
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  tag_theme

# Save Panel D as vector to Final_Figures
ggsave(file.path(FINAL_FIG1, "Figure_1D_PlantPCA_Panel.svg"),
       p_plant + theme(legend.position = "right"), width = 8, height = 6)
ggsave(file.path(FINAL_FIG1, "Figure_1D_PlantPCA_Panel.pdf"),
       p_plant + theme(legend.position = "right"), width = 8, height = 6)
cat("  Saved: Figure_1D_PlantPCA_Panel.svg/pdf\n")

# --- PANEL E: KO Composition PCoA (live ggplot object) ---
cat("  Building Panel E: KO Composition PCoA...\n")

# 2-axis PCoA for plotting (ko_bray_p already computed above)
ko_pcoa2_p <- cmdscale(ko_bray_p, k = 2, eig = TRUE)
bray_var_p <- ko_pcoa2_p$eig / sum(ko_pcoa2_p$eig) * 100

perm_bray_p   <- adonis2(ko_bray_p ~ Treatment,
                          data = as.data.frame(plant_matched_p),
                          permutations = 9999)
perm_bray_r2  <- round(perm_bray_p$R2[1], 3)
perm_bray_pv  <- perm_bray_p$`Pr(>F)`[1]
perm_bray_lbl <- ifelse(perm_bray_pv < 0.001, "p < 0.001",
                        sprintf("p = %.3f", perm_bray_pv))
cat(sprintf("  KO Bray-Curtis PERMANOVA: R\u00B2 = %.3f, %s\n\n",
            perm_bray_r2, perm_bray_lbl))

plot_data_bray <- data.frame(
  PCo1      = ko_pcoa2_p$points[, 1],
  PCo2      = ko_pcoa2_p$points[, 2],
  treatment = factor(plant_matched_p$Treatment, levels = treatment_order)
)

bray_hulls <- plot_data_bray %>%
  group_by(treatment) %>%
  slice(chull(PCo1, PCo2)) %>%
  ungroup()

p_bray <- ggplot(plot_data_bray, aes(x = PCo1, y = PCo2, color = treatment,
                                      shape = treatment)) +
  geom_polygon(data = bray_hulls, aes(fill = treatment), alpha = 0.15,
               show.legend = FALSE) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = site_colors_raw, name = "Treatment",
                     labels = treatment_labels) +
  scale_fill_manual(values  = site_colors_raw, guide = "none") +
  scale_shape_manual(values = site_shapes_raw, name = "Treatment",
                     labels = treatment_labels) +
  labs(
    title    = "KO Composition PCoA",
    subtitle = sprintf("Hellinger + Bray-Curtis | PERMANOVA: R\u00B2 = %.3f, %s",
                       perm_bray_r2, perm_bray_lbl),
    x        = sprintf("PCo1 (%.1f%% variance)", bray_var_p[1]),
    y        = sprintf("PCo2 (%.1f%% variance)", bray_var_p[2]),
    tag      = "E"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",   # legend extracted into shared panel below
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle    = element_text(size = 11, hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  tag_theme

# Output directory for Figure 2 panels (vector files)
FINAL_FIG2 <- "G:/My Drive/Moshe/Efrat_Guy_Project/Final_Figures/Figure_2"
dir.create(FINAL_FIG2, showWarnings = FALSE, recursive = TRUE)

# Save Panel E standalone — used as a panel in remake_taxonomy_barplot_all_levels.R
ggsave("ko_pcoa_panel_e.png",
       plot  = p_bray + theme(legend.position = "right"),
       width = 6, height = 5, dpi = 300)
ggsave("ko_pcoa_panel_e.svg",
       plot  = p_bray + theme(legend.position = "right"),
       width = 6, height = 5)
ggsave(file.path(FINAL_FIG2, "Figure_2A_KO_PCoA_Panel.svg"),
       plot  = p_bray + theme(legend.position = "right"),
       width = 6, height = 5)
ggsave(file.path(FINAL_FIG2, "Figure_2A_KO_PCoA_Panel.pdf"),
       plot  = p_bray + theme(legend.position = "right"),
       width = 6, height = 5)
cat("  Saved: ko_pcoa_panel_e.png/svg\n")
cat("  Saved: Figure_2A_KO_PCoA_Panel.svg/pdf\n")

# --- Extract shared D/E legend (same color scale) into its own grob ---
p_legend_src <- ggplot(plot_data_bray, aes(x = PCo1, y = PCo2, color = treatment,
                                            shape = treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = site_colors_raw, name = "Treatment",
                     labels = treatment_labels) +
  scale_shape_manual(values = site_shapes_raw, name = "Treatment",
                     labels = treatment_labels) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 13),
    legend.text     = element_text(size = 12)
  )

shared_legend <- {
  tmp     <- ggplot_gtable(ggplot_build(p_legend_src))
  leg_idx <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg_idx]]
}

# --- Save standalone combined ordinations (side output) ---
p_ord_combined <- p_plant + p_bray +
  plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "Ordination Analysis: Plant Traits and KO Composition",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

png("original_ordinations_combined.png", width = 16, height = 5, units = "in", res = 300)
print(p_ord_combined)
dev.off()
pdf("original_ordinations_combined.pdf", width = 16, height = 5)
print(p_ord_combined)
dev.off()
cat("  Saved: original_ordinations_combined.png / .pdf\n\n")

# ══ SECTION 5: PANEL A + FINAL ASSEMBLY ═══════════════════════════════════════
cat("================================================================================\n")
cat("SECTION 5: PANEL A + FINAL ASSEMBLY\n")
cat("================================================================================\n\n")

# --- Panel A: Inoculation success bar chart ---
data_inoc <- read_excel("Inoculation_Percent.xlsx")

treatment_display_map <- c(
  "carK"  = "CAR-C",
  "carR"  = "CAR-G",
  "mtz"   = "MTZ",
  "DGN"   = "DGN",
  "RH"    = "RH",
  "YATIR" = "YATIR",
  "hok"   = "HUK",
  "ces"   = "CES"
)

inoculation_summary <- data_inoc %>%
  rename(percent_inoculated = Inoculation_Sucess) %>%
  dplyr::select(treatment, percent_inoculated) %>%
  filter(!is.na(treatment)) %>%
  mutate(
    treatment_display = dplyr::recode(treatment, !!!treatment_display_map),
    treatment_display = factor(treatment_display, levels = sort(unique(treatment_display)))
  )

panel_a <- ggplot(inoculation_summary,
                  aes(x = treatment_display,
                      y = percent_inoculated,
                      fill = treatment_display)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", percent_inoculated)),
            vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = site_colors) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  labs(
    y     = "Nodulation success rate (%)",
    title = "Inoculation Success",
    tag   = "A"
  ) +
  theme_bw() +
  theme(
    axis.title.x       = element_blank(),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1,
                                      face = "bold", size = 12),
    axis.text.y        = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title.y       = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 15),
    panel.grid.major.x = element_blank(),
    plot.margin        = margin(10, 5, 2, 5)
  ) +
  tag_theme

cat("  Panel A created\n")

# Save Panel A as vector to Final_Figures
ggsave(file.path(FINAL_FIG1, "Figure_1A_Nodulation_Panel.svg"), panel_a, width = 6, height = 6)
ggsave(file.path(FINAL_FIG1, "Figure_1A_Nodulation_Panel.pdf"), panel_a, width = 6, height = 6)
cat("  Saved: Figure_1A_Nodulation_Panel.svg/pdf\n")

# --- 5-panel layout assembly ---
cat("  Assembling 5-panel layout...\n\n")

# Layout (8 columns):
#   Row 1 (shorter):    A (cols 1-3)    | B (cols 4-8)
#   Rows 2-3 (taller): C (cols 1-4)    | D (cols 5-7) | Legend (col 8)
# Panel E has been moved to remake_taxonomy_barplot_all_levels.R.
layout_matrix <- rbind(c(1, 1, 1, 2, 2, 2, 2, 2),
                       c(3, 3, 3, 3, 4, 4, 4, 5),
                       c(3, 3, 3, 3, 4, 4, 4, 5))

grob_a <- ggplotGrob(panel_a)
grob_b <- panel_b_grob           # direct grob — no PNG round-trip
grob_c <- ggplotGrob(panel_c)
grob_d <- ggplotGrob(p_plant)   # live ggplot — no PNG round-trip, no legend

final_grob <- gridExtra::arrangeGrob(
  grob_a, grob_b, grob_c, grob_d, shared_legend,
  layout_matrix = layout_matrix,
  heights       = c(1, 1.3, 1.3),      # bottom rows taller for ordinations
  widths        = c(1, 1, 1, 1, 1, 1, 1, 1)   # 8 equal columns
)

cat("  Panels combined\n")

# --- Save final figure ---
png("combined_four_panel_egg.png", width = 16, height = 12,
    units = "in", res = 300, bg = "white")
grid.draw(final_grob)
dev.off()
cat("  Saved: combined_four_panel_egg.png\n")

tryCatch({
  pdf("combined_four_panel_egg.pdf", width = 16, height = 12, bg = "white")
  grid.draw(final_grob)
  dev.off()
  cat("  Saved: combined_four_panel_egg.pdf\n\n")
}, error = function(e) {
  message("  NOTE: could not write combined_four_panel_egg.pdf (file may be open). ",
          "Close the PDF viewer and re-run.")
})

cat("================================================================================\n")
cat("FIGURE CREATION COMPLETE\n")
cat("================================================================================\n\n")
cat("Output files:\n")
cat("  plant_traits_final_formatted.png  (Panel B standalone)\n")
cat("  mantel_linkET_figure.png          (Panel C standalone)\n")
cat("  original_ordinations_combined.png (Panels D+E standalone)\n")
cat("  procrustes_main_result.png        (Procrustes superimposition)\n")
cat("  combined_four_panel_egg.pdf       (MAIN FIGURE — vector)\n")
cat("  combined_four_panel_egg.png       (MAIN FIGURE — raster)\n")
cat("  trait_correlation_matrix.csv\n")
cat("  mantel_test_results.csv\n\n")
cat("Done!\n")
