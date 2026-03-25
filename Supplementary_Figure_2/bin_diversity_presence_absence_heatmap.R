################################################################################
# PRESENCE / ABSENCE HEATMAP — BORUTA ALPHA-DIVERSITY KOs
# Three separate plots, one per diversity metric (richness / shannon / simpson)
#
# Layout matches nitrogen_metabolism_symbiont_vs_endophyte.svg:
#   rows   = bins (44), split Symbiont | Endophyte, clustered within groups
#   cols   = KOs confirmed for that metric, split by Functional_Category
#   fill   = white (absent) / #3C5488 dark blue (present)
#   left   = rowAnnotation: Organism (red/cyan) + Treatment
#   top    = HeatmapAnnotation: Functional_Category + Correlation Direction
#   bottom = prevalence barplots (Symbiont, Endophyte)
#
# Output files (per metric):
#   bin_diversity_presence_absence_heatmap_<metric>.svg / .pdf
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(svglite)
})

cat("================================================================================\n")
cat("BORUTA ALPHA-DIVERSITY KO  —  PRESENCE / ABSENCE HEATMAP  (per metric)\n")
cat("================================================================================\n\n")

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------

BASE_DIR   <- "/mnt/c/Users/owner/My Drive (moshe.alon@mail.huji.ac.il)/Moshe/Efrat_Guy_Project"
OUTPUT_DIR <- file.path(BASE_DIR, "Boruta_New/bin_diversity_boruta")

boruta_csv    <- file.path(OUTPUT_DIR, "bin_diversity_KO_all_metrics_annotated.csv")
ko_table_xlsx <- file.path(BASE_DIR, "Boruta_New/ko_presence_absence_table.xlsx")
abundance_tsv <- file.path(BASE_DIR, "New_Binning_Pipeline_Coassembly/combined_abundance_long_renormalized.tsv")

stopifnot("Boruta CSV not found"               = file.exists(boruta_csv))
stopifnot("KO presence/absence table not found" = file.exists(ko_table_xlsx))
stopifnot("Abundance TSV not found"             = file.exists(abundance_tsv))

################################################################################
# STEP 1 — Load data
################################################################################

cat("Step 1: Loading data...\n")

boruta_df <- read_csv(boruta_csv, show_col_types = FALSE)
ko_pa     <- read_excel(ko_table_xlsx)
abundance <- read_tsv(abundance_tsv, show_col_types = FALSE)

cat(sprintf("  Boruta CSV: %d rows | metrics: %s\n",
            nrow(boruta_df),
            paste(sort(unique(boruta_df$Metric)), collapse = ", ")))
cat(sprintf("  KO table:   %d KOs × %d bins\n", nrow(ko_pa), ncol(ko_pa) - 1L))

################################################################################
# STEP 2 — Bin metadata  (shared across all three plots)
################################################################################

cat("\nStep 2: Building bin metadata...\n")

extract_readable_name <- function(classification) {
  if (is.na(classification) || grepl("^Unclassified", trimws(classification)))
    return(trimws(classification))

  species <- str_extract(classification, "s__([^;]+)$")
  species <- gsub("s__", "", species)
  genus   <- str_extract(classification, "g__([^;]+);")
  if (is.na(genus)) genus <- str_extract(classification, "g__([^;]+)$")
  genus <- gsub("g__|;", "", genus)

  if (!is.na(species) && nchar(trimws(species)) > 0) {
    sc <- gsub("_", " ", species)
    if (grepl("^[A-Z][a-z]+\\s", sc)) return(sc)
    if (!is.na(genus) && nchar(trimws(genus)) > 0) return(paste(genus, sc))
    return(sc)
  }
  if (!is.na(genus) && nchar(trimws(genus)) > 0) return(genus)

  for (pat in c("f__([^;]+);", "o__([^;]+);", "c__([^;]+);",
                "p__([^;]+);", "d__([^;]+);")) {
    x <- str_extract(classification, pat)
    x <- gsub("[fgocdp]__|;", "", x)
    if (!is.na(x) && nchar(trimws(x)) > 0) return(x)
  }
  return("Unknown")
}

treatment_names <- c(
  "ces"  = "CES",
  "RH"   = "RH",
  "carR" = "CAR-G",
  "carK" = "CAR-C",
  "hok"  = "HUK",
  "mtz"  = "MTZ"
)
treatment_order <- c("MTZ", "HUK", "CES", "RH", "CAR-C", "CAR-G")

bin_meta <- abundance %>%
  distinct(Bin_Prefixed, Organism, Treatment, classification) %>%
  mutate(
    organism_type  = Organism,
    Treatment_Full = recode(Treatment, !!!treatment_names),
    Treatment_Full = factor(Treatment_Full, levels = treatment_order),
    Readable_name  = vapply(classification, extract_readable_name, character(1))
  ) %>%
  arrange(Treatment_Full,
          factor(organism_type, levels = c("Symbiont", "Endophyte")),
          Readable_name) %>%
  group_by(Treatment_Full, organism_type, Readable_name) %>%
  mutate(facet_n = n(), facet_rank = row_number()) %>%
  ungroup() %>%
  mutate(bin_label_v1 = if_else(facet_n > 1,
                                paste0(Readable_name, " (", facet_rank, ")"),
                                Readable_name)) %>%
  group_by(bin_label_v1) %>%
  mutate(global_n = n(), global_rank = row_number()) %>%
  ungroup() %>%
  mutate(bin_label = if_else(global_n > 1,
                             paste0(bin_label_v1, " [", Treatment, "]"),
                             bin_label_v1)) %>%
  select(-facet_n, -facet_rank, -bin_label_v1, -global_n, -global_rank)

# Fixed row order: Symbiont first, then Endophyte; within each group by treatment
row_order_all <- bin_meta %>%
  arrange(factor(organism_type, levels = c("Symbiont", "Endophyte")),
          Treatment_Full, bin_label) %>%
  pull(Bin_Prefixed)

cat(sprintf("  Unique bins: %d\n", nrow(bin_meta)))

################################################################################
# STEP 3 — Functional categories  (shared)
################################################################################

cat("\nStep 3: Assigning functional categories...\n")

brite_to_func <- c(
  "Energy metabolism"                                  = "Energy & Carbon metabolism",
  "Carbohydrate metabolism"                            = "Energy & Carbon metabolism",
  "Metabolism of cofactors and vitamins"               = "Cofactor & Vitamin metabolism",
  "Amino acid metabolism"                              = "Amino acid metabolism",
  "Xenobiotics biodegradation and metabolism"          = "Xenobiotics degradation",
  "Lipid metabolism"                                   = "Lipid metabolism",
  "Protein families: signaling and cellular processes" = "Signaling & Cell regulation",
  "Cellular community - prokaryotes"                   = "Signaling & Cell regulation",
  "Quorum sensing"                                     = "Signaling & Cell regulation",
  "Nucleotide metabolism"                              = "Nucleotide metabolism",
  "Transport and catabolism"                           = "Transport",
  "Not Included in Pathway or Brite"                   = "Unknown / Poorly characterized",
  "Poorly characterized"                               = "Unknown / Poorly characterized"
)

func_order <- c(
  "Cofactor & Vitamin metabolism",
  "Energy & Carbon metabolism",
  "Xenobiotics degradation",
  "Amino acid metabolism",
  "Signaling & Cell regulation",
  "Lipid metabolism",
  "Transport",
  "Nucleotide metabolism",
  "Other metabolism",
  "Unknown / Poorly characterized"
)

# Add Functional_Category to every row of boruta_df (keeps all metric rows)
boruta_df <- boruta_df %>%
  mutate(
    BRITE_L2_first      = trimws(sub(";.*", "", BRITE_L2)),
    Functional_Category = if_else(
      BRITE_L2_first %in% names(brite_to_func),
      brite_to_func[BRITE_L2_first],
      "Other metabolism"
    ),
    Functional_Category = factor(Functional_Category, levels = func_order)
  )

################################################################################
# STEP 4 — Full presence/absence matrix  (subset per metric later)
################################################################################

cat("\nStep 4: Building full presence/absence matrix...\n")

all_kos       <- unique(boruta_df$KO)
kos_available <- intersect(all_kos, ko_pa$KO)
kos_missing   <- setdiff(all_kos, ko_pa$KO)
if (length(kos_missing) > 0)
  cat(sprintf("  NOTE: %d KOs absent from table (skipped): %s\n",
              length(kos_missing), paste(kos_missing, collapse = ", ")))

pa_wide <- ko_pa %>%
  filter(KO %in% kos_available) %>%
  column_to_rownames("KO") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Bin_Prefixed") %>%
  semi_join(bin_meta, by = "Bin_Prefixed")

pa_matrix_full <- as.matrix(pa_wide[, -1])
rownames(pa_matrix_full) <- pa_wide$Bin_Prefixed
storage.mode(pa_matrix_full) <- "numeric"
pa_matrix_full[pa_matrix_full > 1]   <- 1
pa_matrix_full[is.na(pa_matrix_full)] <- 0

# Apply fixed row order
row_order <- intersect(row_order_all, rownames(pa_matrix_full))
pa_matrix_full <- pa_matrix_full[row_order, ]

# Aligned bin annotation (fixed for all three plots)
bin_annot <- bin_meta[match(rownames(pa_matrix_full), bin_meta$Bin_Prefixed), ]
is_sym    <- bin_annot$organism_type == "Symbiont"
is_end    <- bin_annot$organism_type == "Endophyte"
n_sym     <- sum(is_sym)
n_end     <- sum(is_end)

cat(sprintf("  Full matrix: %d bins × %d KOs\n",
            nrow(pa_matrix_full), ncol(pa_matrix_full)))

################################################################################
# STEP 5 — Shared color palettes
################################################################################

organism_colors <- c("Symbiont" = "#E64B35", "Endophyte" = "#4DBBD5")

treatment_colors <- c(
  "MTZ"   = "#FFD92F",
  "HUK"   = "#A6D854",
  "CES"   = "#8DA0CB",
  "RH"    = "#E5C494",
  "CAR-C" = "#66C2A5",
  "CAR-G" = "#FC8D62"
)

func_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(length(func_order)),
  func_order
)

direction_colors <- c("Positive" = "#2ca25f", "Negative" = "#de2d26")

# Fixed left row annotation (same for all three plots)
ha_row <- rowAnnotation(
  Organism  = bin_annot$organism_type,
  Treatment = as.character(bin_annot$Treatment_Full),
  col = list(Organism = organism_colors, Treatment = treatment_colors),
  annotation_name_gp = gpar(fontsize = 9),
  annotation_legend_param = list(
    Organism  = list(title = "Type",
                     title_gp = gpar(fontsize = 10, fontface = "bold"),
                     labels_gp = gpar(fontsize = 9)),
    Treatment = list(title = "Treatment",
                     title_gp = gpar(fontsize = 10, fontface = "bold"),
                     labels_gp = gpar(fontsize = 9))
  ),
  width = unit(8, "mm")
)

################################################################################
# STEP 6 — Per-metric heatmap loop
################################################################################

metric_titles <- c(
  richness = "Species Richness",
  shannon  = "Shannon Diversity",
  simpson  = "Simpson Diversity"
)

for (metric_name in c("richness", "shannon", "simpson")) {

  cat(sprintf("\n--- Metric: %s ---\n", metric_name))

  # KOs confirmed for this metric (from boruta_df, exactly this metric's rows)
  metric_df <- boruta_df %>%
    filter(Metric == metric_name, KO %in% kos_available) %>%
    arrange(Functional_Category, desc(MeanImp))

  n_kos_metric <- nrow(metric_df)
  cat(sprintf("  KOs: %d\n", n_kos_metric))

  if (n_kos_metric == 0) {
    cat("  Skipping — no KOs available.\n")
    next
  }

  # Subset matrix to this metric's KOs (in category + MeanImp order)
  col_order_m  <- metric_df$KO
  pa_m         <- pa_matrix_full[, col_order_m, drop = FALSE]
  ko_annot_m   <- metric_df   # already ordered

  # Prevalence for this KO subset
  sym_prev_m <- if (n_sym > 0)
    colSums(pa_m[is_sym, , drop = FALSE]) / n_sym * 100
  else rep(0, ncol(pa_m))

  end_prev_m <- if (n_end > 0)
    colSums(pa_m[is_end, , drop = FALSE]) / n_end * 100
  else rep(0, ncol(pa_m))

  # ---- CSV export ----
  csv_df <- metric_df %>%
    select(KO, Metric, MeanImp, SpearmanRho, Pvalue,
           Description, Pathways, Modules,
           BRITE_L1, BRITE_L2, BRITE_L3,
           in_symbiont, symbiont_label,
           Functional_Category) %>%
    mutate(
      Symbiont_prevalence_pct  = round(sym_prev_m[KO], 2),
      Endophyte_prevalence_pct = round(end_prev_m[KO], 2),
      Total_prevalence_pct     = round(
        colSums(pa_m) / nrow(pa_m) * 100, 2)[KO],
      Correlation_direction    = if_else(SpearmanRho >= 0, "Positive", "Negative")
    ) %>%
    arrange(Functional_Category, desc(MeanImp))

  csv_out <- file.path(OUTPUT_DIR,
                       sprintf("bin_diversity_KO_%s_annotated_with_prevalence.csv",
                               metric_name))
  write_csv(csv_df, csv_out)
  cat(sprintf("  Saved CSV: %s\n", basename(csv_out)))

  # ---- Top annotation (Category + Direction) ----
  ha_top_m <- HeatmapAnnotation(
    Category  = as.character(ko_annot_m$Functional_Category),
    Direction = if_else(ko_annot_m$SpearmanRho >= 0, "Positive", "Negative"),
    col = list(Category = func_colors, Direction = direction_colors),
    annotation_name_gp = gpar(fontsize = 9),
    annotation_legend_param = list(
      Category  = list(title = "Functional\nCategory",
                       title_gp = gpar(fontsize = 10, fontface = "bold"),
                       labels_gp = gpar(fontsize = 8)),
      Direction = list(title = "Correlation\nDirection",
                       title_gp = gpar(fontsize = 10, fontface = "bold"),
                       labels_gp = gpar(fontsize = 9))
    ),
    height = unit(10, "mm")
  )

  # ---- Bottom annotation: Spearman rho barplot with p-value asterisks ----
  pval_to_star <- function(p) {
    dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  }

  rho_vec    <- ko_annot_m$SpearmanRho
  pval_vec   <- ko_annot_m$Pvalue
  star_vec   <- pval_to_star(pval_vec)
  bar_colors <- ifelse(rho_vec >= 0,
                       direction_colors["Positive"],
                       direction_colors["Negative"])

  ha_bottom_m <- HeatmapAnnotation(
    Stars = anno_text(
      star_vec,
      gp = gpar(fontsize = 7),
      height = unit(6, "mm")
    ),
    `Spearman rho` = anno_barplot(
      rho_vec,
      gp = gpar(fill = bar_colors, col = bar_colors),
      baseline = 0,
      ylim = c(min(c(rho_vec, -0.05)), max(c(rho_vec, 0.05))),
      height = unit(3, "cm"),
      axis_param = list(side = "left", gp = gpar(fontsize = 8))
    ),
    annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
    annotation_name_side = "left",
    gap = unit(1, "mm")
  )

  # ---- Column and row labels ----
  col_labels_m <- paste0(
    str_trunc(sub(";.*", "", ko_annot_m$Description), 22),
    "\n", ko_annot_m$KO
  )
  row_labels_m <- bin_annot$bin_label

  # ---- Main heatmap ----
  ht_m <- Heatmap(
    pa_m,
    name = "Gene\nPresence",

    col = c("0" = "white", "1" = "#3C5488"),

    # Rows ordered by treatment within Symbiont / Endophyte groups
    cluster_rows  = FALSE,
    row_split     = bin_annot$organism_type,

    # Column ordering — no clustering, respect category order
    cluster_columns = FALSE,
    column_split    = as.character(ko_annot_m$Functional_Category),
    column_title    = NULL,

    row_labels            = row_labels_m,
    row_names_gp          = gpar(fontsize = 8),
    column_labels         = col_labels_m,
    column_names_gp       = gpar(fontsize = 7),
    column_names_rot      = 90,
    column_names_centered = FALSE,

    left_annotation   = ha_row,
    top_annotation    = ha_top_m,
    bottom_annotation = ha_bottom_m,

    border     = TRUE,
    row_gap    = unit(3, "mm"),
    column_gap = unit(2, "mm"),
    row_title  = NULL,
    show_column_dend = FALSE,

    heatmap_legend_param = list(
      title     = "Gene\nPresence",
      title_gp  = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9),
      at        = c(0, 1),
      labels    = c("Absent", "Present")
    ),

    width  = unit(0.7 * ncol(pa_m), "cm"),
    height = unit(0.45 * nrow(pa_m), "cm")
  )

  # ---- Dimensions & save ----
  fig_w <- min(36, max(16, 0.7 * ncol(pa_m) + 10))
  fig_h <- min(24, max(14, 0.45 * nrow(pa_m) + 12))
  cat(sprintf("  Canvas: %.1f × %.1f inches\n", fig_w, fig_h))

  stem    <- sprintf("bin_diversity_presence_absence_heatmap_%s", metric_name)
  pdf_out <- file.path(OUTPUT_DIR, paste0(stem, ".pdf"))
  svg_out <- file.path(OUTPUT_DIR, paste0(stem, ".svg"))

  pdf(pdf_out, width = fig_w, height = fig_h)
  draw(ht_m,
       merge_legend        = TRUE,
       column_title        = sprintf("Boruta KOs - %s", metric_titles[metric_name]),
       column_title_gp     = gpar(fontsize = 13, fontface = "bold"),
       padding             = unit(c(2, 2, 2, 12), "mm"),
       heatmap_legend_side = "right")
  dev.off()
  cat(sprintf("  Saved PDF: %s\n", basename(pdf_out)))

  svglite(svg_out, width = fig_w, height = fig_h)
  draw(ht_m,
       merge_legend        = TRUE,
       column_title        = sprintf("Boruta KOs - %s", metric_titles[metric_name]),
       column_title_gp     = gpar(fontsize = 13, fontface = "bold"),
       padding             = unit(c(2, 2, 2, 12), "mm"),
       heatmap_legend_side = "right")
  dev.off()
  cat(sprintf("  Saved SVG: %s\n", basename(svg_out)))
}

cat("\n================================================================================\n")
cat("Done. Output files per metric:\n")
for (m in c("richness", "shannon", "simpson")) {
  cat(sprintf("  bin_diversity_presence_absence_heatmap_%s.svg / .pdf\n", m))
  cat(sprintf("  bin_diversity_KO_%s_annotated_with_prevalence.csv\n", m))
}
cat("================================================================================\n")
