#!/usr/bin/env Rscript
# ============================================================================
# T6SS KO ABUNDANCE vs. BIN-BASED ALPHA DIVERSITY — REGRESSION SCATTER PLOTS
# ============================================================================
# 15 T6SS core KOs regressed against richness, Shannon, Simpson.
# One scatter+regression plot per diversity metric (all genes overlaid).
# vasL (K11911) excluded — all-zero abundance.
# NOTE: No KEGG module (M-number) exists for T6SS; classification source is
#       KEGG BRITE ko02044 (Secretion system).
# ============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(readxl)

select <- dplyr::select
filter <- dplyr::filter

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR <- "/mnt/c/Users/owner/My Drive (moshe.alon@mail.huji.ac.il)/Moshe/Efrat_Guy_Project"

KO_FILE        <- file.path(BASE_DIR, "Boruta_New/normalized_kegg_results.tsv")
BIN_ABUND_FILE <- file.path(BASE_DIR,
  "New_Binning_Pipeline_Coassembly/combined_abundance_long_renormalized.tsv")
KO_PA_FILE     <- file.path(BASE_DIR, "Boruta_New/ko_presence_absence_table.xlsx")
OUTPUT_DIR     <- file.path(BASE_DIR, "Boruta_New/bin_diversity_boruta")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUTPUT_PREFIX <- file.path(OUTPUT_DIR, "T6SS_bin_diversity_regression")

# ── T6SS core KO definitions (KEGG BRITE ko02044, Imp/Vas core) ─────────────
# K11911 (vasL) excluded: all-zero abundance across all samples.
T6SS_KOS <- data.frame(
  KO = c("K11891","K11892","K11893","K11895","K11896","K11897",
          "K11900","K11901","K11902","K11903","K11904","K11905",
          "K11906","K11907","K11910"),
  gene_label = c(
    "impL/vasK (K11891)",
    "impK/vasF (K11892)",
    "impJ/vasE (K11893)",
    "impH/vasB (K11895)",
    "impG/vasA (K11896)",
    "impF (K11897)",
    "impC/tssC (K11900)",
    "impB (K11901)",
    "impA (K11902)",
    "hcp/tssB (K11903)",
    "vgrG/tssA (K11904)",
    "tssE (K11905)",
    "vasD/lip (K11906)",
    "vasG/clpV (K11907)",
    "vasJ (K11910)"
  ),
  stringsAsFactors = FALSE
)

# ── T6SS "related" KOs present in dataset but NOT in canonical core ──────────
# Classified as "Imp/Vas related proteins" in BRITE ko02044.
T6SS_RELATED <- data.frame(
  KO = c("K11890","K11894","K11898","K11899","K11908","K11909"),
  gene_label = c(
    "impM (K11890)",
    "impI/vasC (K11894)",
    "impE (K11898)",
    "impD (K11899)",
    "vasH (K11908)",
    "vasI (K11909)"
  ),
  stringsAsFactors = FALSE
)

cat("\n")
cat("========================================================================\n")
cat("T6SS KO ABUNDANCE vs. BIN-BASED ALPHA DIVERSITY — REGRESSION\n")
cat("========================================================================\n\n")
cat("T6SS KOs to analyse:", nrow(T6SS_KOS), "(vasL/K11911 excluded: all-zero)\n")
cat("NOTE: No KEGG module (M-number) exists for T6SS.\n")
cat("      Classification source: KEGG BRITE ko02044 (Secretion system).\n")
cat("      Related KOs in dataset (NOT in core):", paste(T6SS_RELATED$KO, collapse=", "), "\n\n")

# ============================================================================
# PHASE 1: COMPUTE BIN DIVERSITY PER SAMPLE
# ============================================================================

cat("=== PHASE 1: Compute Bin Alpha Diversity ===\n")

bin_abund <- read_tsv(BIN_ABUND_FILE, show_col_types = FALSE)
cat("  Abundance file loaded:", nrow(bin_abund), "rows\n")

ko_pa <- read_excel(KO_PA_FILE)
cat("  KO presence/absence table loaded:", nrow(ko_pa), "KOs x",
    ncol(ko_pa) - 1, "bins\n")

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
  rename(sample_id = Sample)

cat("  Diversity computed for", nrow(bin_diversity), "samples\n")

# ============================================================================
# PHASE 2: LOAD KO DATA AND MERGE TECHNICAL REPLICATES
# ============================================================================

cat("\n=== PHASE 2: Load KO Data & Merge Technical Replicates ===\n")

ko_data <- read_tsv(KO_FILE, show_col_types = FALSE)
cat("  KO file loaded:", nrow(ko_data), "rows x", ncol(ko_data) - 1, "raw columns\n")

sample_cols <- colnames(ko_data)[-1]

extract_sample_id <- function(col_name) {
  clean <- sub("_kegg_hits\\.tsv$", "", col_name)
  parts <- strsplit(clean, "_")[[1]]
  paste(parts[1], parts[2], sep = "_")
}

sample_info <- data.frame(
  original_name = sample_cols,
  sample_id     = sapply(sample_cols, extract_sample_id),
  stringsAsFactors = FALSE
)

ko_matrix <- as.matrix(ko_data[, -1])
rownames(ko_matrix) <- ko_data$kegg_number

unique_samples <- unique(sample_info$sample_id)
merged_matrix  <- matrix(
  NA,
  nrow = nrow(ko_matrix),
  ncol = length(unique_samples),
  dimnames = list(rownames(ko_matrix), unique_samples)
)
for (sid in unique_samples) {
  cols <- sample_info$original_name[sample_info$sample_id == sid]
  merged_matrix[, sid] <- if (length(cols) == 1) {
    ko_matrix[, cols]
  } else {
    rowMeans(ko_matrix[, cols], na.rm = TRUE)
  }
}
cat("  Merged to", ncol(merged_matrix), "unique samples\n")

# ============================================================================
# PHASE 3: FILTER TO T6SS KOs
# ============================================================================

cat("\n=== PHASE 3: Filter to T6SS KOs ===\n")

present_kos <- T6SS_KOS$KO[T6SS_KOS$KO %in% rownames(merged_matrix)]
missing_kos <- T6SS_KOS$KO[!T6SS_KOS$KO %in% rownames(merged_matrix)]

cat("  Requested:", nrow(T6SS_KOS), "core KOs\n")
cat("  Present in data:", length(present_kos), "\n")
if (length(missing_kos) > 0)
  cat("  MISSING from data:", paste(missing_kos, collapse = ", "), "\n")

# Build long-format table of which bins carry each T6SS KO
t6ss_kos_present <- present_kos   # alias used in carrying-fraction steps below
t6ss_ko_pa_long <- ko_pa %>%
  filter(KO %in% t6ss_kos_present) %>%
  pivot_longer(-KO, names_to = "Bin_Prefixed", values_to = "has_gene") %>%
  filter(has_gene == 1) %>%
  select(KO, Bin_Prefixed)
cat("  Carrying-bin table: ", nrow(t6ss_ko_pa_long),
    "KO × bin pairs across", length(unique(t6ss_ko_pa_long$KO)), "KOs\n")

# Compute carrying-bin fraction (0–100 scale) per KO × sample
carrying_fraction <- t6ss_ko_pa_long %>%
  left_join(
    bin_abund %>% select(Bin_Prefixed, Sample, Relative_Abundance),
    by = "Bin_Prefixed"
  ) %>%
  rename(sample_id = Sample) %>%
  group_by(KO, sample_id) %>%
  summarise(carrying_pct = sum(Relative_Abundance, na.rm = TRUE), .groups = "drop")
cat("  Carrying-fraction table:", nrow(carrying_fraction), "rows\n")

related_present <- T6SS_RELATED$KO[T6SS_RELATED$KO %in% rownames(merged_matrix)]
related_missing <- T6SS_RELATED$KO[!T6SS_RELATED$KO %in% rownames(merged_matrix)]
cat("  Related KOs present in dataset (excluded):", paste(related_present, collapse = ", "), "\n")
if (length(related_missing) > 0)
  cat("  Related KOs absent from dataset:", paste(related_missing, collapse = ", "), "\n")

t6ss_matrix <- merged_matrix[present_kos, , drop = FALSE]

# ============================================================================
# PHASE 4: ALIGN SAMPLES
# ============================================================================

cat("\n=== PHASE 4: Align Samples ===\n")

common_samples <- intersect(colnames(t6ss_matrix), bin_diversity$sample_id)
cat("  Samples in KO table:        ", ncol(t6ss_matrix), "\n")
cat("  Samples in diversity table: ", nrow(bin_diversity), "\n")
cat("  Samples in common:          ", length(common_samples), "\n")

if (length(common_samples) == 0)
  stop("No matching sample IDs between KO and diversity tables.")

t6ss_matched <- t6ss_matrix[, common_samples, drop = FALSE]
div_matched  <- bin_diversity %>%
  filter(sample_id %in% common_samples) %>%
  arrange(match(sample_id, common_samples))

# ============================================================================
# PHASE 5: BUILD LONG-FORMAT DATA FRAME
# ============================================================================

cat("\n=== PHASE 5: Build Long-Format Data Frame ===\n")

t6ss_long <- as.data.frame(t(t6ss_matched)) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(
    cols      = -sample_id,
    names_to  = "KO",
    values_to = "abundance"
  ) %>%
  left_join(T6SS_KOS, by = "KO") %>%
  left_join(div_matched, by = "sample_id") %>%
  mutate(treatment = sub("_\\d+$", "", sample_id))

# Re-weight abundance: multiply by fraction of community that carries the gene
t6ss_long <- t6ss_long %>%
  left_join(carrying_fraction, by = c("KO", "sample_id")) %>%
  mutate(
    carrying_pct = replace_na(carrying_pct, 0),
    abundance    = abundance * (carrying_pct / 100)
  ) %>%
  select(-carrying_pct)

cat("  Long-format rows:", nrow(t6ss_long), "\n")
cat("  Unique samples:", length(unique(t6ss_long$sample_id)), "\n")
cat("  Unique KOs:", length(unique(t6ss_long$KO)), "\n")

# ============================================================================
# PHASE 6: DROP ZERO-VARIANCE KOs
# ============================================================================

cat("\n=== PHASE 6: Check Zero-Variance KOs ===\n")

ko_var <- t6ss_long %>%
  group_by(KO, gene_label) %>%
  summarise(
    max_abund = max(abundance, na.rm = TRUE),
    n_nonzero = sum(abundance > 0, na.rm = TRUE),
    .groups = "drop"
  )

zero_var_kos <- ko_var %>% filter(max_abund == 0)
valid_kos    <- ko_var %>% filter(max_abund > 0)

if (nrow(zero_var_kos) > 0) {
  cat("  WARNING: All-zero KOs (excluded from plots):\n")
  for (i in seq_len(nrow(zero_var_kos)))
    cat("    -", zero_var_kos$KO[i], "(", zero_var_kos$gene_label[i], ")\n")
} else {
  cat("  All", nrow(ko_var), "KOs have non-zero values — none excluded.\n")
}

t6ss_plot <- t6ss_long %>% filter(KO %in% valid_kos$KO)
n_genes   <- length(unique(t6ss_plot$KO))
cat("  KOs retained for plotting:", n_genes, "\n")

# ============================================================================
# PHASE 7: COLOUR PALETTE
# ============================================================================

cat("\n=== PHASE 7: Build Colour Palette ===\n")

gene_levels  <- unique(t6ss_plot$gene_label)
n_genes_plot <- length(gene_levels)

t6ss_colors <- setNames(
  colorRampPalette(brewer.pal(12, "Paired"))(n_genes_plot),
  gene_levels
)
cat("  Colours assigned:", n_genes_plot, "\n")

# ============================================================================
# PHASE 8: SPEARMAN CORRELATIONS
# ============================================================================

cat("\n=== PHASE 8: Spearman Correlations ===\n")

metrics     <- c("richness", "shannon", "simpson")
metric_labs <- c(richness = "Bin Richness", shannon = "Shannon Index (H')",
                 simpson  = "Simpson Index (1 - D)")

spearman_all <- expand.grid(KO = valid_kos$KO, metric = metrics,
                             stringsAsFactors = FALSE) %>%
  left_join(T6SS_KOS, by = "KO") %>%
  rowwise() %>%
  mutate(
    vals_abund = list(t6ss_plot$abundance[t6ss_plot$KO == KO]),
    vals_div   = list(t6ss_plot[[metric]][t6ss_plot$KO == KO]),
    n          = length(vals_abund),
    ct         = list(tryCatch(
      cor.test(vals_abund, vals_div, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )),
    rho  = if (!is.null(ct)) unname(ct$estimate) else NA_real_,
    pval = if (!is.null(ct)) ct$p.value           else NA_real_
  ) %>%
  select(KO, gene_label, metric, n, rho, pval) %>%
  ungroup()

cat("\n  Spearman correlation table (sorted by |rho|):\n")
cat(sprintf("  %-22s %-10s %7s %10s %5s\n",
            "gene_label", "metric", "rho", "p-value", "n"))
cat("  ", paste(rep("-", 60), collapse=""), "\n")
spearman_all %>%
  arrange(desc(abs(rho))) %>%
  rowwise() %>%
  group_walk(~ {
    cat(sprintf("  %-22s %-10s %+7.4f %10.4g %5d\n",
                substr(.x$gene_label, 1, 22), .x$metric,
                .x$rho, .x$pval, .x$n))
  })

# ============================================================================
# PHASE 9: PLOTTING FUNCTION
# ============================================================================

make_regression_plot <- function(data, metric, y_label, title_label,
                                 colors, spearman_df) {

  top3 <- spearman_df %>%
    filter(metric == !!metric) %>%
    arrange(desc(abs(rho))) %>%
    slice_head(n = 3) %>%
    mutate(lab = sprintf("%s: ρ=%+.2f, p=%.3g", substr(gene_label, 1, 18), rho, pval))
  caption_txt <- paste0("Top |ρ| genes:\n", paste(top3$lab, collapse = "\n"))

  p <- ggplot(data, aes(x = abundance, y = .data[[metric]],
                        color = gene_label, fill = gene_label)) +
    geom_point(size = 2, alpha = 0.65, shape = 16) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.10, linewidth = 0.85) +
    scale_color_manual(values = colors, name = "T6SS gene") +
    scale_fill_manual( values = colors, name = "T6SS gene") +
    labs(
      title   = title_label,
      x       = "KO Abundance (carrying-bin weighted)",
      y       = y_label,
      caption = caption_txt
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.caption    = element_text(size = 7.5, color = "grey40", hjust = 0),
      legend.position = "right",
      legend.text     = element_text(size = 7.5),
      legend.key.size = unit(0.55, "lines"),
      legend.title    = element_text(size = 8.5, face = "bold")
    )

  return(p)
}

# ============================================================================
# PHASE 10: BUILD AND SAVE INDIVIDUAL PLOTS
# ============================================================================

cat("\n=== PHASE 10: Build & Save Plots ===\n")

plot_list <- list()

for (metric in metrics) {
  cat("  Building plot for:", metric, "\n")

  p <- make_regression_plot(
    data        = t6ss_plot,
    metric      = metric,
    y_label     = metric_labs[[metric]],
    title_label = paste("T6SS Genes vs.", metric_labs[[metric]]),
    colors      = t6ss_colors,
    spearman_df = spearman_all
  )

  plot_list[[metric]] <- p

  pdf(paste0(OUTPUT_PREFIX, "_", metric, ".pdf"), width = 8, height = 6)
  print(p)
  dev.off()

  png(paste0(OUTPUT_PREFIX, "_", metric, ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(p)
  dev.off()

  cat("    Saved:", metric, "(PDF + PNG)\n")
}

# ============================================================================
# PHASE 11: COMBINED 3-PANEL FIGURE
# ============================================================================

cat("\n=== PHASE 11: Combined 3-Panel Figure ===\n")

p_combined <- (plot_list$richness + plot_list$shannon + plot_list$simpson) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title    = "T6SS Gene Cluster Abundance vs. Bin Alpha Diversity",
    subtitle = paste0(
      "n = ", length(common_samples), " samples | ",
      n_genes_plot, " T6SS KOs | Spearman ρ per gene"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  ) &
  theme(legend.position = "right")

pdf(paste0(OUTPUT_PREFIX, "_combined.pdf"), width = 20, height = 6)
print(p_combined)
dev.off()

png(paste0(OUTPUT_PREFIX, "_combined.png"),
    width = 20, height = 6, units = "in", res = 300)
print(p_combined)
dev.off()

svg(paste0(OUTPUT_PREFIX, "_combined.svg"), width = 20, height = 6)
print(p_combined)
dev.off()

cat("  Saved: combined 3-panel (PDF + PNG + SVG)\n")

# ============================================================================
# FULL SPEARMAN TABLE TO CONSOLE
# ============================================================================

cat("\n=== Full Spearman Correlation Table (all", nrow(spearman_all), "pairs) ===\n")
cat(sprintf("  %-28s %-10s %7s %10s %5s\n",
            "gene_label", "metric", "rho", "p-value", "n"))
cat("  ", paste(rep("-", 66), collapse=""), "\n")
spearman_all %>%
  arrange(desc(abs(rho))) %>%
  rowwise() %>%
  group_walk(~ {
    sig <- if (!is.na(.x$pval) && .x$pval < 0.05) " *" else ""
    cat(sprintf("  %-28s %-10s %+7.4f %10.4g %5d%s\n",
                substr(.x$gene_label, 1, 28), .x$metric,
                .x$rho, .x$pval, .x$n, sig))
  })

# ============================================================================
# PHASE 12: OPERON SUBSET PLOTS
# ============================================================================
# Six KOs forming a complete operon identified within the bins.
# ============================================================================

cat("\n=== PHASE 12: Operon Subset Plots ===\n")

OPERON_KOS <- c("K11905","K11900","K11896","K11895","K11910","K11901")

t6ss_operon <- t6ss_plot %>% filter(KO %in% OPERON_KOS)
cat("  Operon KOs:", paste(OPERON_KOS, collapse = ", "), "\n")
cat("  Rows in subset:", nrow(t6ss_operon), "\n")

operon_levels <- unique(t6ss_operon$gene_label)
operon_colors <- setNames(
  brewer.pal(length(operon_levels), "Dark2"),
  operon_levels
)

spearman_operon <- spearman_all %>% filter(KO %in% OPERON_KOS)

OPERON_PREFIX <- file.path(OUTPUT_DIR, "T6SS_operon_regression")

make_operon_plot <- function(data, metric, y_label, title_label,
                             colors, spearman_df) {

  rho_labels <- spearman_df %>%
    filter(metric == !!metric) %>%
    arrange(desc(abs(rho))) %>%
    mutate(lab = sprintf("%s\nρ=%+.2f, p=%.3g",
                         substr(gene_label, 1, 20), rho, pval))

  p <- ggplot(data, aes(x = abundance, y = .data[[metric]],
                        color = gene_label, fill = gene_label)) +
    geom_point(size = 2.5, alpha = 0.70, shape = 16) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 1.0) +
    scale_color_manual(
      values = colors,
      name   = "Operon gene",
      labels = setNames(rho_labels$lab, rho_labels$gene_label)
    ) +
    scale_fill_manual(
      values = colors,
      name   = "Operon gene",
      labels = setNames(rho_labels$lab, rho_labels$gene_label)
    ) +
    labs(
      title = title_label,
      x     = "KO Abundance (carrying-bin weighted)",
      y     = y_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "right",
      legend.text     = element_text(size = 8),
      legend.key.size = unit(0.7, "lines"),
      legend.title    = element_text(size = 9, face = "bold")
    )

  return(p)
}

operon_plot_list <- list()

for (metric in metrics) {
  cat("  Building operon plot for:", metric, "\n")

  p <- make_operon_plot(
    data        = t6ss_operon,
    metric      = metric,
    y_label     = metric_labs[[metric]],
    title_label = paste("T6SS Operon vs.", metric_labs[[metric]]),
    colors      = operon_colors,
    spearman_df = spearman_operon
  )

  operon_plot_list[[metric]] <- p

  pdf(paste0(OPERON_PREFIX, "_", metric, ".pdf"), width = 8, height = 6)
  print(p)
  dev.off()

  png(paste0(OPERON_PREFIX, "_", metric, ".png"),
      width = 8, height = 6, units = "in", res = 300)
  print(p)
  dev.off()

  cat("    Saved:", metric, "(PDF + PNG)\n")
}

# Combined 3-panel
p_operon_combined <- (operon_plot_list$richness +
                      operon_plot_list$shannon  +
                      operon_plot_list$simpson) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(
    title    = "T6SS Complete Operon (bin-identified) vs. Bin Alpha Diversity",
    subtitle = paste0(
      "KOs: ", paste(OPERON_KOS, collapse = ", "),
      " | n = ", length(common_samples), " samples"
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
    )
  ) &
  theme(legend.position = "right")

pdf(paste0(OPERON_PREFIX, "_combined.pdf"), width = 20, height = 6)
print(p_operon_combined)
dev.off()

png(paste0(OPERON_PREFIX, "_combined.png"),
    width = 20, height = 6, units = "in", res = 300)
print(p_operon_combined)
dev.off()

svg(paste0(OPERON_PREFIX, "_combined.svg"), width = 20, height = 6)
print(p_operon_combined)
dev.off()

cat("  Saved: operon combined 3-panel (PDF + PNG + SVG)\n")

cat("\n  Operon Spearman correlations:\n")
cat(sprintf("  %-28s %-10s %7s %10s %5s\n",
            "gene_label", "metric", "rho", "p-value", "n"))
cat("  ", paste(rep("-", 66), collapse=""), "\n")
spearman_operon %>%
  arrange(desc(abs(rho))) %>%
  rowwise() %>%
  group_walk(~ {
    sig <- if (!is.na(.x$pval) && .x$pval < 0.05) " *" else ""
    cat(sprintf("  %-28s %-10s %+7.4f %10.4g %5d%s\n",
                substr(.x$gene_label, 1, 28), .x$metric,
                .x$rho, .x$pval, .x$n, sig))
  })

# ============================================================================
# PHASE 13: T6SS GENES vs. ENDOPHYTE ABUNDANCE PROPORTION
# ============================================================================

cat("\n=== PHASE 13: T6SS Genes vs. Endophyte Abundance Proportion ===\n")

# Compute endophyte proportion per sample:
# sum of Relative_Abundance (0-100 scale) for Endophyte bins / 100
endophyte_prop <- bin_abund %>%
  group_by(Sample) %>%
  summarise(
    endophyte_prop = sum(Relative_Abundance[Organism == "Endophyte"], na.rm = TRUE) / 100,
    .groups = "drop"
  ) %>%
  rename(sample_id = Sample)

cat("  Endophyte proportion computed for", nrow(endophyte_prop), "samples\n")
cat("  Range: [", round(min(endophyte_prop$endophyte_prop), 4), ",",
    round(max(endophyte_prop$endophyte_prop), 4), "]\n")
print(endophyte_prop %>% arrange(sample_id))

# Merge with T6SS long data
t6ss_endo <- t6ss_plot %>%
  select(sample_id, KO, gene_label, abundance, treatment) %>%
  left_join(endophyte_prop, by = "sample_id")

cat("  Merged rows:", nrow(t6ss_endo), "\n")

# Spearman correlations: each KO vs endophyte proportion
spearman_endo <- data.frame(KO = valid_kos$KO, stringsAsFactors = FALSE) %>%
  left_join(T6SS_KOS, by = "KO") %>%
  rowwise() %>%
  mutate(
    vals_abund = list(t6ss_endo$abundance[t6ss_endo$KO == KO]),
    vals_endo  = list(t6ss_endo$endophyte_prop[t6ss_endo$KO == KO]),
    n          = length(vals_abund),
    ct         = list(tryCatch(
      cor.test(vals_abund, vals_endo, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )),
    rho  = if (!is.null(ct)) unname(ct$estimate) else NA_real_,
    pval = if (!is.null(ct)) ct$p.value           else NA_real_
  ) %>%
  select(KO, gene_label, n, rho, pval) %>%
  ungroup()

cat("\n  Spearman: T6SS KO vs Endophyte Proportion (sorted by |rho|):\n")
cat(sprintf("  %-28s %7s %10s %5s\n", "gene_label", "rho", "p-value", "n"))
cat("  ", paste(rep("-", 55), collapse=""), "\n")
spearman_endo %>%
  arrange(desc(abs(rho))) %>%
  rowwise() %>%
  group_walk(~ {
    sig <- if (!is.na(.x$pval) && .x$pval < 0.05) " *" else ""
    cat(sprintf("  %-28s %+7.4f %10.4g %5d%s\n",
                substr(.x$gene_label, 1, 28), .x$rho, .x$pval, .x$n, sig))
  })

ENDO_PREFIX <- file.path(OUTPUT_DIR, "T6SS_endophyte_regression")

# ── All 15 KOs plot ──────────────────────────────────────────────────────────
top3_endo <- spearman_endo %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = 3) %>%
  mutate(lab = sprintf("%s: ρ=%+.2f, p=%.3g", substr(gene_label, 1, 18), rho, pval))
caption_endo <- paste0("Top |ρ| genes:\n", paste(top3_endo$lab, collapse = "\n"))

p_endo_all <- ggplot(t6ss_endo,
                     aes(x = abundance, y = endophyte_prop,
                         color = gene_label, fill = gene_label)) +
  geom_point(size = 2, alpha = 0.65, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.10, linewidth = 0.85) +
  scale_color_manual(values = t6ss_colors, name = "T6SS gene") +
  scale_fill_manual( values = t6ss_colors, name = "T6SS gene") +
  labs(
    title   = "T6SS Genes vs. Endophyte Abundance Proportion",
    x       = "KO Abundance (carrying-bin weighted)",
    y       = "Endophyte Proportion (of total community)",
    caption = caption_endo
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.caption    = element_text(size = 7.5, color = "grey40", hjust = 0),
    legend.position = "right",
    legend.text     = element_text(size = 7.5),
    legend.key.size = unit(0.55, "lines"),
    legend.title    = element_text(size = 8.5, face = "bold")
  )

pdf(paste0(ENDO_PREFIX, "_all_KOs.pdf"), width = 9, height = 6)
print(p_endo_all)
dev.off()
png(paste0(ENDO_PREFIX, "_all_KOs.png"),
    width = 9, height = 6, units = "in", res = 300)
print(p_endo_all)
dev.off()
cat("  Saved: all-KO endophyte plot (PDF + PNG)\n")

# ── Operon subset plot ───────────────────────────────────────────────────────
t6ss_endo_operon <- t6ss_endo %>% filter(KO %in% OPERON_KOS)
spearman_endo_operon <- spearman_endo %>% filter(KO %in% OPERON_KOS)

rho_endo_labels <- spearman_endo_operon %>%
  arrange(desc(abs(rho))) %>%
  mutate(lab = sprintf("%s\nρ=%+.2f, p=%.3g",
                       substr(gene_label, 1, 20), rho, pval))

p_endo_operon <- ggplot(t6ss_endo_operon,
                        aes(x = abundance, y = endophyte_prop,
                            color = gene_label, fill = gene_label)) +
  geom_point(size = 2.5, alpha = 0.70, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.12, linewidth = 1.0) +
  scale_color_manual(
    values = operon_colors,
    name   = "Operon gene",
    labels = setNames(rho_endo_labels$lab, rho_endo_labels$gene_label)
  ) +
  scale_fill_manual(
    values = operon_colors,
    name   = "Operon gene",
    labels = setNames(rho_endo_labels$lab, rho_endo_labels$gene_label)
  ) +
  labs(
    title = "T6SS Operon vs. Endophyte Abundance Proportion",
    x     = "KO Abundance (carrying-bin weighted)",
    y     = "Endophyte Proportion (of total community)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.7, "lines"),
    legend.title    = element_text(size = 9, face = "bold")
  )

pdf(paste0(ENDO_PREFIX, "_operon.pdf"), width = 9, height = 6)
print(p_endo_operon)
dev.off()
png(paste0(ENDO_PREFIX, "_operon.png"),
    width = 9, height = 6, units = "in", res = 300)
print(p_endo_operon)
dev.off()
cat("  Saved: operon endophyte plot (PDF + PNG)\n")

# ── Side-by-side combined ────────────────────────────────────────────────────
p_endo_combined <- (p_endo_all + p_endo_operon) +
  plot_layout(ncol = 2, guides = "keep") +
  plot_annotation(
    title = "T6SS Abundance vs. Endophyte Community Proportion",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
  )

pdf(paste0(ENDO_PREFIX, "_combined.pdf"), width = 18, height = 6)
print(p_endo_combined)
dev.off()
png(paste0(ENDO_PREFIX, "_combined.png"),
    width = 18, height = 6, units = "in", res = 300)
print(p_endo_combined)
dev.off()
svg(paste0(ENDO_PREFIX, "_combined.svg"), width = 18, height = 6)
print(p_endo_combined)
dev.off()
cat("  Saved: combined endophyte panel (PDF + PNG + SVG)\n")

# ============================================================================
# DONE
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("DONE\n")
cat("========================================================================\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
cat("Files produced:\n")
for (f in sort(list.files(OUTPUT_DIR, pattern = "T6SS_bin_diversity_regression"))) {
  cat(" -", f, "\n")
}
cat("========================================================================\n")
