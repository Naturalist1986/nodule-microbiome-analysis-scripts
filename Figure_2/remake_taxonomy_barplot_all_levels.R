#!/usr/bin/env Rscript

setwd("G:/My Drive/Moshe/Efrat_Guy_Project/Singlem")

# ============================================================================
# Script to Remake taxonomy_barplot_all_levels.pdf
# ============================================================================
# This script reads SingleM profile TSV files and creates a stacked barplot
# showing taxonomic composition at all taxonomic levels (species, genus,
# family, order, class, phylum, domain) for each sample.
# ============================================================================

library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(png)
library(grid)
library(rsvg)

# ===================================================================
# CONFIGURATION
# ===================================================================

# Directory containing profile TSV files
profiles_dir <- "profiles"

# Treatment definitions (abbreviated, matching Procrustes figure)
treatment_labels <- c(
  "carK" = "CAR-C",
  "carR" = "CAR-G",
  "ces"  = "CES",
  "hok"  = "HUK",
  "mtz"  = "MTZ",
  "RH"   = "RH"
)

treatment_colors <- c(
  "CAR-C" = "#3C5488",
  "CAR-G" = "#00A087",
  "CES"   = "#E64B35",
  "HUK"   = "#F39B7F",
  "MTZ"   = "#8491B4",
  "RH"    = "#4DBBD5"
)

# ===================================================================
# HELPER FUNCTIONS
# ===================================================================

# Extract treatment code from sample name
extract_treatment <- function(sample_name) {
  treatment_code <- str_extract(sample_name, "^[^_]+")
  return(treatment_code)
}

# Extract the most specific taxonomic level from a taxonomy string
extract_most_specific_taxon <- function(taxonomy_string) {
  # Split by semicolon and get individual levels
  levels <- str_split(taxonomy_string, "; ")[[1]]

  # Remove "Root" if present
  levels <- levels[levels != "Root"]

  # Get the last (most specific) level
  if (length(levels) > 0) {
    last_level <- levels[length(levels)]
    # Remove rank prefix (e.g., "s__", "g__", "f__", etc.)
    taxon <- str_replace(last_level, "^[a-z]__", "")
    return(taxon)
  } else {
    return("Unknown")
  }
}

# ===================================================================
# LOAD AND PROCESS DATA
# ===================================================================

message("\n", paste(rep("=", 70), collapse=""))
message("LOADING SINGLEM PROFILE DATA")
message(paste(rep("=", 70), collapse=""))

# Get list of all profile TSV files
profile_files <- list.files(profiles_dir, pattern = "\\.tsv$", full.names = TRUE)
profile_files <- profile_files[!str_detect(profile_files, "taxonomy_barplot")]

message("Found ", length(profile_files), " profile files")

# Read and combine all profiles
all_data <- map_dfr(profile_files, function(file) {
  message("Reading: ", basename(file))

  # Read the file
  data <- read_tsv(file, show_col_types = FALSE)

  # Ensure lowercase column names for consistency
  colnames(data) <- tolower(colnames(data))

  # Extract sample name from filename (remove extension and path)
  sample_name <- str_replace(basename(file), "\\.tsv$", "")
  # Extract the base sample name (e.g., "carK_1" from "carK_1_DKDN250001546-1A_22KVLLLT4_L6")
  sample_base <- str_extract(sample_name, "^[^_]+_[0-9]+")

  data %>%
    mutate(
      sample_name = sample_base,
      file_source = basename(file)
    )
})

message("Total rows loaded: ", nrow(all_data))
message("Column names: ", paste(colnames(all_data), collapse=", "))

# ===================================================================
# EXTRACT MOST SPECIFIC TAXON FOR EACH ENTRY
# ===================================================================

message("\n", paste(rep("=", 70), collapse=""))
message("EXTRACTING TAXONOMIC INFORMATION")
message(paste(rep("=", 70), collapse=""))

# Verify required columns exist
if (!"taxonomy" %in% colnames(all_data)) {
  stop("ERROR: 'taxonomy' column not found in data. Available columns: ",
       paste(colnames(all_data), collapse=", "))
}
if (!"coverage" %in% colnames(all_data)) {
  stop("ERROR: 'coverage' column not found in data. Available columns: ",
       paste(colnames(all_data), collapse=", "))
}

# For each taxonomy string, extract the most specific taxon
# Use .data$ to avoid scoping issues when running interactively
taxonomy_data <- all_data %>%
  mutate(
    most_specific_taxon = map_chr(.data$taxonomy, extract_most_specific_taxon),
    treatment_code = map_chr(.data$sample_name, extract_treatment),
    treatment_label = treatment_labels[treatment_code]
  ) %>%
  # Remove entries that are just domain-level or too generic
  filter(!is.na(most_specific_taxon))

message("Total unique taxa: ", n_distinct(taxonomy_data$most_specific_taxon))
message("Total samples: ", n_distinct(taxonomy_data$sample_name))

# Aggregate coverage by sample and taxon (sum across all taxonomy strings that end in the same taxon)
taxon_summary <- taxonomy_data %>%
  group_by(sample_name, most_specific_taxon, treatment_code, treatment_label) %>%
  summarise(total_coverage = sum(coverage), .groups = "drop")

# ===================================================================
# SELECT TOP TAXA AND CALCULATE RELATIVE ABUNDANCES
# ===================================================================

message("\n", paste(rep("=", 70), collapse=""))
message("CALCULATING RELATIVE ABUNDANCES")
message(paste(rep("=", 70), collapse=""))

# Calculate total coverage per sample
sample_totals <- taxon_summary %>%
  group_by(sample_name) %>%
  summarise(sample_total = sum(total_coverage), .groups = "drop")

# Join back and calculate relative abundance
taxon_rel_abundance <- taxon_summary %>%
  left_join(sample_totals, by = "sample_name") %>%
  mutate(rel_abundance = (total_coverage / sample_total) * 100)

# Find top taxa across all samples (by mean relative abundance)
top_taxa <- taxon_rel_abundance %>%
  group_by(most_specific_taxon) %>%
  summarise(mean_abundance = mean(rel_abundance), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  head(30) %>%
  pull(most_specific_taxon)

message("Selected top ", length(top_taxa), " taxa by mean abundance")

# Create plot data with "Other" category
plot_data <- taxon_rel_abundance %>%
  mutate(
    display_taxon = if_else(most_specific_taxon %in% top_taxa,
                            most_specific_taxon,
                            "Other (low abundance)")
  ) %>%
  group_by(sample_name, display_taxon, treatment_code, treatment_label) %>%
  summarise(abundance = sum(rel_abundance), .groups = "drop")

# ===================================================================
# ASSIGN COLORS
# ===================================================================

message("\n", paste(rep("=", 70), collapse=""))
message("ASSIGNING COLORS")
message(paste(rep("=", 70), collapse=""))

# Get all unique taxa
all_taxa <- unique(plot_data$display_taxon)
all_taxa <- all_taxa[all_taxa != "Other (low abundance)"]

# Sort taxa so Bradyrhizobium species come first
bradyrhizobium_taxa <- all_taxa[str_detect(all_taxa, "Bradyrhizobium")]
other_taxa <- all_taxa[!str_detect(all_taxa, "Bradyrhizobium")]

sorted_taxa <- c(sort(bradyrhizobium_taxa), sort(other_taxa), "Other (low abundance)")

# Generate color palette - using original figure colors
# Define exact colors from the original figure
original_colors <- c(
  "#984EA3", # Bright pink/red
  "#FF7F00", # Orange
  "#FF6600", # Dark orange
  "#4DAF4A", # Teal/green
  "#FFFF33", # Yellow
  "#E41A1C", # Purple
  "#F781BF", # Light pink/magenta
  "#FDBF6F", # Light orange/peach
  "#A65628", # Brown/tan
  "#6BAED6", # Light blue
  "#377EB8", # Blue
  "#FE9929", # Bright orange
  "#FB9A99", # Salmon pink
  "#CAB2D6", # Light purple
  "#B2DF8A", # Light green
  "#33A02C", # Green
  "#1F78B4"  # Dark blue
)

n_brady <- length(bradyrhizobium_taxa)
n_other <- length(other_taxa)

# Use the original color scheme
if (n_brady > 0) {
  # Use original colors for Bradyrhizobium species
  brady_palette <- colorRampPalette(original_colors)
  brady_colors <- brady_palette(n_brady)
  names(brady_colors) <- sort(bradyrhizobium_taxa)
} else {
  brady_colors <- c()
}

if (n_other > 0) {
  # Use extended original palette for other taxa
  other_palette <- colorRampPalette(c(
    "#B3DE69", # Light yellow-green
    "#FDBF6F", # Light orange
    "#CAB2D6", # Light purple
    "#FB8072", # Light red/salmon
    "#80B1D3", # Light blue
    "#BEBADA", # Lavender
    "#8DD3C7", # Teal
    "#FFFFB3", # Pale yellow
    "#FDB462", # Orange
    "#BC80BD", # Orchid
    "#CCEBC5", # Mint
    "#FFED6F", # Yellow
    "#E78AC3"  # Pink
  ))
  other_colors <- other_palette(n_other)
  names(other_colors) <- sort(other_taxa)
} else {
  other_colors <- c()
}

# Combine colors
all_colors <- c(brady_colors, other_colors, "Other (low abundance)" = "#CCCCCC")

message("Assigned colors to ", length(all_colors), " taxa")

# ===================================================================
# CREATE PLOT
# ===================================================================

message("\n", paste(rep("=", 70), collapse=""))
message("CREATING BARPLOT")
message(paste(rep("=", 70), collapse=""))

# Order samples by treatment and sample number; extract replicate number for x-axis
plot_data <- plot_data %>%
  mutate(
    sample_num    = str_extract(as.character(sample_name), "[0-9]+$"),
    sample_name   = factor(sample_name),
    display_taxon = factor(display_taxon, levels = sorted_taxa),
    treatment_label = factor(treatment_label, levels = treatment_labels)
  ) %>%
  arrange(treatment_label, sample_name)

# Get sample order
sample_order <- unique(plot_data$sample_name)

# Reorder samples and propagate replicate-number ordering to sample_num factor
plot_data <- plot_data %>%
  mutate(
    sample_name = factor(sample_name, levels = sample_order),
    sample_num  = factor(sample_num,  levels = unique(sample_num[order(sample_name)]))
  )

# Create the barplot (x = replicate number only; treatment shown in facet strip)
p <- ggplot(plot_data, aes(x = sample_num, y = abundance, fill = display_taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.95) +
  facet_grid(~ treatment_label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = all_colors, name = "Taxon") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    panel.spacing = unit(0.2, "lines"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    title = "Taxonomic Composition: All Levels",
    x = "Replicate",
    y = "Relative Abundance (%)"
  )

# ===================================================================
# SAVE OUTPUT
# ===================================================================

message("\n", paste(rep("=", 70), collapse=""))
message("SAVING OUTPUT FILES")
message(paste(rep("=", 70), collapse=""))

# Save plot
ggsave("taxonomy_barplot_all_levels_new.pdf", p, width = 16, height = 9)
ggsave("taxonomy_barplot_all_levels_new.png", p, width = 16, height = 9, dpi = 300)

message("✓ Saved: taxonomy_barplot_all_levels_new.pdf")
message("✓ Saved: taxonomy_barplot_all_levels_new.png")

# Save Panel B as true vector to Final_Figures
FINAL_FIG2 <- "G:/My Drive/Moshe/Efrat_Guy_Project/Final_Figures/Figure_2"
dir.create(FINAL_FIG2, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(FINAL_FIG2, "Figure_2B_Taxonomy_Barplot.svg"), p, width = 16, height = 9)
ggsave(file.path(FINAL_FIG2, "Figure_2B_Taxonomy_Barplot.pdf"), p, width = 16, height = 9)
message("✓ Saved: Figure_2B_Taxonomy_Barplot.svg/pdf")

# Save summary table
summary_table <- plot_data %>%
  select(sample_name, treatment_label, display_taxon, abundance) %>%
  arrange(sample_name, desc(abundance))

write_csv(summary_table, "taxonomy_summary_all_levels_new.csv")
message("✓ Saved: taxonomy_summary_all_levels_new.csv")

# ===================================================================
# ASSEMBLE COMBINED FIGURE
# ===================================================================

message("\n", paste(rep("=", 70), collapse=""))
message("ASSEMBLING COMBINED FIGURE")
message(paste(rep("=", 70), collapse=""))

tag_theme <- theme(
  plot.tag          = element_text(face = "bold", size = 16),
  plot.tag.position = "topleft"
)

# Helper: load a PNG file as a ggplot panel
png_to_panel <- function(path, tag_label) {
  img  <- png::readPNG(path)
  grob <- grid::rasterGrob(img, interpolate = TRUE)
  ggplot() +
    annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    labs(tag = tag_label) + theme_void() + tag_theme
}

# Helper: load an SVG file as a ggplot panel
svg_to_panel <- function(path, tag_label) {
  raw_png <- rsvg::rsvg_png(path, width = 1500)
  img     <- png::readPNG(raw_png)
  grob    <- grid::rasterGrob(img, interpolate = TRUE)
  ggplot() +
    annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    labs(tag = tag_label) + theme_void() + tag_theme
}

# Panel A — KO Composition PCoA (produced by make_figure.R in Procrustes/)
panel_a <- png_to_panel("../Procrustes/ko_pcoa_panel_e.png", "A")

# Panel B — taxonomy barplot (produced above)
panel_b <- p + labs(tag = "B") + tag_theme

# Panel C — Species-level PCoA
panel_c <- svg_to_panel("taxonomy_barplot_species_top20_PCoA.svg", "C")

# Panel D — Genus-level PCoA
panel_d <- svg_to_panel("taxonomy_barplot_genus_top20_PCoA.svg", "D")

# Copy species and genus PCoA SVGs to Final_Figures/Figure_2/ as individual vector panels
if (file.exists("taxonomy_barplot_species_top20_PCoA.svg"))
  file.copy("taxonomy_barplot_species_top20_PCoA.svg",
            file.path(FINAL_FIG2, "Figure_2C_Species_PCoA.svg"), overwrite = TRUE)
if (file.exists("taxonomy_barplot_genus_top20_PCoA.svg"))
  file.copy("taxonomy_barplot_genus_top20_PCoA.svg",
            file.path(FINAL_FIG2, "Figure_2D_Genus_PCoA.svg"), overwrite = TRUE)
message("✓ Copied: Figure_2C_Species_PCoA.svg and Figure_2D_Genus_PCoA.svg")

# Layout: A (KO PCoA) | B (barplot) on top; C (species PCoA) | D (genus PCoA) below
combined_fig <- (panel_a | panel_b) /
  (panel_c | panel_d)

ggsave("taxonomy_combined_figure.pdf", combined_fig, width = 18, height = 14)
ggsave("taxonomy_combined_figure.png", combined_fig, width = 18, height = 14, dpi = 300)

message("✓ Saved: taxonomy_combined_figure.pdf")
message("✓ Saved: taxonomy_combined_figure.png")

message("\n", paste(rep("=", 70), collapse=""))
message("ANALYSIS COMPLETE!")
message(paste(rep("=", 70), collapse=""))
message("\nSummary:")
message("  - Samples analyzed: ", n_distinct(plot_data$sample_name))
message("  - Unique taxa displayed: ", length(sorted_taxa))
message("  - Bradyrhizobium taxa: ", n_brady)
message("  - Other taxa: ", n_other)
message("\nOutput files:")
message("  - taxonomy_barplot_all_levels_new.pdf  (standalone barplot)")
message("  - taxonomy_barplot_all_levels_new.png  (standalone barplot)")
message("  - taxonomy_combined_figure.pdf         (A: barplot | B: KO PCoA | C: species PCoA | D: genus PCoA)")
message("  - taxonomy_combined_figure.png")
message("  - taxonomy_summary_all_levels_new.csv")
