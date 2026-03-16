#!/usr/bin/env Rscript

# ============================================================================
# CONFIGURATION - SET THESE VARIABLES BEFORE RUNNING
# ============================================================================

# Choose your input mode (change TRUE/FALSE):
include_unmapped <- TRUE  # Set to TRUE for unmapped, FALSE for renormalized

# Or directly specify the input file:
# input_file <- "combined_abundance_long_with_unmapped.tsv"
# input_file <- "combined_abundance_long_renormalized.tsv"

# Automatic file selection based on include_unmapped flag
if (include_unmapped) {
  input_file <- "combined_abundance_long_with_unmapped.tsv"
} else {
  input_file <- "combined_abundance_long_renormalized.tsv"
}

# ============================================================================
# MAIN SCRIPT STARTS HERE
# ============================================================================

cat("=====================================\n")
cat("Configuration\n")
cat("=====================================\n")
cat("Input file:", input_file, "\n")
cat("Include unmapped:", include_unmapped, "\n\n")

# Check if file exists
if (!file.exists(input_file)) {
  stop("Error: Input file '", input_file, "' not found!\n",
       "Please make sure the file is in your working directory.\n",
       "Current working directory: ", getwd())
}

# Load required libraries
cat("Loading required libraries...\n")
required_packages <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer", "scales", "viridis")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("All libraries loaded successfully!\n\n")

# Read data
data <- read.delim(input_file, stringsAsFactors = FALSE)

cat("Data loaded successfully\n")
cat("Total rows:", nrow(data), "\n")
cat("Unique samples:", length(unique(data$Sample)), "\n")
cat("Unique bins:", length(unique(data$Bin)), "\n")

# Prepare taxonomy labels
data <- data %>%
  mutate(
    Best_Taxonomy_Label = case_when(
      Bin == "unmapped" ~ "Unmapped",
      is.na(classification) | classification == "" ~ "Unclassified",
      TRUE ~ classification
    )
  )

# Function to extract the most specific taxonomic rank
extract_taxonomy <- function(tax_string) {
  if (is.na(tax_string) || tax_string == "" || tax_string == "Unmapped") {
    return(tax_string)
  }
  
  # Split by semicolon
  parts <- strsplit(tax_string, ";")[[1]]
  
  # Find the last part that has a name
  for (i in length(parts):1) {
    part <- trimws(parts[i])
    if (grepl("__", part)) {
      name <- sub("^[a-z]__", "", part)
      if (name != "") {
        rank <- sub("__.*", "", part)
        rank_label <- switch(rank,
                             "s" = "",  # species
                             "g" = "",  # genus  
                             "f" = "f. ", # family
                             "o" = "o. ", # order
                             "c" = "c. ", # class
                             "p" = "p. ", # phylum
                             "d" = "d. ", # domain
                             "")
        return(paste0(rank_label, name))
      }
    }
  }
  return(tax_string)
}

# Extract simplified taxonomy labels
data <- data %>%
  rowwise() %>%
  mutate(Taxonomy_Only = extract_taxonomy(Best_Taxonomy_Label)) %>%
  ungroup()

cat("Unique taxonomy labels:", length(unique(data$Taxonomy_Only)), "\n")
cat("Unique treatments:", length(unique(data$Treatment)), "\n\n")

# Optionally exclude unmapped
if (!include_unmapped) {
  n_before <- nrow(data)
  data <- data %>% filter(Bin != "unmapped")
  n_after <- nrow(data)
  cat("Excluded unmapped entries:", n_before - n_after, "rows\n\n")
}

# Identify bins that have duplicate taxonomy assignments
bin_taxonomy_map <- data %>%
  select(Bin, Taxonomy_Only) %>%
  distinct()

taxonomy_counts <- bin_taxonomy_map %>%
  group_by(Taxonomy_Only) %>%
  summarise(n_bins = n(), bins = paste(Bin, collapse = ", "), .groups = 'drop')

duplicate_taxonomies <- taxonomy_counts %>%
  filter(n_bins > 1) %>%
  pull(Taxonomy_Only)

cat("\nFound", length(duplicate_taxonomies), "taxonomies assigned to multiple bins:\n")
if (length(duplicate_taxonomies) > 0) {
  for (tax in duplicate_taxonomies) {
    bins_with_tax <- bin_taxonomy_map %>%
      filter(Taxonomy_Only == tax) %>%
      pull(Bin)
    cat(sprintf("  %s: %s\n", tax, paste(bins_with_tax, collapse = ", ")))
  }
  cat("\nThese will be labeled with bin names to distinguish them.\n")
}

# Create unique labels
data <- data %>%
  mutate(
    Taxonomy_Unique = if_else(
      Taxonomy_Only %in% duplicate_taxonomies & Bin != "unmapped",
      paste0(Taxonomy_Only, " (", Bin, ")"),
      Taxonomy_Only
    )
  )

data <- data %>%
  select(-Taxonomy_Only) %>%
  rename(Taxonomy_Only = Taxonomy_Unique)

# Map treatment codes to full site names
treatment_names <- c(
  "ces" = "Caesarea",
  "RH" = "Ramat-Hanadiv",
  "carR" = "Carmel (Genista fasselata)",
  "carK" = "Carmel (C. villosa)",
  "hok" = "Hukok",
  "mtz" = "Metzer"
)

data <- data %>%
  mutate(Treatment_Full = if_else(Treatment %in% names(treatment_names),
                                  treatment_names[Treatment],
                                  Treatment))

# Calculate totals before normalization
sample_totals <- data %>%
  group_by(Sample) %>%
  summarise(Total_Abundance = sum(Relative_Abundance), .groups = 'drop')

cat("\nSample totals before normalization:\n")
print(sample_totals)

# Determine if normalization is needed
needs_normalization <- any(abs(sample_totals$Total_Abundance - 100) > 0.1)

if (needs_normalization) {
  cat("\nNormalizing to 100%...\n")
  
  data_clean <- data %>%
    mutate(
      Sample = trimws(Sample),
      Treatment_Full = trimws(Treatment_Full),
      Taxonomy_Only = trimws(Taxonomy_Only),
      Relative_Abundance = as.numeric(Relative_Abundance)
    ) %>%
    replace_na(list(Relative_Abundance = 0)) %>%
    group_by(Sample, Treatment_Full, Taxonomy_Only) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance), .groups = "drop")
  
  data_normalized <- data_clean %>%
    group_by(Sample) %>%
    group_modify(~{
      df <- .x
      total <- sum(df$Relative_Abundance)
      if (total == 0) {
        df$Relative_Abundance <- 0
        return(df)
      }
      
      vals <- df$Relative_Abundance / total * 100
      vals_round <- round(vals, 4)
      diff <- 100 - sum(vals_round)
      
      if (abs(diff) > 1e-8) {
        i <- which.max(df$Relative_Abundance)
        vals_round[i] <- vals_round[i] + diff
      }
      
      df$Relative_Abundance <- vals_round
      df
    }) %>%
    ungroup()
} else {
  cat("\nData already normalized to 100%, skipping normalization...\n")
  data_normalized <- data %>%
    mutate(
      Sample = trimws(Sample),
      Treatment_Full = trimws(Treatment_Full),
      Taxonomy_Only = trimws(Taxonomy_Only),
      Relative_Abundance = as.numeric(Relative_Abundance)
    ) %>%
    replace_na(list(Relative_Abundance = 0)) %>%
    group_by(Sample, Treatment_Full, Taxonomy_Only) %>%
    summarise(Relative_Abundance = sum(Relative_Abundance), .groups = "drop")
}

# Aggregate
data_normalized <- data_normalized %>%
  group_by(Sample, Treatment_Full, Taxonomy_Only) %>%
  summarise(Relative_Abundance = sum(Relative_Abundance), .groups = 'drop')

# Verify normalization
verification <- data_normalized %>%
  group_by(Sample) %>%
  summarise(Total = sum(Relative_Abundance), .groups = 'drop')

cat("\nVerification - all samples should sum to 100%:\n")
print(verification)

not_hundred <- verification %>% filter(abs(Total - 100) > 0.01)
if (nrow(not_hundred) > 0) {
  cat("\nWARNING: Some samples don't sum to 100%:\n")
  print(not_hundred)
} else {
  cat("\n✓ All samples sum to 100%!\n")
}

# Get unique taxonomy labels sorted by abundance
taxonomy_order <- data_normalized %>%
  group_by(Taxonomy_Only) %>%
  summarise(Mean_Abundance = mean(Relative_Abundance), .groups = 'drop') %>%
  arrange(desc(Mean_Abundance)) %>%
  pull(Taxonomy_Only)

# Force "Unmapped" to be first (will be at top of stacked bars after reversal)
if ("Unmapped" %in% taxonomy_order) {
  taxonomy_order <- c("Unmapped", setdiff(taxonomy_order, "Unmapped"))
}

n_taxa <- length(taxonomy_order)
cat("\nNumber of unique taxa:", n_taxa, "\n")

# Generate colors
generate_scientific_colors <- function(n) {
  if (n <= 8) {
    colors <- brewer.pal(min(n, 8), "Dark2")
  } else if (n <= 12) {
    colors <- brewer.pal(min(n, 12), "Paired")
  } else if (n <= 20) {
    colors <- c(
      brewer.pal(8, "Dark2"),
      brewer.pal(8, "Set2"),
      brewer.pal(min(n - 16, 8), "Accent")
    )
  } else if (n <= 30) {
    colors <- viridis(n, option = "turbo")
  } else {
    base_colors <- c(
      brewer.pal(8, "Dark2"),
      brewer.pal(8, "Set2"),
      brewer.pal(8, "Accent"),
      viridis(n - 24, option = "turbo")
    )
    colors <- base_colors[1:n]
  }
  return(colors[1:n])
}

taxa_colors <- generate_scientific_colors(n_taxa)
names(taxa_colors) <- taxonomy_order

# IMPORTANT: Only "Unmapped" gets grey color - all other taxa get scientific colors
if ("Unmapped" %in% taxonomy_order) {
  taxa_colors["Unmapped"] <- "#CCCCCC"  # Light grey ONLY for unmapped reads
  cat("\n✓ Unmapped will be shown in grey at the top of bars\n")
}

# Order samples
sample_levels <- data_normalized %>%
  select(Sample, Treatment_Full) %>%
  distinct() %>%
  arrange(Treatment_Full, Sample) %>%
  pull(Sample)

data_normalized$Sample <- factor(data_normalized$Sample, levels = sample_levels)

# Create sample labels
data_normalized <- data_normalized %>%
  mutate(Sample_Label = sub("^[^_]+_", "", as.character(Sample)))

data_normalized$Taxonomy_Only <- factor(data_normalized$Taxonomy_Only, levels = rev(taxonomy_order))

# Create unified plot
p <- ggplot(data_normalized, aes(x = Sample, y = Relative_Abundance, fill = Taxonomy_Only)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = taxa_colors, breaks = taxonomy_order) +
  scale_x_discrete(labels = setNames(data_normalized$Sample_Label, data_normalized$Sample)) +
  labs(
    title = "Taxonomic Composition - All Samples",
    x = "Replicate",
    y = "Relative Abundance (%)",
    fill = "Taxonomy"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100.01), expand = c(0, 0))

# Add treatment separators
sample_positions <- data.frame(
  Sample = sample_levels,
  Position = seq_along(sample_levels)
)

treatment_positions <- data_normalized %>%
  select(Sample, Treatment_Full) %>%
  distinct() %>%
  left_join(sample_positions, by = "Sample") %>%
  group_by(Treatment_Full) %>%
  summarise(
    min_pos = min(Position),
    max_pos = max(Position),
    mid_pos = (min(Position) + max(Position)) / 2,
    .groups = 'drop'
  )

for (i in 1:(nrow(treatment_positions) - 1)) {
  p <- p + geom_vline(xintercept = treatment_positions$max_pos[i] + 0.5, 
                      color = "gray30", linetype = "solid", linewidth = 0.5)
}

# Faceted plot
p_faceted <- ggplot(data_normalized, aes(x = Sample, y = Relative_Abundance, fill = Taxonomy_Only)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = taxa_colors, breaks = taxonomy_order) +
  scale_x_discrete(labels = setNames(data_normalized$Sample_Label, data_normalized$Sample)) +
  facet_grid(. ~ Treatment_Full, scales = "free_x", space = "free_x") +
  labs(
    title = "Taxonomic Composition - All Samples by Treatment",
    x = "Replicate",
    y = "Relative Abundance (%)",
    fill = "Taxonomy"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100.01), expand = c(0, 0))

# Simplified plot
top_n_taxa <- 15
top_taxa <- taxonomy_order[1:min(top_n_taxa, n_taxa)]

data_simplified <- data_normalized %>%
  mutate(Taxonomy_Simplified = if_else(Taxonomy_Only %in% top_taxa, as.character(Taxonomy_Only), "Other (< 2%)"))

data_simplified_agg <- data_simplified %>%
  group_by(Sample, Treatment_Full, Taxonomy_Simplified) %>%
  summarise(Relative_Abundance = sum(Relative_Abundance), .groups = 'drop')

sample_label_map <- data_normalized %>%
  select(Sample, Sample_Label) %>%
  distinct()

data_simplified_agg <- data_simplified_agg %>%
  left_join(sample_label_map, by = "Sample")

taxonomy_order_simplified <- c(top_taxa, "Other (< 2%)")
data_simplified_agg$Taxonomy_Simplified <- factor(data_simplified_agg$Taxonomy_Simplified, 
                                                  levels = rev(taxonomy_order_simplified))
data_simplified_agg$Sample <- factor(data_simplified_agg$Sample, levels = sample_levels)

taxa_colors_simplified <- c(taxa_colors[top_taxa], "Other (< 2%)" = "#999999")

p_simplified <- ggplot(data_simplified_agg, 
                       aes(x = Sample, y = Relative_Abundance, fill = Taxonomy_Simplified)) +
  geom_bar(stat = "identity", width = 0.9, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = taxa_colors_simplified, breaks = taxonomy_order_simplified) +
  scale_x_discrete(labels = setNames(data_simplified_agg$Sample_Label, data_simplified_agg$Sample)) +
  facet_grid(. ~ Treatment_Full, scales = "free_x", space = "free_x") +
  labs(
    title = "Taxonomic Composition - All Samples by Treatment (Top Taxa)",
    x = "Replicate",
    y = "Relative Abundance (%)",
    fill = "Taxonomy"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dashed"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  scale_y_continuous(breaks = seq(0, 100, 25), limits = c(0, 100.01), expand = c(0, 0))

# Output suffix
suffix <- ifelse(include_unmapped, "_with_unmapped", "_renormalized")

# Save plots
cat("\nSaving plots...\n")
ggsave(paste0("unified_abundance_plot", suffix, ".pdf"), p, width = 20, height = 6, dpi = 300)
ggsave(paste0("unified_abundance_plot", suffix, ".png"), p, width = 20, height = 6, dpi = 300, bg = "white")
ggsave(paste0("unified_abundance_plot", suffix, ".svg"), p, width = 20, height = 6)

ggsave(paste0("unified_abundance_plot_faceted", suffix, ".pdf"), p_faceted, width = 20, height = 6, dpi = 300)
ggsave(paste0("unified_abundance_plot_faceted", suffix, ".png"), p_faceted, width = 20, height = 6, dpi = 300, bg = "white")
ggsave(paste0("unified_abundance_plot_faceted", suffix, ".svg"), p_faceted, width = 20, height = 6)

ggsave(paste0("unified_abundance_plot_simplified", suffix, ".pdf"), p_simplified, width = 20, height = 6, dpi = 300)
ggsave(paste0("unified_abundance_plot_simplified", suffix, ".png"), p_simplified, width = 20, height = 6, dpi = 300, bg = "white")
ggsave(paste0("unified_abundance_plot_simplified", suffix, ".svg"), p_simplified, width = 20, height = 6)

# Summary statistics
cat("\n=====================================\n")
cat("Summary Statistics\n")
cat("=====================================\n")
cat("Total samples:", length(unique(data_normalized$Sample)), "\n")
cat("Total treatments:", length(unique(data_normalized$Treatment_Full)), "\n")
cat("Total unique taxa:", n_taxa, "\n")

if (include_unmapped && "Unmapped" %in% data_normalized$Taxonomy_Only) {
  unmapped_stats <- data_normalized %>%
    filter(Taxonomy_Only == "Unmapped") %>%
    summarise(
      Mean = mean(Relative_Abundance),
      Median = median(Relative_Abundance),
      Min = min(Relative_Abundance),
      Max = max(Relative_Abundance)
    )
  cat("\nUnmapped reads statistics (%):\n")
  cat("  Mean:", round(unmapped_stats$Mean, 2), "\n")
  cat("  Median:", round(unmapped_stats$Median, 2), "\n")
  cat("  Range:", round(unmapped_stats$Min, 2), "-", round(unmapped_stats$Max, 2), "\n")
}

cat("\nSamples per treatment:\n")
samples_per_treatment <- data_normalized %>%
  select(Sample, Treatment_Full) %>%
  distinct() %>%
  group_by(Treatment_Full) %>%
  summarise(Count = n(), .groups = 'drop')
print(samples_per_treatment)

cat("\n=====================================\n")
cat("All plots generated successfully!\n")
cat("=====================================\n")
cat("\nGenerated files (in working directory):\n")
cat("Working directory:", getwd(), "\n\n")
cat("1. unified_abundance_plot", suffix, ".* - All samples in one plot\n", sep = "")
cat("2. unified_abundance_plot_faceted", suffix, ".* - Samples grouped by treatment\n", sep = "")
cat("3. unified_abundance_plot_simplified", suffix, ".* - Top taxa only\n", sep = "")
cat("\nEach in formats: .pdf, .png, .svg\n")