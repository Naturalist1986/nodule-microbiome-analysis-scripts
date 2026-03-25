################################################################################
# COMPREHENSIVE NITROGEN METABOLISM GENE-TRAIT ASSOCIATION ANALYSIS
# With Built-in Alternative Multiple Testing Corrections
################################################################################
#
# This script performs end-to-end analysis of nitrogen metabolism genes:
# 1. Load and merge all data
# 2. Perform KO-level regressions
# 3. Apply multiple correction methods (standard FDR, q-value, trait-stratified)
# 4. Analyze symbiont vs endophyte distribution
# 5. Create publication-ready visualizations
# 6. Generate comprehensive reports
#
# Required input files:
# - normalized_kegg_results.xlsx (KO abundances)
# - Guy_measurements.xlsx (plant traits)
# - ko_presence_absence_table.xlsx (KO presence in bins)
# - bin_classification.csv (bin taxonomy and type)
# - List_of_KOs.csv (gene annotations)
#
################################################################################

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(readr)
library(readxl)
library(ggplot2)
library(broom)
library(viridis)
library(patchwork)  # For combining plots
library(ComplexHeatmap)  # For Tier1-style heatmap
library(circlize)
library(grid)
library(svglite)  # For SVG output

# For alternative corrections
library(qvalue)
library(poolr)  # For Meff

setwd("C:/Temp/KO_Analysis")

cat("================================================================================\n")
cat("NITROGEN METABOLISM GENE-TRAIT ASSOCIATION ANALYSIS\n")
cat("With Alternative Multiple Testing Corrections\n")
cat("================================================================================\n\n")

################################################################################
# TREATMENT NAME MAPPING
################################################################################

# Create treatment name mapping
treatment_names <- c(
  "ces"  = "CES",
  "RH"   = "RH",
  "carR" = "CAR-G",
  "carK" = "CAR-C",
  "hok"  = "HUK",
  "mtz"  = "MTZ"
)

# Define treatment order for plots
treatment_order <- c("CAR-C", "CAR-G", "CES", "HUK", "MTZ", "RH")

################################################################################
# 1. LOAD DATA
################################################################################

cat("STEP 1: Loading data...\n")

# Load gene list
gene_list <- read.csv("List_of_KOs.csv", skip = 2)
colnames(gene_list) <- c("Gene", "KEGG_KO", "Category", "Function", "EC_Number", "KEGG_Pathway")

# Load KO abundances (try both possible names)
if (file.exists("normalized_kegg_results.xlsx")) {
  ko_abundance <- read_excel("normalized_kegg_results.xlsx")
  cat("✓ Loaded: normalized_kegg_results.xlsx\n")
} else if (file.exists("normalized_kegg_result.xlsx")) {
  ko_abundance <- read_excel("normalized_kegg_result.xlsx")
  cat("✓ Loaded: normalized_kegg_result.xlsx\n")
} else {
  stop("Error: Could not find normalized_kegg_results.xlsx or normalized_kegg_result.xlsx")
}

# Load plant measurements
measurements <- read_excel("Guy_measurements.xlsx")
cat("✓ Loaded: Guy_measurements.xlsx\n")

# Load KO presence in bins
ko_bins_raw <- read_excel("ko_presence_absence_table.xlsx")

# Check data structure - KOs might be in rows instead of columns
if ("KO" %in% colnames(ko_bins_raw)) {
  # Data has KOs as rows, bins as columns - need to transpose
  cat("  Detected KO presence table with KOs as rows (transposing...)\n")
  
  ko_col <- ko_bins_raw$KO
  ko_bins_matrix <- as.matrix(ko_bins_raw[, -1])  # Remove KO column
  rownames(ko_bins_matrix) <- ko_col
  
  # Transpose so bins are rows, KOs are columns
  ko_bins_t <- t(ko_bins_matrix)
  ko_bins <- as.data.frame(ko_bins_t)
  ko_bins$bin <- rownames(ko_bins_t)
  rownames(ko_bins) <- NULL
  
  # Reorder columns to have bin first
  ko_bins <- ko_bins %>% select(bin, everything())
  
  cat(sprintf("  Transposed: %d bins × %d KOs\n", nrow(ko_bins), ncol(ko_bins) - 1))
} else {
  # Data already has bins as rows
  ko_bins <- ko_bins_raw
  cat(sprintf("  Loaded: %d bins × %d KOs\n", nrow(ko_bins), ncol(ko_bins) - 1))
}
ko_bins <- ko_bins %>%
  mutate(across(-bin, ~ifelse(. > 0, 1, 0)))
cat("✓ Loaded: ko_presence_absence_table.xlsx\n")

# Load bin classification from combined abundance file
cat("  Loading bin data from combined_abundance_long_renormalized.tsv...\n")
abundance_data <- read_tsv("combined_abundance_long_renormalized.tsv", show_col_types = FALSE)

# Create bin_class from abundance data
bin_class <- abundance_data %>%
  distinct(Bin_Prefixed, .keep_all = TRUE) %>%
  transmute(
    bin = Bin_Prefixed,
    treatment = Treatment,
    organism_type = tolower(Organism),  # Convert to lowercase (Symbiont -> symbiont)
    classification = classification,
    genus = str_extract(classification, "g__([^;]+)") %>% str_remove("g__")
  )
cat("✓ Loaded: combined_abundance_long_renormalized.tsv\n")
cat(sprintf("  Created bin classification for %d bins\n", nrow(bin_class)))

cat(sprintf("\nData loaded successfully:\n"))
cat(sprintf("  - %d nitrogen metabolism genes to analyze\n", nrow(gene_list)))
cat(sprintf("  - %d samples with KO abundances\n", ncol(ko_abundance) - 1))
cat(sprintf("  - %d plant measurements\n", nrow(measurements)))
cat(sprintf("  - %d metagenome-assembled genomes (bins)\n", nrow(bin_class)))

################################################################################
# 2. PROCESS AND MERGE DATA
################################################################################

cat("\nSTEP 2: Processing and merging data...\n")

# Extract sample info from column names
sample_names <- colnames(ko_abundance)[-1]
extract_sample_info <- function(sample_name) {
  parts <- str_split(sample_name, "_")[[1]]
  treatment <- parts[1]
  replicate <- as.numeric(parts[2])
  return(c(treatment, replicate))
}

sample_info <- t(sapply(sample_names, extract_sample_info))
colnames(sample_info) <- c("Treatment", "Replicate")
sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
sample_info$Replicate <- as.numeric(sample_info$Replicate)
sample_info$Sample <- sample_names

# Create clean KO abundance matrix
ko_abund_matrix <- ko_abundance %>%
  select(-kegg_number) %>%
  t() %>%
  as.data.frame()
colnames(ko_abund_matrix) <- ko_abundance$kegg_number
ko_abund_matrix$Sample <- rownames(ko_abund_matrix)
rownames(ko_abund_matrix) <- NULL

# Merge all data
ko_abund_clean <- merge(sample_info, ko_abund_matrix, by = "Sample")
data_merged <- merge(ko_abund_clean, measurements, 
                     by = c("Treatment", "Replicate"), 
                     all = FALSE)

# Subset to genes of interest
genes_of_interest <- gene_list$KEGG_KO
ko_cols <- intersect(genes_of_interest, colnames(data_merged))

cat(sprintf("✓ Merged %d samples with complete data\n", nrow(data_merged)))
cat(sprintf("✓ Found %d/%d genes of interest in the data\n", 
            length(ko_cols), length(genes_of_interest)))

################################################################################
# 3. KO-LEVEL REGRESSION ANALYSIS
################################################################################

cat("\nSTEP 3: Running KO-level regressions...\n")

trait_cols <- c("percent_N", "Plant_Biomass", "Fixation_per_Nodule", "NMF")

# Run regressions for all KO-trait combinations
regression_results <- list()

for (ko in ko_cols) {
  for (trait in trait_cols) {
    tryCatch({
      # Linear model
      formula_str <- paste(trait, "~", ko)
      model <- lm(as.formula(formula_str), data = data_merged)
      
      # Extract results
      summary_model <- summary(model)
      coef_data <- coef(summary_model)[2, ]  # Coefficient for the KO
      
      regression_results[[length(regression_results) + 1]] <- data.frame(
        KO = ko,
        Trait = trait,
        estimate = coef_data[1],
        std.error = coef_data[2],
        statistic = coef_data[3],
        p.value = coef_data[4],
        R_squared = summary_model$r.squared,
        Adj_R_squared = summary_model$adj.r.squared
      )
    }, error = function(e) {
      # Skip if model fails
    })
  }
}

# Combine results
ko_results <- do.call(rbind, regression_results)
rownames(ko_results) <- NULL

# Add gene annotations
ko_results <- ko_results %>%
  left_join(gene_list %>% 
              distinct(KEGG_KO, .keep_all = TRUE), 
            by = c("KO" = "KEGG_KO"),
            relationship = "many-to-one") %>%
  select(KO, Gene, Category, Function, Trait, estimate, std.error, 
         statistic, p.value, R_squared, Adj_R_squared)

cat(sprintf("✓ Completed %d regression tests (%d KOs × %d traits)\n", 
            nrow(ko_results), length(ko_cols), length(trait_cols)))

################################################################################
# 4. MULTIPLE TESTING CORRECTIONS
################################################################################

cat("\nSTEP 4: Applying multiple testing corrections...\n")

# 4.1 Standard Benjamini-Hochberg FDR
ko_results$p_adj_standard <- p.adjust(ko_results$p.value, method = "BH")
n_sig_standard <- sum(ko_results$p_adj_standard < 0.05)
cat(sprintf("  Standard BH-FDR (< 0.05): %d significant\n", n_sig_standard))

# 4.2 Storey's Q-value
cat("  Calculating Storey's q-value...\n")
tryCatch({
  qobj <- qvalue(ko_results$p.value, lambda = seq(0.05, 0.95, 0.05))
  ko_results$q_value <- qobj$qvalues
  ko_results$local_fdr <- qobj$lfdr
  
  n_sig_qvalue <- sum(ko_results$q_value < 0.05)
  cat(sprintf("    π₀ (proportion of true nulls): %.3f\n", qobj$pi0))
  cat(sprintf("    True alternatives: %.1f%%\n", 100 * (1 - qobj$pi0)))
  cat(sprintf("    Q-value (< 0.05): %d significant (+%d vs standard)\n", 
              n_sig_qvalue, n_sig_qvalue - n_sig_standard))
}, error = function(e) {
  ko_results$q_value <- NA
  ko_results$local_fdr <- NA
  cat("    ⚠ Could not calculate q-values (p-value distribution not suitable)\n")
})

# 4.3 Trait-Stratified FDR
cat("  Calculating trait-stratified FDR...\n")
ko_results <- ko_results %>%
  group_by(Trait) %>%
  mutate(p_adj_trait_stratified = p.adjust(p.value, method = "BH")) %>%
  ungroup()

n_sig_stratified <- sum(ko_results$p_adj_trait_stratified < 0.05)
cat(sprintf("    Trait-stratified FDR (< 0.05): %d significant (+%d vs standard)\n", 
            n_sig_stratified, n_sig_stratified - n_sig_standard))

# 4.4 Effective Number of Tests (Meff) - if possible
cat("  Calculating effective number of tests (Meff)...\n")
tryCatch({
  # Get correlation matrix
  ko_cols_in_data <- intersect(ko_cols, colnames(data_merged))
  if (length(ko_cols_in_data) > 10) {
    ko_matrix <- data_merged %>%
      select(all_of(ko_cols_in_data)) %>%
      as.matrix() %>%
      t()
    
    cor_matrix <- cor(t(ko_matrix), use = "pairwise.complete.obs", method = "spearman")
    meff_galwey <- meff(cor_matrix, method = "galwey")
    
    n_traits <- length(unique(ko_results$Trait))
    effective_total_tests <- meff_galwey * n_traits
    
    ko_results$p_adj_meff <- p.adjust(ko_results$p.value, method = "BH", 
                                      n = round(effective_total_tests))
    
    n_sig_meff <- sum(ko_results$p_adj_meff < 0.05)
    cat(sprintf("    Effective tests: %.1f (vs %d actual)\n", 
                effective_total_tests, nrow(ko_results)))
    cat(sprintf("    Meff-adjusted FDR (< 0.05): %d significant (+%d vs standard)\n",
                n_sig_meff, n_sig_meff - n_sig_standard))
  } else {
    ko_results$p_adj_meff <- NA
    cat("    ⚠ Not enough KOs for Meff calculation\n")
  }
}, error = function(e) {
  ko_results$p_adj_meff <- NA
  cat("    ⚠ Could not calculate Meff\n")
})

# 4.5 Identify which methods call each significant
# Check which columns exist
has_meff <- "p_adj_meff" %in% colnames(ko_results)
has_qvalue <- "q_value" %in% colnames(ko_results)

ko_results <- ko_results %>%
  mutate(
    Sig_standard = p_adj_standard < 0.05,
    Sig_meff = if(has_meff) !is.na(p_adj_meff) & p_adj_meff < 0.05 else FALSE,
    Sig_qvalue = if(has_qvalue) !is.na(q_value) & q_value < 0.05 else FALSE,
    Sig_trait_stratified = p_adj_trait_stratified < 0.05,
    N_methods_significant = as.integer(Sig_standard) + as.integer(Sig_meff) + 
      as.integer(Sig_qvalue) + as.integer(Sig_trait_stratified)
  )

# Save complete results
write.csv(ko_results, "KO_results_all_methods.csv", row.names = FALSE)
cat("\n✓ Saved: KO_results_all_methods.csv\n")

################################################################################
# 5. IDENTIFY ROBUST ASSOCIATIONS
################################################################################

cat("\nSTEP 5: Identifying robust associations...\n")

robust_associations <- ko_results %>%
  filter(N_methods_significant >= 2) %>%
  arrange(desc(N_methods_significant), p.value)

cat(sprintf("✓ Found %d robust associations (≥2 methods significant):\n", 
            nrow(robust_associations)))
cat(sprintf("    - 4 methods: %d\n", sum(robust_associations$N_methods_significant == 4)))
cat(sprintf("    - 3 methods: %d\n", sum(robust_associations$N_methods_significant == 3)))
cat(sprintf("    - 2 methods: %d\n", sum(robust_associations$N_methods_significant == 2)))

write.csv(robust_associations, "Robust_associations.csv", row.names = FALSE)
cat("✓ Saved: Robust_associations.csv\n")

################################################################################
# 6. SYMBIONT VS ENDOPHYTE DISTRIBUTION
################################################################################

cat("\nSTEP 6: Analyzing symbiont vs endophyte distribution...\n")

# Get bin names by type
symbiont_bins <- bin_class %>% filter(organism_type == "symbiont") %>% pull(bin)
endophyte_bins <- bin_class %>% filter(organism_type == "endophyte") %>% pull(bin)

cat(sprintf("  Symbionts (Bradyrhizobium): %d bins\n", length(symbiont_bins)))
cat(sprintf("  Endophytes (other bacteria): %d bins\n", length(endophyte_bins)))

# Function to analyze KO distribution
analyze_ko_distribution <- function(ko_list, ko_bins_df, bin_class_df) {
  
  results <- list()
  
  for (ko in ko_list) {
    # Check if KO exists in ko_bins
    if (ko %in% colnames(ko_bins_df)) {
      # Get bins containing this KO
      ko_data <- ko_bins_df %>% select(bin, all_of(ko))
      bins_with_ko <- ko_data %>% filter(.data[[ko]] == 1) %>% pull(bin)
      
      if (length(bins_with_ko) > 0) {
        # Match with classification
        bin_info <- bin_class_df %>% filter(bin %in% bins_with_ko)
        
        # Count by type
        n_symbiont <- sum(bin_info$organism_type == "symbiont", na.rm = TRUE)
        n_endophyte <- sum(bin_info$organism_type == "endophyte", na.rm = TRUE)
        
        # Get treatments
        treatments_symbiont <- bin_info %>% 
          filter(organism_type == "symbiont") %>% 
          pull(treatment) %>% 
          unique() %>% 
          paste(collapse = ", ")
        
        treatments_endophyte <- bin_info %>% 
          filter(organism_type == "endophyte") %>% 
          pull(treatment) %>% 
          unique() %>% 
          paste(collapse = ", ")
        
        # Classify primary location
        total_bins <- n_symbiont + n_endophyte
        primarily_in <- case_when(
          total_bins == 0 ~ "Not found",
          n_symbiont > n_endophyte * 2 ~ "Mainly Symbiont",
          n_endophyte > n_symbiont * 2 ~ "Mainly Endophyte",
          n_symbiont > 0 & n_endophyte > 0 ~ "Both",
          n_symbiont > 0 ~ "Symbiont only",
          n_endophyte > 0 ~ "Endophyte only",
          TRUE ~ "Not found"
        )
        
        results[[ko]] <- data.frame(
          KO = ko,
          N_bins_symbiont = n_symbiont,
          N_bins_endophyte = n_endophyte,
          Total_bins = total_bins,
          Treatments_symbiont = treatments_symbiont,
          Treatments_endophyte = treatments_endophyte,
          Primarily_in = primarily_in
        )
      } else {
        results[[ko]] <- data.frame(
          KO = ko,
          N_bins_symbiont = 0,
          N_bins_endophyte = 0,
          Total_bins = 0,
          Treatments_symbiont = "",
          Treatments_endophyte = "",
          Primarily_in = "Not in bins"
        )
      }
    } else {
      results[[ko]] <- data.frame(
        KO = ko,
        N_bins_symbiont = NA,
        N_bins_endophyte = NA,
        Total_bins = NA,
        Treatments_symbiont = "",
        Treatments_endophyte = "",
        Primarily_in = "Not found in data"
      )
    }
  }
  
  do.call(rbind, results)
}

# Analyze robust KOs
robust_kos <- unique(robust_associations$KO)
distribution_robust <- analyze_ko_distribution(robust_kos, ko_bins, bin_class)

# Diagnostic output
n_found_in_bins <- sum(distribution_robust$Total_bins > 0, na.rm = TRUE)
cat(sprintf("  Checked %d robust KOs for bin presence\n", length(robust_kos)))
cat(sprintf("  Found %d KOs in at least one bin\n", n_found_in_bins))

if (n_found_in_bins == 0) {
  cat("\n  ⚠ WARNING: No robust genes found in bins!\n")
  cat("  Checking KO column names in ko_presence_absence_table...\n")
  ko_bin_cols <- colnames(ko_bins)[!colnames(ko_bins) %in% c("bin")]
  cat(sprintf("    First 5 KO columns in bins: %s\n", 
              paste(head(ko_bin_cols, 5), collapse = ", ")))
  cat(sprintf("    First 5 robust KOs: %s\n", 
              paste(head(robust_kos, 5), collapse = ", ")))
  cat("  → Check if column names match!\n\n")
}

# Merge with robust associations
robust_with_distribution <- robust_associations %>%
  left_join(distribution_robust, by = "KO") %>%
  mutate(
    Effect_direction = ifelse(estimate > 0, "Positive", "Negative"),
    Effect_size = case_when(
      abs(estimate) > 0.5 ~ "Large",
      abs(estimate) > 0.2 ~ "Medium",
      TRUE ~ "Small"
    )
  )

write.csv(robust_with_distribution, "Robust_associations_with_distribution.csv", 
          row.names = FALSE)
cat("✓ Saved: Robust_associations_with_distribution.csv\n")

# Summary by location
location_summary <- robust_with_distribution %>%
  group_by(Primarily_in) %>%
  summarise(
    N_associations = n(),
    Genes = paste(unique(Gene), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(N_associations))

cat("\n  Distribution summary:\n")
print(location_summary)

################################################################################
# 7. CREATE VISUALIZATIONS
################################################################################

cat("\nSTEP 7: Creating visualizations...\n")

dir.create("plots", showWarnings = FALSE)

# Plot 1: Robust associations by trait
p1 <- ggplot(robust_with_distribution, 
             aes(x = reorder(Gene, estimate), y = estimate, 
                 size = R_squared, color = Category)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() +
  facet_wrap(~Trait, scales = "free") +
  scale_size_continuous(name = "R²", range = c(3, 8)) +
  scale_color_viridis_d() +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom") +
  labs(
    title = "Robust Gene-Trait Associations",
    subtitle = sprintf("%d associations significant by ≥2 methods", 
                       nrow(robust_associations)),
    x = "Gene",
    y = "Effect Size (β coefficient)"
  )

ggsave("plots/robust_associations_by_trait.svg", p1, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Plot 2: Symbiont vs Endophyte distribution
genes_in_bins <- sum(robust_with_distribution$Total_bins > 0, na.rm = TRUE)

if (genes_in_bins > 0) {
  
  plot_data <- robust_with_distribution %>%
    filter(Total_bins > 0) %>%
    mutate(Gene_label = paste0(Gene, "\n(", KO, ")"))
  
  p2 <- ggplot(plot_data, 
               aes(x = reorder(Gene_label, N_bins_symbiont))) +
    geom_col(aes(y = N_bins_symbiont), fill = "#1f77b4", alpha = 0.7) +
    geom_col(aes(y = -N_bins_endophyte), fill = "#ff7f0e", alpha = 0.7) +
    geom_hline(yintercept = 0, color = "black") +
    coord_flip() +
    scale_y_continuous(labels = abs) +
    theme_minimal(base_size = 11) +
    labs(
      title = "Distribution of Robust Genes: Symbionts vs Endophytes",
      subtitle = "Blue = Bradyrhizobium symbionts, Orange = Other bacteria (endophytes)",
      x = "Gene (KO)",
      y = "Number of Bins"
    ) +
    annotate("text", x = 1, y = max(plot_data$N_bins_symbiont, na.rm = TRUE) * 0.7, 
             label = "Symbionts →", color = "#1f77b4", fontface = "bold", size = 4) +
    annotate("text", x = 1, y = -max(plot_data$N_bins_endophyte, na.rm = TRUE) * 0.7, 
             label = "← Endophytes", color = "#ff7f0e", fontface = "bold", size = 4)
  
  ggsave("plots/symbiont_endophyte_distribution.svg", p2, 
         width = 10, height = 8, dpi = 300, bg = "white")
  
  cat(sprintf("✓ Saved symbiont/endophyte plot (%d genes found in bins)\n", genes_in_bins))
  
} else {
  cat("⚠ Skipping symbiont/endophyte plot - no robust genes found in bins\n")
}

# Plot 3: Method agreement
methods_data <- robust_with_distribution %>%
  select(Gene, Trait, Sig_standard, Sig_qvalue, Sig_trait_stratified) %>%
  mutate(Gene_Trait = paste0(Gene, " - ", Trait)) %>%
  pivot_longer(cols = starts_with("Sig_"), 
               names_to = "Method", values_to = "Significant") %>%
  mutate(
    Method = recode(Method,
                    "Sig_standard" = "Standard FDR",
                    "Sig_qvalue" = "Q-value",
                    "Sig_trait_stratified" = "Trait-Stratified"),
    Significant = replace_na(Significant, FALSE)
  )

p3 <- ggplot(methods_data, aes(x = Method, y = Gene_Trait, fill = Significant)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "gray90"),
                    labels = c("Not Sig", "Significant")) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    title = "Method Agreement for Robust Associations",
    x = "Correction Method",
    y = "Gene - Trait",
    fill = "Status"
  )

ggsave("plots/method_agreement.svg", p3, 
       width = 8, height = max(6, nrow(robust_associations) * 0.15), 
       dpi = 300, bg = "white")

# Plot 4: Complete N-Metabolism Gene Distribution (ALL genes, not just robust)
cat("\nCreating complete N-metabolism gene distribution heatmap...\n")

# Get all N-metabolism genes that are in the ko_bins data
all_n_genes <- gene_list$KEGG_KO
genes_in_bins_data <- intersect(all_n_genes, colnames(ko_bins))

if (length(genes_in_bins_data) > 0) {
  
  cat(sprintf("  Found %d N-metabolism genes in bins data\n", length(genes_in_bins_data)))
  
  # --- PREPARE BIN ANNOTATION ---
  
  extract_readable_name <- function(classification) {
    # Species
    species <- str_extract(classification, "s__([^;]+)$")
    species <- gsub("s__", "", species)
    
    # Genus
    genus <- str_extract(classification, "g__([^;]+);")
    if (is.na(genus)) genus <- str_extract(classification, "g__([^;]+)$")
    genus <- gsub("g__|;", "", genus)
    
    # If species exists and is not empty
    if (!is.na(species) && nchar(trimws(species)) > 0) {
      species_clean <- gsub("_", " ", species)
      if (grepl("^[A-Z][a-z]+\\s", species_clean)) return(species_clean)
      if (!is.na(genus) && nchar(trimws(genus)) > 0) return(paste(genus, species_clean))
      return(species_clean)
    } 
    
    # If genus exists
    if (!is.na(genus) && nchar(trimws(genus)) > 0) return(genus)
    
    # Family
    family <- str_extract(classification, "f__([^;]+);")
    family <- gsub("f__|;", "", family)
    if (!is.na(family) && nchar(trimws(family)) > 0) return(family)
    
    # Order
    order <- str_extract(classification, "o__([^;]+);")
    order <- gsub("o__|;", "", order)
    if (!is.na(order) && nchar(trimws(order)) > 0) return(order)
    
    # Class
    class <- str_extract(classification, "c__([^;]+);")
    class <- gsub("c__|;", "", class)
    if (!is.na(class) && nchar(trimws(class)) > 0) return(class)
    
    # Phylum
    phylum <- str_extract(classification, "p__([^;]+);")
    phylum <- gsub("p__|;", "", phylum)
    if (!is.na(phylum) && nchar(trimws(phylum)) > 0) return(phylum)
    
    # Domain
    domain <- str_extract(classification, "d__([^;]+);")
    domain <- gsub("d__|;", "", domain)
    if (!is.na(domain) && nchar(trimws(domain)) > 0) return(domain)
    
    return("Unknown")
  }
  
  # Map treatment codes to full names
  bin_class_complete <- bin_class %>%
    mutate(Treatment_Full = if_else(treatment %in% names(treatment_names),
                                    treatment_names[treatment],
                                    treatment))
  
  bin_annotation_complete <- bin_class_complete %>%
    mutate(
      Treatment_Full = factor(Treatment_Full, levels = treatment_order),
      Readable_name = sapply(classification, extract_readable_name),
      Readable_label = Readable_name
    ) %>%
    arrange(Treatment_Full, organism_type, genus)
  
  # --- PREPARE COMPLETE HEATMAP DATA ---
  
  ko_presence_complete <- ko_bins %>%
    filter(bin %in% bin_annotation_complete$bin) %>%
    select(bin, all_of(genes_in_bins_data)) %>%
    pivot_longer(cols = -bin, names_to = "KO", values_to = "Present") %>%
    left_join(bin_annotation_complete %>% 
                select(bin, treatment, Treatment_Full, organism_type, Readable_label), 
              by = "bin") %>%
    left_join(gene_list %>% 
                select(KEGG_KO, Gene, Category) %>% 
                distinct(KEGG_KO, .keep_all = TRUE), 
              by = c("KO" = "KEGG_KO"))
  
  # Create gene labels with KO IDs
  gene_labels_complete <- ko_presence_complete %>%
    distinct(KO, Gene, Category) %>%
    mutate(Gene_label = paste0(Gene, "\n(", KO, ")"))
  
  # Order genes by category and gene name
  gene_order <- gene_labels_complete %>%
    arrange(Category, Gene) %>%
    pull(Gene_label)
  
  ko_presence_plot_complete <- ko_presence_complete %>%
    left_join(gene_labels_complete, by = c("KO", "Gene", "Category")) %>%
    mutate(
      Gene_label = factor(Gene_label, levels = gene_order),
      Present = factor(Present, levels = c(0, 1), labels = c("Absent", "Present")),
      Treatment_Full = factor(Treatment_Full, levels = treatment_order)
    ) %>%
    arrange(Treatment_Full, organism_type, bin) %>%
    mutate(bin = factor(bin, levels = unique(bin)))
  
  # Bin label mapping
  bin_label_map_complete <- ko_presence_plot_complete %>%
    distinct(bin, Readable_label) %>%
    arrange(match(bin, levels(ko_presence_plot_complete$bin)))
  
  # --- CREATE COMPLETE HEATMAP ---
  
  p_complete <- ggplot(ko_presence_plot_complete, 
                       aes(x = bin, y = Gene_label)) +
    geom_tile(aes(fill = interaction(Present, organism_type, sep = ".")), 
              color = "white", size = 0.5) +
    scale_fill_manual(
      values = c(
        "Absent.symbiont" = "white",
        "Absent.endophyte" = "white",
        "Present.symbiont" = "#0571b0",  # Blue for symbiont
        "Present.endophyte" = "#ff7f0e"  # Orange for endophyte
      ),
      breaks = c("Present.symbiont", "Present.endophyte"),
      labels = c("Symbiont", "Endophyte"),
      name = "Gene Status"
    ) +
    facet_grid(. ~ Treatment_Full, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = setNames(bin_label_map_complete$Readable_label, 
                                       bin_label_map_complete$bin)) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray90", color = "black"),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(
      title = "Complete N-Metabolism Gene Distribution Across All Bins",
      subtitle = "Blue = Symbiont, Orange = Endophyte. Includes ALL genes present in bins, not just statistically robust associations.",
      x = "Organism",
      y = "Gene (KO)"
    )
  
  # Calculate dimensions
  n_genes_complete <- length(genes_in_bins_data)
  n_bins_complete <- nrow(bin_annotation_complete)
  
  width_complete <- max(16, n_bins_complete * 0.35 + 4)
  height_complete <- max(10, n_genes_complete * 0.25 + 3)
  
  ggsave("plots/gene_bin_heatmap_complete.svg", p_complete, 
         width = width_complete, height = height_complete, dpi = 300, bg = "white")
  ggsave("plots/gene_bin_heatmap_complete.jpg", p_complete, 
         width = width_complete, height = height_complete, dpi = 300, bg = "white")
  
  cat(sprintf("✓ Saved complete gene-bin heatmap (%d genes × %d bins)\n", 
              n_genes_complete, n_bins_complete))
  
} else {
  cat("⚠ Skipping complete gene-bin heatmap - no N-metabolism genes found in bins\n")
}

# Plot 5: Gene presence heatmap + Effect Size Barplots (ROBUST GENES ONLY)
cat("\nCreating gene-bin presence heatmap with effect size bars...\n")

# Get robust KOs that are in bins
robust_kos_in_bins <- robust_with_distribution %>%
  filter(Total_bins > 0) %>%
  pull(KO) %>%
  unique()

if (length(robust_kos_in_bins) > 0) {
  
  # --- 1. PREPARE BIN ANNOTATION ---
  
  extract_readable_name <- function(classification) {
    # Try to extract taxonomy at each level from specific to general
    # GTDB format: d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
    
    # Species
    species <- str_extract(classification, "s__([^;]+)$")
    species <- gsub("s__", "", species)
    
    # Genus
    genus <- str_extract(classification, "g__([^;]+);")
    if (is.na(genus)) genus <- str_extract(classification, "g__([^;]+)$")
    genus <- gsub("g__|;", "", genus)
    
    # If species exists and is not empty
    if (!is.na(species) && nchar(trimws(species)) > 0) {
      species_clean <- gsub("_", " ", species)
      if (grepl("^[A-Z][a-z]+\\s", species_clean)) return(species_clean)
      if (!is.na(genus) && nchar(trimws(genus)) > 0) return(paste(genus, species_clean))
      return(species_clean)
    } 
    
    # If genus exists and is not empty
    if (!is.na(genus) && nchar(trimws(genus)) > 0) {
      return(genus)
    }
    
    # Family
    family <- str_extract(classification, "f__([^;]+);")
    family <- gsub("f__|;", "", family)
    if (!is.na(family) && nchar(trimws(family)) > 0) {
      return(family)
    }
    
    # Order
    order <- str_extract(classification, "o__([^;]+);")
    order <- gsub("o__|;", "", order)
    if (!is.na(order) && nchar(trimws(order)) > 0) {
      return(order)
    }
    
    # Class
    class <- str_extract(classification, "c__([^;]+);")
    class <- gsub("c__|;", "", class)
    if (!is.na(class) && nchar(trimws(class)) > 0) {
      return(class)
    }
    
    # Phylum
    phylum <- str_extract(classification, "p__([^;]+);")
    phylum <- gsub("p__|;", "", phylum)
    if (!is.na(phylum) && nchar(trimws(phylum)) > 0) {
      return(phylum)
    }
    
    # Domain
    domain <- str_extract(classification, "d__([^;]+);")
    domain <- gsub("d__|;", "", domain)
    if (!is.na(domain) && nchar(trimws(domain)) > 0) {
      return(domain)
    }
    
    # If all else fails
    return("Unknown")
  }
  
  # Map treatment codes to full names in bin_class
  bin_class <- bin_class %>%
    mutate(Treatment_Full = if_else(treatment %in% names(treatment_names),
                                    treatment_names[treatment],
                                    treatment))
  
  bin_annotation <- bin_class %>%
    mutate(
      Treatment_Full = factor(Treatment_Full, levels = treatment_order),
      Readable_name = sapply(classification, extract_readable_name),
      Readable_label = Readable_name,  # Just use taxonomy, no (SYM)/(END) suffix
      Organism_label = case_when(
        organism_type == "symbiont" ~ paste0(genus, " (SYM)"),
        organism_type == "endophyte" ~ paste0(genus, " (END)"),
        TRUE ~ paste0(genus, " (UNK)")
      )
    ) %>%
    arrange(Treatment_Full, organism_type, genus)
  
  # --- 2. PREPARE HEATMAP DATA ---
  
  ko_presence_long <- ko_bins %>%
    filter(bin %in% bin_annotation$bin) %>%
    select(bin, all_of(robust_kos_in_bins)) %>%
    pivot_longer(cols = -bin, names_to = "KO", values_to = "Present") %>%
    left_join(bin_annotation %>% select(bin, treatment, Treatment_Full, organism_type, 
                                        Organism_label, Readable_label), 
              by = "bin") %>%
    left_join(gene_list %>% select(KEGG_KO, Gene, Category) %>% distinct(KEGG_KO, .keep_all = TRUE), 
              by = c("KO" = "KEGG_KO"))
  
  # Create gene labels
  gene_labels <- ko_presence_long %>%
    distinct(KO, Gene, Category) %>%
    mutate(
      Category_short = case_when(
        str_detect(Category, "Core Nitrogenase") ~ "Core Nif",
        str_detect(Category, "Cofactor") ~ "Nif Cofactor",
        str_detect(Category, "Denitrification") ~ "Denitrif",
        str_detect(Category, "Ammonia") ~ "NH4 Transport",
        TRUE ~ Category
      ),
      Gene_label = paste0(Gene, "\n(", KO, ")")
    )
  
  ko_presence_plot <- ko_presence_long %>%
    left_join(gene_labels, by = c("KO", "Gene", "Category")) %>%
    mutate(
      Gene_label = factor(Gene_label, levels = gene_labels$Gene_label),
      Present = factor(Present, levels = c(0, 1), labels = c("Absent", "Present")),
      Treatment_Full = factor(Treatment_Full, levels = treatment_order)
    ) %>%
    arrange(Treatment_Full, organism_type, bin) %>%
    mutate(bin = factor(bin, levels = unique(bin)))
  
  # Calculate endophyte-only
  endophyte_only_genes <- ko_presence_plot %>%
    group_by(Treatment_Full, KO, Gene_label) %>%
    summarise(
      has_symbiont = any(Present == "Present" & organism_type == "symbiont"),
      has_endophyte = any(Present == "Present" & organism_type == "endophyte"),
      .groups = "drop"
    ) %>%
    filter(has_endophyte & !has_symbiont) %>%
    mutate(is_endophyte_only = TRUE)
  
  ko_presence_plot_with_flag <- ko_presence_plot %>%
    left_join(endophyte_only_genes %>% select(Treatment_Full, KO, is_endophyte_only),
              by = c("Treatment_Full", "KO")) %>%
    mutate(is_endophyte_only = replace_na(is_endophyte_only, FALSE))
  
  # --- 3. PREPARE BARPLOT DATA WITH SIGNIFICANCE FOR ALL TRAITS ---
  
  # Get correlation data for ALL traits using trait-stratified p-values
  correlation_annotation <- robust_with_distribution %>%
    filter(KO %in% robust_kos_in_bins) %>%
    select(KO, Gene, Trait, estimate, p.value, p_adj_trait_stratified) %>%
    mutate(
      sig_symbol = case_when(
        p_adj_trait_stratified < 0.001 ~ "***",
        p_adj_trait_stratified < 0.01 ~ "**", 
        p_adj_trait_stratified < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    left_join(gene_labels %>% select(KO, Gene_label), by = "KO") %>%
    mutate(
      Gene_label = factor(Gene_label, levels = levels(ko_presence_plot_with_flag$Gene_label)),
      Direction = ifelse(estimate > 0, "Positive", "Negative")
    )
  
  # --- 4. CREATE PLOTS ---
  
  # A. Bin Label Mapping
  bin_label_map <- ko_presence_plot %>%
    distinct(bin, Readable_label) %>%
    arrange(match(bin, levels(ko_presence_plot$bin)))
  
  # B. Main Heatmap
  p4_enhanced_main <- ggplot(ko_presence_plot_with_flag, 
                             aes(x = bin, y = Gene_label)) +
    geom_tile(aes(fill = interaction(Present, organism_type, sep = ".")), 
              color = "white", size = 0.5) +
    # Add black borders for endophyte-only genes
    geom_tile(data = filter(ko_presence_plot_with_flag, 
                            is_endophyte_only & Present == "Present" & organism_type == "endophyte"),
              aes(x = bin, y = Gene_label),
              fill = NA, color = "black", linewidth = 1.2) +
    scale_fill_manual(
      values = c(
        "Absent.symbiont" = "white",
        "Absent.endophyte" = "white",
        "Present.symbiont" = "#0571b0",  # Blue for symbiont
        "Present.endophyte" = "#ff7f0e"  # Orange for endophyte
      ),
      breaks = c("Present.symbiont", "Present.endophyte"),
      labels = c("Symbiont", "Endophyte"),
      name = "Gene Status"
    ) +
    facet_grid(. ~ Treatment_Full, scales = "free_x", space = "free_x") +
    scale_x_discrete(labels = setNames(bin_label_map$Readable_label, bin_label_map$bin)) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray90", color = "black"),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "bottom",
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(
      x = "Organism",
      y = "Gene (KO)"
    )
  
  # C. Diverging Bar Plots (Effect Sizes) - WITH SIGNIFICANCE FOR ALL GENES
  p_effect_bars <- ggplot(correlation_annotation, aes(x = estimate, y = Gene_label, fill = Direction)) +
    geom_col(width = 0.7) + 
    geom_vline(xintercept = 0, color = "black", size = 0.5) +
    # Add significance stars for ALL genes
    geom_text(aes(label = sig_symbol, 
                  x = estimate,
                  hjust = ifelse(estimate > 0, -0.2, 1.2)), 
              vjust = 0.75, size = 3.5, fontface = "bold") +
    # Use manual discrete colors instead of gradient
    scale_fill_manual(
      values = c("Negative" = "#d01c8b", "Positive" = "#4dac26"),
      name = "Effect Direction",
      na.value = "transparent"
    ) +
    facet_wrap(~ Trait, nrow = 1, scales = "free_x") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(), # Hide Y labels (shared with heatmap)
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "#e6f2ff", color = "black"),
      legend.position = "bottom"
    ) +
    labs(x = "Effect Size (β)")
  
  # --- 5. COMBINE PLOTS ---
  
  p_combined <- p4_enhanced_main + p_effect_bars + 
    plot_layout(widths = c(3, 1.2), guides = "collect") +
    plot_annotation(
      title = "Robust N-Metabolism Genes: Distribution and Trait Associations",
      subtitle = paste0(
        "LEFT: Gene presence (Blue=Symbiont, Orange=Endophyte). Black border = Endophyte-only in treatment.\n",
        "RIGHT: Effect size barplots. Bars extend Right (Green) for positive correlation, Left (Pink) for negative.\n",
        "Stars indicate trait-stratified significance: *p<0.05, **p<0.01, ***p<0.001"
      ),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10),
        legend.position = "bottom"
      )
    )
  
  # Calculate dimensions
  n_genes <- length(robust_kos_in_bins)
  n_bins <- nrow(bin_annotation)
  
  width <- max(14, n_bins * 0.4 + 5)
  height <- max(8, n_genes * 0.4 + 3)
  
  ggsave("plots/gene_bin_heatmap_enhanced.svg", p_combined, 
         width = 22, height = height, dpi = 300, bg = "white")
  ggsave("plots/gene_bin_heatmap_enhanced.jpg", p_combined, 
         width = 22, height = height, dpi = 300, bg = "white")
  cat("✓ Saved enhanced gene-bin heatmap with effect size barplots\n")
  
} else {
  cat("⚠ Skipping gene-bin heatmap - no robust genes found in bins\n")
}

# Plot 6: Tier1-style ComplexHeatmap for Robust N-Metabolism Genes
cat("\nCreating Tier1-style ComplexHeatmap for robust N-metabolism genes...\n")

if (length(robust_kos_in_bins) > 0) {
  
  cat(sprintf("  Generating heatmap for %d robust N-metabolism genes\n", 
              length(robust_kos_in_bins)))
  
  # --- 1. PREPARE DATA FOR COMPLEXHEATMAP ---
  
  # Get presence data for robust KOs
  ko_presence_for_heatmap <- ko_bins %>%
    filter(bin %in% bin_annotation$bin) %>%
    select(bin, all_of(robust_kos_in_bins)) %>%
    pivot_longer(cols = -bin, names_to = "KO", values_to = "present") %>%
    filter(present == 1)
  
  # Add treatment and organism type
  ko_presence_with_meta <- ko_presence_for_heatmap %>%
    left_join(bin_annotation %>% 
                select(bin, treatment, Treatment_Full, organism_type), 
              by = "bin")
  
  # Summarize presence per KO × treatment × organism_type
  pres_summary <- ko_presence_with_meta %>%
    distinct(KO, treatment, organism_type, bin) %>%
    count(KO, treatment, organism_type, name = "n_bins")
  
  pres_wide_heatmap <- pres_summary %>%
    mutate(present_flag = 1L) %>%
    select(-n_bins) %>%
    pivot_wider(
      names_from = organism_type,
      values_from = present_flag,
      values_fill = 0L
    )
  
  # Create complete KO × treatment grid
  treat_short <- c("ces", "RH", "carR", "carK", "hok", "mtz")
  
  grid_heatmap <- expand_grid(
    KO = robust_kos_in_bins,
    treatment = treat_short
  )
  
  pres_full_heatmap <- grid_heatmap %>%
    left_join(pres_wide_heatmap, by = c("KO", "treatment")) %>%
    mutate(
      endophyte = replace_na(endophyte, 0L),
      symbiont = replace_na(symbiont, 0L),
      State = case_when(
        endophyte == 1L & symbiont == 1L ~ "Both",
        endophyte == 1L & symbiont == 0L ~ "Endophyte_only",
        endophyte == 0L & symbiont == 1L ~ "Symbiont_only",
        TRUE ~ "Absent"
      )
    )
  
  # --- 2. ADD METADATA ---
  
  # Get gene annotations - handle multiple traits per KO
  ko_meta_heatmap <- robust_with_distribution %>%
    filter(KO %in% robust_kos_in_bins) %>%
    # For KOs with multiple traits, keep the one with smallest p-value
    group_by(KO) %>%
    slice_min(p_adj_trait_stratified, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    distinct(KO, Gene, Trait, Category) %>%
    left_join(gene_list %>% 
                distinct(KEGG_KO, Category, .keep_all = TRUE) %>%
                select(KEGG_KO, Category), 
              by = c("KO" = "KEGG_KO"),
              suffix = c("_robust", "_list"))
  
  # Define trait labels
  trait_labels_heatmap <- c(
    "Plant_Biomass" = "Plant\nbiomass (g)",
    "percent_N" = "Leaf N\n(%)",
    "NMF" = "Nodule mass\nfraction",
    "Fixation_per_Nodule" = "Fixation rate\nper nodule"
  )
  
  # Add metadata to presence data
  pres_full_heatmap <- pres_full_heatmap %>%
    mutate(Treatment_Full = treatment_names[treatment]) %>%
    left_join(ko_meta_heatmap, by = "KO", relationship = "many-to-one") %>%
    mutate(
      Trait = factor(Trait, 
                     levels = c("Plant_Biomass", "percent_N", 
                                "NMF", "Fixation_per_Nodule")),
      Trait_label = trait_labels_heatmap[as.character(Trait)],
      
      # Define functional categories based on gene category
      Functional_Category = case_when(
        !is.na(Category_list) & str_detect(Category_list, "Nitrogenase|Nif") ~ "Nitrogenase & N-fixation",
        !is.na(Category_list) & str_detect(Category_list, "Denitrification") ~ "Denitrification",
        !is.na(Category_list) & str_detect(Category_list, "Ammonia|Nitrate|Nitrite") ~ "Inorganic N metabolism",
        !is.na(Category_list) & str_detect(Category_list, "Assimilation|Transport") ~ "N transport & assimilation",
        !is.na(Category_robust) & str_detect(Category_robust, "Nitrogenase|Nif") ~ "Nitrogenase & N-fixation",
        !is.na(Category_robust) & str_detect(Category_robust, "Denitrification") ~ "Denitrification",
        !is.na(Category_robust) & str_detect(Category_robust, "Ammonia|Nitrate|Nitrite") ~ "Inorganic N metabolism",
        !is.na(Category_robust) & str_detect(Category_robust, "Assimilation|Transport") ~ "N transport & assimilation",
        TRUE ~ "Other N-metabolism"
      )
    )
  
  # --- 3. ORDER GENES ---
  
  # Get unique KO order (one entry per KO, ordered by primary trait assignment)
  ko_order_heatmap <- pres_full_heatmap %>%
    distinct(KO, Gene, Trait, Functional_Category) %>%
    arrange(Trait, Functional_Category, Gene) %>%
    distinct(KO, .keep_all = TRUE) %>%  # Ensure unique KOs
    pull(KO)
  
  # Verify no duplicates
  if (any(duplicated(ko_order_heatmap))) {
    cat("  WARNING: Duplicate KOs found, removing duplicates\n")
    ko_order_heatmap <- unique(ko_order_heatmap)
  }
  
  cat(sprintf("  Ordered %d unique KOs for heatmap\n", length(ko_order_heatmap)))
  
  pres_full_heatmap <- pres_full_heatmap %>%
    mutate(
      KO = factor(KO, levels = ko_order_heatmap),
      State = factor(State, 
                     levels = c("Absent", "Symbiont_only", 
                                "Endophyte_only", "Both")),
      Treatment_Full = factor(Treatment_Full, levels = treatment_order)
    )
  
  # --- 4. BUILD MATRIX ---
  
  mat_state_df_heatmap <- pres_full_heatmap %>%
    select(Treatment_Full, KO, State) %>%
    arrange(Treatment_Full, KO) %>%
    pivot_wider(names_from = KO, values_from = State)
  
  state_mat_heatmap <- as.matrix(mat_state_df_heatmap[, -1])
  rownames(state_mat_heatmap) <- mat_state_df_heatmap$Treatment_Full
  
  # --- 5. CREATE COLUMN ANNOTATIONS ---
  
  # Get unique annotation data (one row per KO)
  col_annot_df_heatmap <- pres_full_heatmap %>%
    distinct(KO, Gene, Trait_label, Functional_Category) %>%
    # Keep only the first instance of each KO (matches ko_order_heatmap)
    group_by(KO) %>%
    slice(1) %>%
    ungroup() %>%
    arrange(KO)
  
  # Ensure annotation order matches matrix columns
  col_annot_df_heatmap <- col_annot_df_heatmap[
    match(colnames(state_mat_heatmap), col_annot_df_heatmap$KO), 
  ]
  
  # Verify match
  if (!all(col_annot_df_heatmap$KO == colnames(state_mat_heatmap))) {
    warning("Column annotation KOs don't match matrix columns")
  }
  
  # Use Gene names as column labels
  colnames(state_mat_heatmap) <- col_annot_df_heatmap$Gene
  
  # --- 6. DEFINE COLORS ---
  
  # State colors (matching Tier1 style)
  state_cols_heatmap <- c(
    "Absent" = "#f0f0f0",          # light grey
    "Symbiont_only" = "#8DD3C7",   # cyan
    "Endophyte_only" = "#FFFFB3",  # yellow
    "Both" = "#BEBADA"             # purple
  )
  
  # Functional category colors
  cat_cols_heatmap <- c(
    "Nitrogenase & N-fixation" = "#FB8072",        # red/orange
    "Denitrification" = "#FDB462",                 # orange
    "Inorganic N metabolism" = "#B3DE69",          # yellow-green
    "N transport & assimilation" = "#80B1D3",      # blue
    "Other N-metabolism" = "#D9D9D9"               # grey
  )
  
  # Primary trait colors for top annotation
  primary_trait_cols <- c(
    "Leaf N\n(%)"             = "#66C2A5",
    "Nodule mass\nfraction"   = "#FC8D62",
    "Plant\nbiomass (g)"      = "#8DA0CB",
    "Fixation rate\nper nodule" = "#A6D854"
  )

  # Create top annotation — Functional Category + Primary Trait indicator
  ha_top_heatmap <- HeatmapAnnotation(
    "Primary trait" = col_annot_df_heatmap$Trait_label,
    "Functional Category" = col_annot_df_heatmap$Functional_Category,
    col = list(
      "Primary trait" = primary_trait_cols,
      "Functional Category" = cat_cols_heatmap
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10),
    height = unit(4, "mm")
  )
  
  # --- 7. CREATE HEATMAP WITH BARPLOT ANNOTATIONS ---
  
  # Prepare effect size data for annotations - one value per gene per trait
  # IMPORTANT: Include ALL trait associations, not just primary
  correlation_for_barplot <- robust_with_distribution %>%
    filter(KO %in% robust_kos_in_bins) %>%
    select(KO, Gene, Trait, estimate, p_adj_trait_stratified) %>%
    mutate(
      sig_symbol = case_when(
        p_adj_trait_stratified < 0.001 ~ "***",
        p_adj_trait_stratified < 0.01 ~ "**",
        p_adj_trait_stratified < 0.05 ~ "*",
        TRUE ~ ""
      ),
      # Convert Trait to Trait_label for matching
      Trait_label = case_when(
        Trait == "Plant_Biomass" ~ "Plant\nbiomass (g)",
        Trait == "percent_N" ~ "Leaf N\n(%)",
        Trait == "NMF" ~ "Nodule mass\nfraction",
        Trait == "Fixation_per_Nodule" ~ "Fixation rate\nper nodule",
        TRUE ~ Trait
      )
    )
  
  # DIAGNOSTIC: Check for genes with multiple trait associations
  multi_trait_genes <- correlation_for_barplot %>%
    group_by(Gene) %>%
    summarise(n_traits = n_distinct(Trait)) %>%
    filter(n_traits > 1)
  
  if (nrow(multi_trait_genes) > 0) {
    cat(sprintf("\n  Found %d genes with multiple trait associations:\n", nrow(multi_trait_genes)))
    for (i in 1:min(3, nrow(multi_trait_genes))) {
      gene <- multi_trait_genes$Gene[i]
      traits <- correlation_for_barplot %>%
        filter(Gene == gene) %>%
        pull(Trait) %>%
        paste(collapse = ", ")
      cat(sprintf("    - %s: %s\n", gene, traits))
    }
  }
  
  # Get gene order from heatmap (this is based on PRIMARY trait assignment)
  genes_in_order <- colnames(state_mat_heatmap)
  
  # Create data for bottom barplot annotations (one per trait)
  trait_list <- unique(col_annot_df_heatmap$Trait_label)
  
  # Create vectors for each trait - showing ALL associations
  trait_effects <- list()
  trait_colors <- list()
  
  for (trait_name in trait_list) {
    # Get ALL genes associated with this trait (not just primary)
    trait_data <- correlation_for_barplot %>%
      filter(Trait_label == trait_name)
    
    cat(sprintf("  Trait '%s': %d gene associations\n", trait_name, nrow(trait_data)))
    
    # Create vector matching all genes in order
    effect_vec <- rep(0, length(genes_in_order))
    names(effect_vec) <- genes_in_order
    
    # Fill in values for genes that have this trait association
    for (i in 1:nrow(trait_data)) {
      gene_name <- trait_data$Gene[i]
      if (gene_name %in% genes_in_order) {
        effect_vec[gene_name] <- trait_data$estimate[i]
      }
    }
    
    # Create color vector (green for positive, pink for negative)
    color_vec <- ifelse(effect_vec > 0, "#4CAF50", "#E91E63")
    color_vec[effect_vec == 0] <- "white"  # No color for absent values
    
    trait_effects[[trait_name]] <- effect_vec
    trait_colors[[trait_name]] <- color_vec
  }
  
  # Create bottom annotation with barplots for each trait
  # Build annotations dynamically - EACH TRAIT GETS ITS OWN SCALE
  ha_bottom_list <- list()
  
  # ============================================================================
  # MANUAL Y-AXIS RANGES - SET YOUR OWN VALUES HERE
  # ============================================================================
  # Option 1: Set specific ranges for each trait
  manual_ranges <- list(
    "Leaf N\n(%)" = c(-0.6, 1.2),              # Change these values as needed
    "Nodule mass\nfraction" = c(-0.03, 0.03),  # Change these values as needed
    "Plant\nbiomass (g)" = NULL,           # NULL = auto-calculate
    "Fixation rate\nper nodule" = NULL     # NULL = auto-calculate
  )
  
  # Option 2: Use one range for ALL traits (uncomment to use)
  # manual_ranges <- c(-1, 1)  # Same range for all traits
  # ============================================================================
  
  for (i in seq_along(trait_list)) {
    trait_name <- trait_list[i]
    
    # Check if manual range is specified for this trait
    if (is.list(manual_ranges) && !is.null(manual_ranges[[trait_name]])) {
      # Use manual range
      axis_range_trait <- manual_ranges[[trait_name]]
      cat(sprintf("  Using manual range for %s: [%.2f, %.2f]\n", 
                  trait_name, axis_range_trait[1], axis_range_trait[2]))
    } else if (!is.list(manual_ranges) && length(manual_ranges) == 2) {
      # Use single manual range for all
      axis_range_trait <- manual_ranges
    } else {
      # Auto-calculate scale for THIS trait only
      trait_values <- trait_effects[[trait_name]]
      trait_nonzero <- trait_values[trait_values != 0]
      
      if (length(trait_nonzero) > 0) {
        max_abs_effect_trait <- max(abs(trait_nonzero))
        axis_range_trait <- c(-max_abs_effect_trait * 1.15, max_abs_effect_trait * 1.15)
      } else {
        # If no values, use default
        axis_range_trait <- c(-1, 1)
      }
      cat(sprintf("  Auto-calculated range for %s: [%.2f, %.2f]\n", 
                  trait_name, axis_range_trait[1], axis_range_trait[2]))
    }
    
    ha_bottom_list[[trait_name]] <- anno_barplot(
      trait_effects[[trait_name]],
      bar_width = 0.9,
      gp = gpar(fill = trait_colors[[trait_name]]),
      baseline = 0,  # Draw baseline at zero
      baseline_gp = gpar(lty = 2, col = "black", lwd = 1),  # Dashed line
      axis = TRUE,
      axis_param = list(
        side = "left",
        at = seq(axis_range_trait[1], axis_range_trait[2], length.out = 5),
        labels_rot = 0,
        gp = gpar(fontsize = 8)
      ),
      ylim = axis_range_trait,
      height = unit(4, "cm")
    )
  }
  
  # Combine into HeatmapAnnotation
  ha_bottom <- do.call(HeatmapAnnotation, c(
    ha_bottom_list,
    list(
      annotation_name_side = "left",
      annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
      annotation_name_rot = 0,
      gap = unit(2, "mm")
    )
  ))
  
  cell_side_heatmap <- unit(6, "mm")
  ht_width_heatmap <- cell_side_heatmap * ncol(state_mat_heatmap)
  ht_height_heatmap <- cell_side_heatmap * nrow(state_mat_heatmap)
  
  ht_complexheatmap <- Heatmap(
    state_mat_heatmap,
    name = "Gene status",
    col = state_cols_heatmap,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_side = "bottom",
    column_names_rot = 90,
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_gp = gpar(fontsize = 9),
    top_annotation = ha_top_heatmap,
    bottom_annotation = ha_bottom,  # ADD BARPLOTS AT BOTTOM
    width = ht_width_heatmap,
    height = ht_height_heatmap,
    row_names_side = "left",
    border = TRUE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  # Calculate output dimensions (taller to accommodate larger barplots)
  barplot_total_height <- length(trait_list) * 4 + (length(trait_list) - 1) * 0.2  # 4cm per trait + gaps
  calc_h_heatmap <- convertHeight(ht_height_heatmap + unit(barplot_total_height + 8, "cm"), 
                                  "in", valueOnly = TRUE)
  calc_w_heatmap <- convertWidth(ht_width_heatmap + unit(14, "cm"), 
                                 "in", valueOnly = TRUE)
  
  # --- 8. SAVE HEATMAP WITH BARPLOTS ---
  
  # SVG
  svglite("plots/nitrogen_genes_complexheatmap_with_bars.svg", 
          width = calc_w_heatmap, height = calc_h_heatmap)
  draw(ht_complexheatmap, 
       merge_legend = TRUE, 
       padding = unit(c(24, 2, 12, 8), "mm"))  # Increased bottom and left padding
  grid.text("Plant Traits", 
            x = unit(25, "mm"), 
            y = unit(1, "npc") - unit(5, "mm"), 
            just = "left", 
            gp = gpar(fontsize = 10, fontface = "bold"))
  dev.off()
  
  # PNG
  png("plots/nitrogen_genes_complexheatmap_with_bars.png", 
      width = calc_w_heatmap, height = calc_h_heatmap, 
      units = "in", res = 300)
  draw(ht_complexheatmap, 
       merge_legend = TRUE, 
       padding = unit(c(24, 2, 12, 8), "mm"))  # Increased bottom and left padding
  grid.text("Plant Traits", 
            x = unit(25, "mm"), 
            y = unit(1, "npc") - unit(5, "mm"), 
            just = "left", 
            gp = gpar(fontsize = 10, fontface = "bold"))
  dev.off()
  
  # PDF
  pdf("plots/nitrogen_genes_complexheatmap_with_bars.pdf", 
      width = calc_w_heatmap, height = calc_h_heatmap)
  draw(ht_complexheatmap, 
       merge_legend = TRUE, 
       padding = unit(c(24, 2, 12, 8), "mm"))  # Increased bottom and left padding
  grid.text("Plant Traits", 
            x = unit(25, "mm"), 
            y = unit(1, "npc") - unit(5, "mm"), 
            just = "left", 
            gp = gpar(fontsize = 10, fontface = "bold"))
  dev.off()
  
  cat("✓ Saved ComplexHeatmap with effect size barplots:\n")
  cat("  - plots/nitrogen_genes_complexheatmap_with_bars.svg\n")
  cat("  - plots/nitrogen_genes_complexheatmap_with_bars.png\n")
  cat("  - plots/nitrogen_genes_complexheatmap_with_bars.pdf\n")
  
  
} else {
  cat("⚠ Skipping ComplexHeatmap - no robust genes found in bins\n")
}

# Count figures saved  
n_figures <- 2 + ifelse(genes_in_bins > 0, 1, 0) + ifelse(length(genes_in_bins_data) > 0, 1, 0) + ifelse(length(robust_kos_in_bins) > 0, 1, 0)
cat(sprintf("\n✓ Saved %d figures to plots/ directory\n", n_figures))
################################################################################
# 8. GENERATE SUMMARY REPORT
################################################################################

cat("\nSTEP 8: Generating summary report...\n")

report <- paste0(
  "================================================================================\n",
  "NITROGEN METABOLISM GENE-TRAIT ASSOCIATION ANALYSIS\n",
  "Summary Report\n",
  "================================================================================\n\n",
  "Analysis Date: ", Sys.Date(), "\n\n",
  "INPUT DATA:\n",
  "- Samples analyzed: ", nrow(data_merged), "\n",
  "- Nitrogen metabolism genes tested: ", length(ko_cols), "\n",
  "- Plant traits: ", paste(trait_cols, collapse = ", "), "\n",
  "- Total tests: ", nrow(ko_results), "\n\n",
  "MULTIPLE TESTING CORRECTION RESULTS:\n",
  "- Standard BH-FDR (< 0.05): ", sum(ko_results$Sig_standard), "\n",
  "- Q-value (< 0.05): ", sum(ko_results$Sig_qvalue, na.rm = TRUE), "\n"
)

if (exists("qobj")) {
  report <- paste0(report,
                   "  π₀ (true nulls): ", sprintf("%.3f", qobj$pi0), " (", 
                   sprintf("%.1f%%", 100 * (1 - qobj$pi0)), " true effects)\n"
  )
}

report <- paste0(report,
                 "- Trait-stratified FDR (< 0.05): ", sum(ko_results$Sig_trait_stratified), "\n"
)

if (!all(is.na(ko_results$p_adj_meff))) {
  report <- paste0(report,
                   "- Meff-adjusted FDR (< 0.05): ", sum(ko_results$Sig_meff, na.rm = TRUE), "\n"
  )
}

report <- paste0(report,
                 "\nROBUST ASSOCIATIONS (≥2 methods):\n",
                 "- Total: ", nrow(robust_associations), "\n",
                 "- By trait:\n"
)

trait_counts <- robust_associations %>% count(Trait) %>% arrange(desc(n))
for (i in 1:nrow(trait_counts)) {
  report <- paste0(report, "  ", trait_counts$Trait[i], ": ", trait_counts$n[i], "\n")
}

report <- paste0(report,
                 "\nFUNCTIONAL CATEGORIES:\n"
)

category_counts <- robust_associations %>% count(Category) %>% arrange(desc(n))
for (i in 1:min(5, nrow(category_counts))) {
  report <- paste0(report, "  ", category_counts$Category[i], ": ", 
                   category_counts$n[i], "\n")
}

report <- paste0(report,
                 "\nSYMBIONT VS ENDOPHYTE DISTRIBUTION:\n"
)

for (i in 1:nrow(location_summary)) {
  report <- paste0(report, "  ", location_summary$Primarily_in[i], ": ", 
                   location_summary$N_associations[i], " associations\n")
}

report <- paste0(report,
                 "\nKEY FINDINGS:\n",
                 "1. Q-value analysis revealed that ", 
                 ifelse(exists("qobj"), 
                        sprintf("%.0f%%", 100 * (1 - qobj$pi0)), 
                        "many"),
                 " of tests represent true associations\n",
                 "2. Alternative corrections identified ", 
                 nrow(robust_associations) - sum(ko_results$Sig_standard),
                 " additional associations missed by standard FDR\n",
                 "3. Most associations involve ",
                 trait_counts$Trait[1], " (", trait_counts$n[1], "/", nrow(robust_associations), ")\n",
                 "4. Distribution analysis shows genes primarily in ",
                 ifelse(nrow(location_summary) > 0, 
                        tolower(location_summary$Primarily_in[1]), 
                        "various locations"), "\n\n",
                 "OUTPUT FILES:\n",
                 "- KO_results_all_methods.csv (all results with all correction methods)\n",
                 "- Robust_associations.csv (associations significant by ≥2 methods)\n",
                 "- Robust_associations_with_distribution.csv (with symbiont/endophyte info)\n",
                 "- plots/robust_associations_by_trait.svg\n",
                 "- plots/symbiont_endophyte_distribution.svg\n",
                 "- plots/method_agreement.svg\n",
                 "- plots/gene_bin_heatmap_complete.svg (ALL N-metabolism genes across all bins)\n",
                 "- plots/nitrogen_genes_complexheatmap_with_bars.svg (ComplexHeatmap + barplots)\n",
                 "- plots/gene_bin_heatmap_enhanced.svg (robust genes with effect size bars)\n\n",
                 "================================================================================\n"
)

cat(report)
writeLines(report, "Analysis_Summary_Report.txt")
cat("\n✓ Saved: Analysis_Summary_Report.txt\n")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("================================================================================\n")
cat("\nKey Output Files:\n")
cat("  Data:\n")
cat("    - Robust_associations_with_distribution.csv (complete results with bin info)\n")
cat("    - KO_results_all_methods.csv (all genes with all correction methods)\n")
cat("  Figures:\n")
cat("    - plots/gene_bin_heatmap_complete.svg (ALL N-genes in all bins)\n")
cat("    - plots/nitrogen_genes_complexheatmap_with_bars.svg (ComplexHeatmap + barplots)\n")
cat("    - plots/gene_bin_heatmap_enhanced.svg (robust genes with trait correlations)\n")
cat("    - plots/robust_associations_by_trait.svg (effect sizes by trait)\n")
cat("    - plots/symbiont_endophyte_distribution.svg (distribution comparison)\n")
cat("    - plots/method_agreement.svg (correction method comparison)\n\n")