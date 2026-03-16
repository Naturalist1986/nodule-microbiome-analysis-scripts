#!/usr/bin/env Rscript

# Combined Four-Panel Figure - EGG VERSION
# Layout: A (square, top left) + B (top middle-right) in upper third
#         C (bottom left) + D (bottom right) in lower two-thirds
# Author: Moshe
# Date: 2025-02-16

# Check and install required packages
required_packages <- c("readxl", "ggplot2", "dplyr", "tidyr", "egg", "grid", "png")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load required libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(egg)
library(grid)
library(png)

cat("================================================================================\n")
cat("CREATING COMBINED FOUR-PANEL FIGURE (EGG METHOD)\n")
cat("================================================================================\n\n")

# ============================================================================
# 1. LOAD EXISTING FIGURES
# ============================================================================

cat("Step 1: Loading existing figures...\n")

mantel_img <- png::readPNG("mantel_linkET_figure.png")
mantel_grob <- rasterGrob(mantel_img, interpolate = TRUE)

traits_img <- png::readPNG("plant_traits_final_formatted.png")
traits_grob <- rasterGrob(traits_img, interpolate = TRUE)

ordinations_img <- png::readPNG("original_ordinations_combined.png")
ordinations_grob <- rasterGrob(ordinations_img, interpolate = TRUE)

cat("  Loaded all PNG files\n\n")

# ============================================================================
# 2. CREATE PERCENT NODULATION PLOT (Panel A - square, top left)
# ============================================================================

cat("Step 2: Creating percent nodulation plot...\n")

# Consistent named color palette (Set2 assigned alphabetically by display name)
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

data <- read_excel("Inoculation_Percent.xlsx")

# Map raw treatment codes to display names
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

inoculation_summary <- data %>%
  rename(percent_inoculated = Inoculation_Sucess) %>%
  select(treatment, percent_inoculated) %>%
  filter(!is.na(treatment)) %>%
  mutate(treatment_display = dplyr::recode(treatment, !!!treatment_display_map)) %>%
  mutate(treatment_display = factor(treatment_display, levels = sort(unique(treatment_display))))

panel_a <- ggplot(inoculation_summary,
                  aes(x = treatment_display,
                      y = percent_inoculated,
                      fill = treatment_display)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", percent_inoculated)),
            vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = site_colors) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  labs(
    y = "Nodulation success rate (%)",
    title = "Inoculation Success"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "bold", size = 11),
    axis.text.y = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.title.y = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 5, 2, 5)
  ) +
  annotate("text", x = -Inf, y = Inf, label = "A", size = 5, fontface = "bold",
           hjust = -0.2, vjust = 1.8)

cat("  Percent inoculation plot created\n\n")

# ============================================================================
# 3. CREATE PANELS B, C, D WITH LABELS
# ============================================================================

cat("Step 3: Creating other panels with labels...\n")

# Panel B: Plant traits
panel_b <- ggplot() +
  annotation_custom(traits_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotate("text", x = 0.01, y = 0.98, label = "B", size = 5, fontface = "bold",
           hjust = 0, vjust = 1) +
  theme_void() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")

# Panel C: Mantel test
panel_c <- ggplot() +
  annotation_custom(mantel_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotate("text", x = 0.02, y = 0.98, label = "C", size = 5, fontface = "bold",
           hjust = 0, vjust = 1) +
  theme_void() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")

# Panel D: Procrustes ordinations
panel_d <- ggplot() +
  annotation_custom(ordinations_grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  annotate("text", x = 0.01, y = 0.98, label = "D", size = 5, fontface = "bold",
           hjust = 0, vjust = 1) +
  theme_void() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off")

cat("  All panels created\n\n")

# ============================================================================
# 4. COMBINE PANELS USING EGG (with gridExtra for layout)
# ============================================================================

cat("Step 4: Combining panels with egg...\n")

# Egg's ggarrange is better suited for simple layouts
# For complex layouts like this, we'll use gridExtra which egg depends on
# Convert to grobs
grob_a <- ggplotGrob(panel_a)
grob_b <- ggplotGrob(panel_b)
grob_c <- ggplotGrob(panel_c)
grob_d <- ggplotGrob(panel_d)

# Create layout matrix
layout_matrix <- rbind(c(1, 2, 2),
                       c(3, 4, 4),
                       c(3, 4, 4))

# Use gridExtra's arrangeGrob (egg extends this)
# Heights: top row (A+B) slightly shorter, bottom two rows (C+D) taller
final_grob <- gridExtra::arrangeGrob(grob_a, grob_b, grob_c, grob_d,
                                     layout_matrix = layout_matrix,
                                     heights = c(1.1, 1, 1),
                                     widths = c(0.9, 1, 1))

cat("  Panels combined successfully\n\n")

# ============================================================================
# 5. SAVE FIGURE
# ============================================================================

cat("Step 5: Saving combined figure...\n")

# Save as PNG
png("combined_four_panel_egg.png", width = 16, height = 12,
    units = "in", res = 300, bg = "white")
grid.draw(final_grob)
dev.off()
cat("  Saved: combined_four_panel_egg.png\n")

# Save as PDF with editable text
cairo_pdf("combined_four_panel_egg.pdf", width = 16, height = 12, bg = "white")
grid.draw(final_grob)
dev.off()
cat("  Saved: combined_four_panel_egg.pdf (editable text)\n\n")

# ============================================================================
# 6. SUMMARY
# ============================================================================

cat("================================================================================\n")
cat("FIGURE CREATION COMPLETE (EGG METHOD)\n")
cat("================================================================================\n\n")

cat("Panel layout:\n")
cat("  Top row (upper 1/3):    A (square, left) + B (wide, middle-right)\n")
cat("  Bottom row (lower 2/3): C (left) + D (right)\n\n")

cat("Output files created:\n")
cat("  - combined_four_panel_egg.png (300 dpi, 16×14 inches)\n")
cat("  - combined_four_panel_egg.pdf (editable text)\n\n")

cat("Done!\n")
