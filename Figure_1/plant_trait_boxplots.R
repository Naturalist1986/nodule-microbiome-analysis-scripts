#!/usr/bin/env Rscript

# Load libraries
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(car)
library(multcomp)
library(emmeans)
library(multcompView)
library(broom)

# 1. DATA PREPARATION
data <- read_excel("Guy_measurements.xlsx")

treatment_map <- c("carK" = "CAR-C", "carR" = "CAR-G",
                   "ces" = "CES", "hok" = "HUK",
                   "mtz" = "MTZ", "RH" = "RH")

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

data <- data %>%
  mutate(
    Treatment_Full = factor(Treatment, levels = names(treatment_map), labels = treatment_map),
    NMF_percent = NMF * 100
  )

# Configuration for the four traits
traits_config <- list(
  list(col = "percent_N", title = "Leaf N", ylab = "Leaf N (%)", tag = "(a)"),
  list(col = "Plant_Biomass", title = "Plant Biomass", ylab = "Plant Biomass (g)", tag = "(b)"),
  list(col = "Fixation_per_Nodule", title = "Fixation rate", 
       ylab = expression(paste(mu, "mol C"[2], "H"[4], " h"^-1, " g"^-1)), tag = "(c)"),
  list(col = "NMF_percent", title = "NMF", ylab = "Nodule mass fraction (%)", tag = "(d)")
)

# 2. PROCESSING FUNCTION
process_trait <- function(config, raw_data) {
  
  trait_col <- config$col
  
  # Statistics
  formula_obj <- as.formula(paste(trait_col, "~ Treatment_Full"))
  model <- aov(formula_obj, data = raw_data)
  emm <- emmeans(model, ~ Treatment_Full)
  cld_obj <- cld(emm, Letters = letters, adjust = "tukey", decreasing = T)
  
  cld_df <- as.data.frame(cld_obj) %>%
    mutate(
      letters = trimws(.group),
      y_pos = max(raw_data[[trait_col]], na.rm = TRUE) * 1.08 # Adjusted for spacing
    )
  
  # Create Plot
  p <- ggplot(raw_data, aes(x = Treatment_Full, y = .data[[trait_col]], fill = Treatment_Full)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7) +
    geom_point(position = position_jitter(width = 0.15), size = 1.2, alpha = 0.4) +
    
    # 1. Tukey Letters (Centered over boxes)
    geom_text(data = cld_df, 
              aes(x = Treatment_Full, y = y_pos, label = letters),
              size = 6, fontface = "bold", inherit.aes = FALSE, vjust = 0) +
    
    # 2. Panel Labels (Extreme top-left)
    # Using -Inf/Inf puts the text exactly at the plot borders
    annotate("text", x = -Inf, y = Inf, label = config$tag, 
             vjust = 1.5, hjust = -0.5, fontface = "bold", size = 8) +
    
    scale_fill_manual(values = site_colors) +
    # Top expansion (0.2) ensures room for the letters
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
    labs(title = config$title, y = config$ylab) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title.y = element_text(hjust = 0.5, face = "bold", size = 15),
      legend.position = "none",
      legend.text = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      panel.grid.major.x = element_blank()
    )
  
  return(list(plot = p))
}



# 3. EXECUTION & LAYOUT
results <- lapply(traits_config, process_trait, raw_data = data)
all_plots <- lapply(results, `[[`, "plot")

# Legend Extraction - Updated to force 1 column and proper alignment
legend_plot <- ggplot(data, aes(x = Treatment_Full, y = percent_N, fill = Treatment_Full)) +
  geom_boxplot() +
  scale_fill_manual(values = site_colors) +
  labs(fill = "Treatment Site") + 
  theme_minimal() + 
  theme(
    legend.position = "right",
    legend.title = element_text(hjust = 0.5, face = "bold", size = 15),
    legend.text = element_text(size = 13)
  ) +
  # This line forces the legend into a single column
  guides(fill = guide_legend(ncol = 1))

get_leg <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  return(tmp$grobs[[leg]])
}

# Combine the 4-panel plot and the legend side-by-side
# widths = c(10, 2) ensures the plots take up most of the space and the legend is on the right
final_grid <- grid.arrange(
  do.call(arrangeGrob, c(all_plots, ncol = 4)),
  get_leg(legend_plot),
  ncol = 2,
  widths = c(10, 2.5) 
)

# 4. SAVE
ggsave("plant_traits_final_formatted.png", final_grid, width = 18, height = 5.5, dpi = 300)
ggsave("plant_traits_final_formatted.svg", final_grid, width = 18, height = 5.5, dpi = 300)
