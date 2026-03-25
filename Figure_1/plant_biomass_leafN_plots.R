#!/usr/bin/env Rscript

# Load libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(emmeans)
library(multcomp)
library(multcompView)

setwd("C:/Temp")

# 1. DATA PREPARATION
data <- read_excel("G:/My Drive/Moshe/Efrat_Guy_Project/Guy_Samples/Guy_data.xlsx", sheet = 1)

treatment_map <- c(
  "CarK" = "CAR-C",
  "CarR" = "CAR-G",
  "Ces"  = "CES",
  "Dgn"  = "DGN",
  "Hok"  = "HUK",
  "Mtz"  = "MTZ",
  "Prl"  = "CTRL (Prl)",
  "RH"   = "RH",
  "Ya"   = "YTR"
)

# 10 visually distinct colors (Set2 + extended)
site_colors <- c(
  "CAR-C"      = "#66C2A5",  # original
  "CAR-G"      = "#FC8D62",  # original
  "CES"        = "#8DA0CB",  # original
  "CTRL (Prl)" = "#80B1D3",  # new
  "DGN"        = "#E78AC3",  # original
  "HUK"        = "#A6D854",  # original
  "MTZ"        = "#FFD92F",  # original
  "RH"         = "#E5C494",  # original
  "YTR"        = "#B3B3B3"   # original (YATIR)
)

treatment_levels <- names(site_colors)

data <- data %>%
  filter(treat %in% names(treatment_map)) %>%
  rename(pctN = `%N`) %>%
  mutate(
    Treatment_Full = factor(treatment_map[treat], levels = treatment_levels),
    biomass = as.numeric(biomass),
    pctN    = as.numeric(pctN)
  )

# 2. PROCESSING FUNCTION (same pattern as original)
process_trait <- function(col, title, ylab, tag, raw_data) {

  df <- raw_data %>% filter(!is.na(.data[[col]]))

  formula_obj <- as.formula(paste0("`", col, "` ~ Treatment_Full"))
  model <- aov(formula_obj, data = df)
  emm   <- emmeans(model, ~ Treatment_Full)
  cld_obj <- cld(emm, Letters = letters, adjust = "tukey", decreasing = TRUE)

  y_max <- max(df[[col]], na.rm = TRUE) * 1.08
  cld_df <- as.data.frame(cld_obj) %>%
    mutate(
      letters = trimws(.group),
      y_pos   = y_max
    )

  p <- ggplot(df, aes(x = Treatment_Full, y = .data[[col]], fill = Treatment_Full)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7) +
    geom_point(position = position_jitter(width = 0.15), size = 1.2, alpha = 0.4) +
    geom_text(data = cld_df,
              aes(x = Treatment_Full, y = y_pos, label = letters),
              size = 6, fontface = "bold", inherit.aes = FALSE, vjust = 0) +
    annotate("text", x = -Inf, y = Inf, label = tag,
             vjust = 1.5, hjust = -0.5, fontface = "bold", size = 8) +
    scale_fill_manual(values = site_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    labs(title = title, y = ylab) +
    theme_bw() +
    theme(
      axis.title.x  = element_blank(),
      axis.text.x   = element_blank(),
      axis.ticks.x  = element_blank(),
      axis.text.y   = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title.y  = element_text(hjust = 0.5, face = "bold", size = 15),
      legend.position  = "none",
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
      panel.grid.major.x = element_blank()
    )

  return(p)
}

# 3. GENERATE PLOTS
p1 <- process_trait("pctN",    "Leaf N",        "Leaf N (%)",        "(a)", data)
p2 <- process_trait("biomass", "Total Biomass", "Total Biomass (g)", "(b)", data)

# 4. EXTRACT SHARED LEGEND
legend_plot <- ggplot(data, aes(x = Treatment_Full, y = pctN, fill = Treatment_Full)) +
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

# 5. COMBINE & SAVE
final_grid <- arrangeGrob(
  p1, p2, get_leg(legend_plot),
  ncol = 3,
  widths = c(5, 5, 2)
)

out_base <- "G:/My Drive/Moshe/Efrat_Guy_Project/Final_Figures/plant_biomass_leafN"

ggsave(paste0(out_base, ".png"), final_grid, width = 14, height = 5.5, dpi = 300)
ggsave(paste0(out_base, ".pdf"), final_grid, width = 14, height = 5.5)
ggsave(paste0(out_base, ".svg"), final_grid, width = 14, height = 5.5)

cat("Done! Saved to", out_base, "\n")
