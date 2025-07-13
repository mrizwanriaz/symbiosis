# Load required libraries
library(tidyverse)   # Includes dplyr, ggplot2, etc.
library(ggpubr)      # Publication-ready plots
library(ggplot2)     # Grammar of graphics (redundant if tidyverse is loaded)
library(agricolae)   # For Duncan's multiple range test

# Load the dataset
biomass <- readRDS("biomass_dataframe.RDS")

# =======================
# Basic Descriptive Stats
# =======================

# Sample variance of ShootBiomass
var_shoot <- var(biomass$ShootBiomass)  # Estimate of spread in the data

# Coefficient of Variation (CV) in percentage
cv_shoot <- sd(biomass$ShootBiomass) / mean(biomass$ShootBiomass) * 100  # Measure of relative variability

# ================================
# Grouped Summary: Mean by Groups
# ================================

mean_shoot_biomass <- biomass %>%
  group_by(Groups) %>%
  summarise(
    MeanShootBiomass = mean(ShootBiomass, na.rm = TRUE),
    .groups = 'drop'  # Ungroup after summarising
  )

# ======================
# Boxplot Visualization
# ======================

p1 <- ggboxplot(
  biomass,
  x = "Groups",
  y = "ShootBiomass",
  fill = "Groups",                  # Fill boxes by group
  bxp.errorbar = TRUE,              # Add error bars to boxplot
  bxp.errorbar.width = 0.4,
  add = "jitter",                   # Add individual points with jitter
  outlier.shape = 16,
  size = 0.4,
  add.params = list(size = 0.4)     # Size of jittered points
) + 
  theme(
    legend.position = "none",
    axis.title = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(
      size = 8,
      face = "italic",
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 8)
  ) +
  scale_fill_manual(values = rev(c(
    "grey", "#419460", "#419462", "#419467", "#419469", "#5F966B", "#65986A",
    "#65986A", "#77996B", "#A59D6A", "#B59D6A", "#BA9D69", "#BE9E69",
    "#EC9E69", "#EEAD6B", "#DE9B78", "#CE8084", "#B15696", "#B15697",
    "#B15698", "#B155A5"
  ))) +
  labs(
    x = "Strains",
    y = "Shoot Biomass"
  ) + 
  labs_pubr()  # Apply publication-style theme

# Display boxplot
print(p1)

# =====================
# One-Way ANOVA Testing
# =====================

anova <- aov(ShootBiomass ~ Strain, data = biomass)

# View ANOVA summary table
anova_result <- summary(anova)
print(anova_result)

# ===========================
# Post-hoc: Duncan's Test
# ===========================

# Duncan's Multiple Range Test (for unequal variances)
duncan_result <- duncan.test(anova, "Strain", alpha = 0.05)
print(duncan_result)

# ================================
# Annotation with Group Letters
# ================================

# Create annotation table: mean + 3rd quartile (for label positioning)
duncan_annot <- biomass %>%
  group_by(Strain) %>%
  summarise(
    mean = mean(ShootBiomass),
    quant = quantile(ShootBiomass, probs = 0.75),
    .groups = "drop"
  ) %>%
  arrange(desc(mean))

# Add Compact Letter Display (CLD) from Duncan's test
duncan_annot$cld <- duncan_result$groups$groups

# Add annotation to the existing boxplot
p2 <- p1 + geom_text(
  data = duncan_annot,
  aes(x = Strain, y = quant, label = cld),
  size = 3,
  vjust = -4.5,
  nudge_x = 0.15
)

# Display the annotated plot
print(p2)
