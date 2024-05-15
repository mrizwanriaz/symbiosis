library(tidyverse)
library(ggpubr)
library(ggplot2)

biomass <- read.csv("DZA_trait_info2.csv")

# sample variance
var(biomass$ShootBiomass) #0.018
# coefficient of variation
sd(biomass$ShootBiomass) / mean(biomass$ShootBiomass) * 100 #44.24

biomass$Groups <- gsub("sm_","",biomass$Groups) 
biomass$Strain <- factor(biomass$Groups)


# Sort the Groups factor by the mean ShootBiomass
biomass$Groups <-
  factor(
    biomass$Groups,
    levels = c(
      "204",
      "206",
      "278",
      "283",
      "141",
      "199",
      "216",
      "144B",
      "176A",
      "182",
      "121",
      "89",
      "282",
      "230",
      "335",
      "724B",
      "25B",
      "748B",
      "31",
      "746B"
    )
  )


p1 <-
  ggboxplot(
    biomass,
    x = "Groups",
    y = "ShootBiomass",
    fill = "Groups",
    bxp.errorbar = TRUE,
    bxp.errorbar.width = 0.4,
    add = c("jitter"),
    outlier.shape = 16,
    size = 0.4,
    add.params = list(size = 0.4),#jitter size
  ) + theme(
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
  scale_fill_manual(values =
                      rev(c("#419460","#419462","#419467","#419469","#5F966B","#65986A",
                        "#65986A","#77996B","#A59D6A","#B59D6A","#BA9D69","#BE9E69",
                        "#EC9E69","#EEAD6B","#DE9B78","#CE8084","#B15696","#B15697",
                        "#B15698","#B155A5"
                      ))) + #rep("#CCCCCC", 20))) + # Set the fill color to neutral grey
  labs(x = "Strains", y = "Shoot Biomass") + 
  labs_pubr() 
p1


# Save the plot as a PNG with 600 dpi resolution
ggsave(
  "figures/shootbiomass_boxplot.tiff",
  plot = p1,
  dpi = 400,
  width = 18,
  height = 10,
  units = "cm"
)


# Perform ANOVA
anova <- aov(ShootBiomass ~ Strain, data = biomass)

# Check ANOVA results
anova_result <- summary(anova)

# Print ANOVA table
print(anova_result)

capture.output(anova_result, file = "shootbiomass_anova.txt")


############### Duncan's test #############
library(agricolae)

# Perform Duncan's test (it assumes variances are not equal among groups)
duncan_result <- duncan.test(anova, "Strain", alpha = 0.05)

# View the Duncan's test results
print(duncan_result)

write.csv(duncan_result$groups,"Groups_DuncanTest.csv")

# table with factors and 3rd quantile
duncan_annot <- group_by(biomass, Strain) %>%
  summarise(mean=mean(ShootBiomass), quant = quantile(ShootBiomass, probs = 0.75)) %>%
  arrange(desc(mean))

# extracting the compact letter display and add it to the figure
duncan_annot$cld <- duncan_result$groups$groups

p2 <- p1 + geom_text(data = duncan_annot, aes(x = Strain, y = quant, label = cld), size = 3, vjust=-4.5, nudge_x = 0.15)

# Save the plot as a TIFF
ggsave(
  "figures/shootbiomass_boxplot_annotDuncan.tiff",
  plot = p2,
  dpi = 400,
  width = 18,
  height = 10,
  units = "cm"
)

