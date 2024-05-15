# Generate PCA plots

############ Data Input #########

pacman::p_load(  # Use p_load function from pacman
  dplyr,   
  ggbiplot,      # Create biplots
  magrittr,      # Pipes
  pacman,        # Load/unload packages
  rio,           # Import/export data
  tidyverse,      # So many reasons
  ggpubr,
  ggfortify,
  BioNERO
)

# P/A of C2 and C3 in each strain
Clusters_present <- c("C2","C2","both","C3","C2","C2","C2","C2","C3","none","C3","C2","C2","C2","C3","none","C3","C3","both","C2")

#label top 4 high and bottom 4 low
high_low <- c("med","med","med","med","med","med","high","high","med","med","low","high","med","high","low","med","med","low","low","med")


#################### PCA of S. meliloti logCPM expression #################

#use strain means (Ensifer)
expData <- readRDS("Rhizobia_MeanExp.RDS")

datTraits <- readRDS("TraitsInfo.RDS")

metaData <- data.frame(Strain = gsub("DZA_","",colnames(expData)),Clusters_Present = Clusters_present, ShootBiomass = datTraits$ShootBiomass, topFour = high_low)

rownames(metaData) <- metaData$Strain

expData.t <- t(expData)

# Principal components model using default method
pc <- tibble(expData.t) %>% 
#  select(-y) %>%    # Exclude variable with class labels
  prcomp(           # Most common method
    center = TRUE,  # Centers means to 0
    scale  = TRUE   # Scale to unit variance (SD = 1)
  )

# Get summary stats
pc %>% summary()

# Screeplot of eigenvalues
pc %>% plot(main = "Eigenvalues")

# Plot the projected data set on the first two principal 
expData.t <- cbind(expData.t,metaData)
rownames(expData.t) <- expData.t$Strain


#biomass + clusters #############
R_biomass_clusters <- pc %>%
  autoplot(
    data = expData.t,
    colour = 'ShootBiomass',
    label = FALSE,
    shape = FALSE
  ) +
  geom_point(aes(shape = Clusters_Present, color = ShootBiomass), size = 2) +
  scale_colour_gradientn(
    colours = rev(c(
      "#be4fa9",
      "#BE4F96",
      "#fab35e",
      "#fa995e",
      "#009680",#"#009668",
      "#00BFA1"#"#00965b"
    )),
    name = "Mean Shoot Biomass"
  ) +
  stat_ellipse(
    data = . %>% filter(topFour != "med"),
    aes(group = topFour),
    geom = "path",
    linetype = 2,
    color = "#7C7878",
    level = 0.7,
    alpha = 0.4
  ) +
  theme_minimal(base_size = 8) + xlab("PC1 (16.51% explained var.)") + ylab("PC2 (11.63% explained var.)") +
  geom_text(aes(label = Strain, colour = ShootBiomass),
            vjust = -.6,
            size = 3) +
  guides(shape = guide_legend(title = "Clusters Present")) +
  labs(title = "S. meliloti") +  
  theme(
    plot.title = element_text(face = "italic", size = 8),
    legend.key.size = unit(0.3, 'cm'),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 8)
  )

ggsave("figures/Rhizobia_PCA.tiff",R_biomass_clusters,dpi = 600,units = "cm", width = 11,height = 9, bg = "white")



# Clear data
rm(metaData,pc, expData, expData.t)  # Removes extra objects from the environment


#################### PCA of Medicago logCPM expression ################# 

#use strain means (Medicago)
expData <- readRDS("Medicago_MeanExp.RDS")

metaData <- data.frame(Strain = gsub("DZA_","",colnames(expData)),Clusters_Present = Clusters_present, ShootBiomass = datTraits$ShootBiomass, topFour = high_low)

rownames(metaData) <- metaData$Strain

#pick top 5000 highly variable genes
#expData.t <- t(BioNERO::filter_by_variance(expData, n = 5000))
expData.t <- readRDS("Medicago_5000Genes_MeanExp.RDS")

##########################################

# Principal components model using default method
pc <- tibble(expData.t) %>% 
  #  select(-y) %>%    # Exclude variable with class labels
  prcomp(           # Most common method
    center = TRUE,  # Centers means to 0
    scale  = TRUE   # Scale to unit variance (SD = 1)
  )

# Get summary stats
pc %>% summary()

# Screeplot of eigenvalues
pc %>% plot(main = "Eigenvalues")

# Plot the projected data set on the first two principal 
expData.t <- cbind(expData.t,metaData)
rownames(expData.t) <- expData.t$Strain

#biomass + clusters #############
M_biomass_clusters <- pc %>%
  autoplot(
    data = expData.t,
    colour = 'ShootBiomass',
    label = FALSE,
    shape = FALSE
  ) +
  geom_point(aes(shape = Clusters_Present, color = ShootBiomass), size = 2) +
  scale_colour_gradientn(
    colours = rev(c(
      "#be4fa9",
      "#BE4F96",
      "#fab35e",
      "#fa995e",
      "#009680",#"#009668",
      "#00BFA1"#"#00965b"
    )),
    name = "Mean Shoot Biomass"
  ) +
  stat_ellipse(
    data = . %>% filter(topFour != "med"),
    aes(group = topFour),
    geom = "path",
    linetype = 2,
    color = "#7C7878",
    level = 0.7,
    alpha = 0.4
  ) +
  theme_minimal(base_size = 8) + xlab("PC1 (35.94% explained var.)") + ylab("PC2 (12.14% explained var.)") +
  geom_text(aes(label = Strain, colour = ShootBiomass),
            vjust = -.6,
            size = 3) +
  guides(shape = guide_legend(title = "Clusters Present")) +
  labs(title = "M. truncatula") +  
  theme(
    plot.title = element_text(face = "italic", size = 8),
    legend.key.size = unit(0.3, 'cm'),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    axis.title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 8)
  )

ggsave("figures/Medicago_top5000Genes_PCA.tiff",M_biomass_clusters,dpi = 600,units = "cm", width = 11,height = 9, bg = "white")



# MULTI-PANEL FIGURE
arng <- ggarrange(R_biomass_clusters, M_biomass_clusters,
          labels = c("A", "B"),
          ncol = 2, nrow = 1, align = "h",font.label = list(size = 8, color = "black", family = "Arial"))

ggsave("figures/multi_panel_PCAs.tiff",arng,dpi = 600,units = "cm", width = 20,height = 9, bg = "white")

