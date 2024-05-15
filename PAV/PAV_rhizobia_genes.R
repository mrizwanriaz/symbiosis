################
# Presence-Absence Variation analysis
#################


library(limma)
library(gplots)
library(Hmisc)
library(tidyverse)
library(WGCNA)
library(grid)
library(ComplexHeatmap)
library(ggplot2)
library(patchwork)
library(psych)
library(ggpubr)

#set working directory


######## STEP 1 : Identify low expression genes #########

expData <- readRDS("expData_rhizobia.RDS")

cpm.DZA <- readRDS("cpm_rhizobia.RDS")
dim(cpm.DZA) #6193 80
 
min_samp <- 3
 
lowExprGenes <- list()

for (i in seq(1, 80, by=4)){

  #DZA
  #median library size
  M <- median(expData$samples[c(i:(i+3)),2]) * 1e-6
  cpm_cutoff.DZA=4/M

  # extract genes with expression below the cutoff in atleast three out of four replicates
  filt.DZA <- rowSums(cpm.DZA[,c(i:(i+3))] <= cpm_cutoff.DZA) >= min_samp

  strainID.DZA <- (str_split(colnames(cpm.DZA)[i],"_"))[[1]][2]

  write.csv(expData[filt.DZA,c(i:(i+3))],paste("low_expression_genelists/",strainID.DZA,".csv",sep = ""))

  lowExprGenes[strainID.DZA] <- list(rownames(expData[filt.DZA,]))


}

View(lowExprGenes)
# 
# 
#UpSet plot 
library(UpSetR)
# 
# 
# #order the sets based on shoot biomass
strains <- c("204","206","278","283","141","199","216","144B","176A","182","121","89","282","230","335","724B","25B","748B","31","746B")
 
jpeg(paste("figs/lowExprGenes_ordered.jpeg",sep=""),width=15,height=9,units="in",res=600,quality=100)

upset(fromList(lowExprGenes), nsets = 20, keep.order = TRUE, sets = rev(strains),
      nintersects = 40, number.angles = 360, order.by = "freq",
      decreasing = T,
      mb.ratio = c(0.7, 0.3),
      text.scale = 1.5,
      point.size = 1,
      line.size = 0.6,
      group.by = "degree",
      sets.bar.color = "#a5aa99",
      main.bar.color ="black",
      matrix.color = "#708090",
      cutoff = 5,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Number of Low Expression Genes",

)
dev.off()

#lowExprGenes <- readRDS("absent_genes/Ensifer/DZA/gene_lists/confirmation/lowExprGenes_rhizobia.RDS")


######### STEP 2 : filtering of low expression genes based on their presence-absence pattern (present in at least three and maximum 17 strains) #######################

GeneList <- Reduce(union,lowExprGenes) # union of all genes
length(GeneList) #1102 genes
# 
#create to dataframe to hold P/A data
#Genes_PA=data.frame(matrix(ncol = 20, nrow = length(GeneList)))
# colnames(Genes_PA) <- strains
# rownames(Genes_PA) <- GeneList
# 
# # #populate data frame with gene presence-absence, 1=absent and 0=present
# Genes_PA[,1] <- ifelse(GeneList %in% as.list(lowExprGenes$`204`),1,0)
# Genes_PA[,2] <- ifelse(GeneList %in% as.list(lowExprGenes$`206`),1,0)
# Genes_PA[,3] <- ifelse(GeneList %in% as.list(lowExprGenes$`278`),1,0)
# Genes_PA[,4] <- ifelse(GeneList %in% as.list(lowExprGenes$`283`),1,0)
# Genes_PA[,5] <- ifelse(GeneList %in% as.list(lowExprGenes$`141`),1,0)
# Genes_PA[,6] <- ifelse(GeneList %in% as.list(lowExprGenes$`199`),1,0)
# Genes_PA[,7] <- ifelse(GeneList %in% as.list(lowExprGenes$`216`),1,0)
# Genes_PA[,8] <- ifelse(GeneList %in% as.list(lowExprGenes$`144B`),1,0)
# Genes_PA[,9] <- ifelse(GeneList %in% as.list(lowExprGenes$`176A`),1,0)
# Genes_PA[,10] <- ifelse(GeneList %in% as.list(lowExprGenes$`182`),1,0)
# Genes_PA[,11] <- ifelse(GeneList %in% as.list(lowExprGenes$`89`),1,0)
# Genes_PA[,12] <- ifelse(GeneList %in% as.list(lowExprGenes$`282`),1,0)
# Genes_PA[,13] <- ifelse(GeneList %in% as.list(lowExprGenes$`230`),1,0)
# Genes_PA[,14] <- ifelse(GeneList %in% as.list(lowExprGenes$`121`),1,0)
# Genes_PA[,15] <- ifelse(GeneList %in% as.list(lowExprGenes$`335`),1,0)
# Genes_PA[,16] <- ifelse(GeneList %in% as.list(lowExprGenes$`724B`),1,0)
# Genes_PA[,17] <- ifelse(GeneList %in% as.list(lowExprGenes$`25B`),1,0)
# Genes_PA[,18] <- ifelse(GeneList %in% as.list(lowExprGenes$`748B`),1,0)
# Genes_PA[,19] <- ifelse(GeneList %in% as.list(lowExprGenes$`31`),1,0)
# Genes_PA[,20] <- ifelse(GeneList %in% as.list(lowExprGenes$`746B`),1,0)
# 
Genes_PA <- readRDS("Genes_PA_1102.RDS")
# 
# #keep genes that are absent (1) in mininum of three and maximum 17 strains
Genes_PA <- Genes_PA[rowSums(Genes_PA >= 1) %in% c(3:17),]
# 
#reduced to 574 gene
dim(Genes_PA)

Genes_PA.filt <- readRDS("Genes_PA_filt_574.RDS")
View(Genes_PA.filt) #574 genes
colnames(Genes_PA.filt) <- paste("DZA_",colnames(Genes_PA.filt),sep = "")

######### STEP 3 : compute the gene PA correlation with plant shoot biomass - cutoff (|r| > 0.5) ################

# read Trait Data
traitData <- readRDS("traitData.RDS")

# correlate P/A variation of 574 genes with plant shoot biomass
#corr_df <- corr.test(t(Genes_PA.filt),traitData[,2], method = "pearson", adjust="fdr")
corr_df <- readRDS("corr_574genes_shootbiomass.RDS")
cor_mat <- data.frame(cor_shootBiomass = corr_df$r,p_value = corr_df$p,fdr=corr_df$p.adj)

highCorr_genes <- rownames(cor_mat[abs(cor_mat$cor_shootBiomass) >= 0.5,]) #82 genes
highCorr_genes

###### STEP 4 : BLAST "highCorr_genes" against long-read sequences and confirmed that 54 have high-confidence presence-absence variation across strains ----- ####


###### STEP 5 : select genes presnt in clusters ##########

#HC_candidates <- readRDS("HighConfidence_59Absentgenes.RDS")

# corrected PA based on alignment for "WP_010967187","WP_013845919"

# re-compute correlation between gene presence-absence and shoot biomass
#cor_df <- data.frame(cor(t(HC_candidates[,c(2:21)]),traitData))
#View(cor_df)

#HC_highCorr_genes <- rownames(cor_df[abs(cor_df$ShootBiomass) >= 0.5,])
#HC_highCorr_genes #resulted in 57 genes (two genes with corrected PA were removed because of low correlation)


#HC_highCorr_57Absentgenes <- data.frame(t(HC_candidates[HC_highCorr_genes,c(2:21)]))
#dim(HC_highCorr_57Absentgenes)
#HC_highCorr_57Absentgenes <- readRDS("HC_highCorr_57Absentgenes.RDS")

#HC_highCor_ClusterGenes <- HC_highCorr_57Absentgenes %>% select(-c("WP_010967204.1", "SM_RS20255", "SM_RS30605", "WP_010975890.1", "WP_015242546.1", "WP_010968014.1")) #include only genes that are part of ASSIGNED cluster of psymA, 6 clusters (51 genes)

HC_highCor_ClusterGenes<-readRDS("HC_highCor_51ClusterGenes.RDS")
View(HC_highCor_ClusterGenes)


rownames(HC_highCor_ClusterGenes) <- gsub("DZA_","",rownames(HC_highCor_ClusterGenes))

####### STEP 6: compute strain-strain correlation, peform t-tests for 20 and 171 strain dataset and generate figures for manuscript ##############

#cor_strains <- rcorr(t(HC_highCor_ClusterGenes))
#View(cor_strains)
cor_strains <- readRDS("strain_strain_correlation.RDS")
cor_strains_mat <- cor_strains$r

# # sometimes, you may have NA in the matrix, and clustering does not play well with it
# # a simple hack is to turn the NA to 0
cor_strains_mat[is.na(cor_strains_mat)]<- 0
# #cor_strains_mat[cor_strains_mat == 1]<- 0
# 
cor_strains_p<- cor_strains$P
cor_strains_p[is.na(cor_strains_p)]<- 1
cor_strains_p[is.nan(cor_strains_p)]<- 1

# 
# #Only label the correlation coefficients with p-values that are <=0.05; add * for p value <=0.05 and ** for p value <=0.01
cell_fun = function(j, i, x, y, w, h, fill){
  if(as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
    grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
  }

  if (cor_strains_p[i, j]  < 0.01 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
    grid.text(paste0(sprintf("%.2f", cor_strains_mat[i, j]),"**"), x, y, gp = gpar(fontsize = 9))
  } else if (cor_strains_p[i, j]  <= 0.05 & as.numeric(x) <= 1 - as.numeric(y) + 1e-6){
    grid.text(paste0(sprintf("%.2f", cor_strains_mat[i, j]),"*"), x, y, gp = gpar(fontsize = 9))
  }
}
 

col_fun<- circlize::colorRamp2(c(-1, 0, 1),c("#33BBee","white","#EE7733"))

hp<- ComplexHeatmap::Heatmap(cor_strains_mat,
                             rect_gp = gpar(type = "none"),
                             column_dend_side = "bottom",
                             #column_title = "Gene PAV correlation",
                             name = "correlation", col = col_fun,
                             cell_fun = cell_fun,
                             cluster_rows = F, cluster_columns = F,
                             row_names_side = "left",
                             row_names_gp = gpar(fontsize = 9.5, fontface = "bold"),
                             column_names_gp = gpar(fontsize = 9.5,fontface = "bold"))


lgd_list = list(
  Legend( labels = c("<0.01", "<0.05"), title = "pvalue",
          graphics = list(
            function(x, y, w, h) grid.text("**", x = x, y = y,
                                           gp = gpar(fill = "black")),
            function(x, y, w, h) grid.text("*", x = x, y = y,
                                           gp = gpar(fill = "black")))
  ))

png("figs/20Strain_Cor.png", res = 600, units = "in", width = 10, height = 4.5)
draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(0.1, "cm") )

dev.off()


###### PAV vs shoot biomass boxplots  ######
rownames(traitData) <- gsub("DZA_","",rownames(traitData))



### CLUSTER 2 (previously cluster 1)
cluster_2 <- c("WP_010967274.1", "WP_010967276.1","WP_010967277.1","SM_RS26160")
cluster_2_df <- data.frame(HC_highCor_ClusterGenes[,cluster_2])
cluster_2_df$rowSum <- rowSums(cluster_2_df)
#keep strains with all or none condition for cluster_1
cluster_2_df$Cluster2 <- ifelse(cluster_2_df$rowSum == 4,"Absent","Present")

cluster_2_df$shootBiomass <- traitData[rownames(cluster_2_df),2]
cluster_2_df$Cluster2 <- factor(cluster_2_df$Cluster2 , levels=c("Present", "Absent"))

t.test(shootBiomass ~ Cluster2, data = cluster_2_df[,c(6,7)], var.equal = FALSE) #p=0.005829

boxplot_colrs <- c("grey","white","#BE4F96","#009680","black")
fill_present = boxplot_colrs[1]
fill_absent = boxplot_colrs[2]
outline_high = boxplot_colrs[5]
outline_low = boxplot_colrs[5]

p<-ggboxplot(cluster_2_df[, c(6, 7)], x = "Cluster2", y = "shootBiomass",
          color = "Cluster2", palette = c(outline_high, outline_low),
          fill = "Cluster2",
          add = "jitter", shape = "Cluster2") +
  xlab("C2") + 
  ylab("Average Shoot Biomass") +
  stat_compare_means(comparisons = list(c("Absent", "Present")), method = "t.test", label.y = 0.49) +
  theme(legend.title = element_blank(), legend.position = "right")+
  scale_fill_manual(values = c(fill_present,fill_absent))+
  scale_shape_manual(values = c(16, 16)) # Using shape 16 for jittered dots

ggsave("figs/cluster2_PAV_20strains.tiff",p, dpi = 600, device = "tiff", units = "in", width = 3.5, height = 3.6)


### CLUSTER 3 (previously cluster 2)
cluster_3 <- c("WP_010967494.1", "WP_010967495.1", "WP_010967496.1", "WP_010967497.1", "WP_010967498.1", "WP_010967499.1", "WP_010967500.1", "WP_010967501.1", "WP_010967502.1")

cluster_3_df <- data.frame(HC_highCor_ClusterGenes[,cluster_3])
cluster_3_df$rowSum <- rowSums(cluster_3_df)
cluster_3_df$Cluster3 <- ifelse(cluster_3_df$rowSum == 9,"Absent","Present")

cluster_3_df$shootBiomass <- traitData[rownames(cluster_3_df),2]
cluster_3_df$Cluster3 <- factor(cluster_3_df$Cluster3 , levels=c("Present", "Absent"))

t.test(shootBiomass ~ Cluster3, data = cluster_3_df[,c(11,12)], var.equal = FALSE) #p=0

p<-ggboxplot(cluster_3_df[, c(11, 12)], x = "Cluster3", y = "shootBiomass",
             color = "Cluster3", palette = c(outline_low,outline_high),
             fill = "Cluster3",
             add = "jitter", shape = "Cluster3") +
  xlab("C3") + 
  ylab("Average Shoot Biomass") +
  stat_compare_means(comparisons = list(c("Absent", "Present")), method = "t.test", label.y = 0.49) +
  theme(legend.title = element_blank(), legend.position = "right")+
  scale_fill_manual(values = c(fill_present,fill_absent))+
  scale_shape_manual(values = c(16, 16)) # Using shape 16 for jittered dots

ggsave("figs/cluster3_PAV_20strains.tiff",p, dpi = 600, device = "tiff", units = "in", width = 3.5, height = 3.6)


##### CORRELATION OF GENE PRESENCE_ABSENCE and SHOOT BIOMASS in 171 strains ##############
HC_highCorr_genes_171strains <- readRDS("genePA_57genes_171strains.RDS")

strains191_traitData <- readRDS("strains191_traitData.RDS")

# t-test of biomass based on presence absence of cluster_2 and cluster_3

### CLUSTER 2 (prev cluster 1)
cluster_2 <- c("WP_010967274.1", "WP_010967276.1","WP_010967277.1","SM_RS26160_WP_234826226.1")
cluster_2_df <- data.frame(t(HC_highCorr_genes_171strains[cluster_2,]))
cluster_2_df$rowSum <- rowSums(cluster_2_df)
#keep strains with all or none condition for cluster_2
cluster_2_df <- filter(cluster_2_df,(cluster_2_df$rowSum == 4 | cluster_2_df$rowSum == 0))
cluster_2_df$Cluster2 <- ifelse(cluster_2_df$rowSum == 4,"Absent","Present")

cluster_2_df$shootBiomass <- strains191_traitData[rownames(cluster_2_df),2]
#remove strain with NA biomass
cluster_2_df <- cluster_2_df[!(row.names(cluster_2_df) %in% c("DZA_26")),]
cluster_2_df$Cluster2 <- factor(cluster_2_df$Cluster2 , levels=c("Present", "Absent"))

t.test(shootBiomass ~ Cluster2, data = cluster_2_df[,c(6,7)], var.equal = FALSE) #p=0

p<-ggboxplot(cluster_2_df[, c(6, 7)], x = "Cluster2", y = "shootBiomass",
             color = "Cluster2", palette = c(outline_high, outline_low),
             fill = "Cluster2",
             add = "jitter", shape = "Cluster2") +
  xlab("C2") + 
  ylab("Average Shoot Biomass") +
  stat_compare_means(comparisons = list(c("Absent", "Present")), method = "t.test") +
  theme(legend.title = element_blank(), legend.position = "right")+
  scale_fill_manual(values = c(fill_present,fill_absent))+
  scale_shape_manual(values = c(16, 16)) # Using shape 16 for jittered dots

ggsave("figs/cluster2_PAV_160strains.tiff",p, dpi = 600, device = "tiff", units = "in", width = 3.5, height = 3.6)



### CLUSTER 3 (prev cluster 2)
cluster_3 <- c("WP_010967494.1", "WP_010967495.1", "WP_010967496.1", "WP_010967497.1", "WP_010967498.1", "WP_010967499.1", "WP_010967500.1", "WP_010967501.1", "WP_010967502.1", "WP_028011575.1")

cluster_3_df <- data.frame(t(HC_highCorr_genes_171strains[cluster_3,]))
cluster_3_df$rowSum <- rowSums(cluster_3_df)
#keep strains with all or none condition for cluster_2
cluster_3_df <- filter(cluster_3_df,(cluster_3_df$rowSum == 10 | cluster_3_df$rowSum == 0))
cluster_3_df$Cluster3 <- ifelse(cluster_3_df$rowSum == 10,"Absent","Present")

cluster_3_df$shootBiomass <- strains191_traitData[rownames(cluster_3_df),2]
#remove strain with NA biomass
cluster_3_df <- cluster_3_df[!(row.names(cluster_3_df) %in% c("DZA_26")),]
cluster_3_df$Cluster3 <- factor(cluster_3_df$Cluster3 , levels=c("Present", "Absent"))

t.test(shootBiomass ~ Cluster3, data = cluster_3_df[,c(12,13)], var.equal = FALSE) #p=0

p<-ggboxplot(cluster_3_df[, c(12, 13)], x = "Cluster3", y = "shootBiomass",
             color = "Cluster3", palette = c(outline_low, outline_high),
             fill = "Cluster3",
             add = "jitter", shape = "Cluster3") +
  xlab("C3") + 
  ylab("Average Shoot Biomass") +
  stat_compare_means(comparisons = list(c("Absent", "Present")), method = "t.test") +
  theme(legend.title = element_blank(), legend.position = "right")+
  scale_fill_manual(values = c(fill_present,fill_absent))+
  scale_shape_manual(values = c(16, 16)) # Using shape 16 for jittered dots

ggsave("figs/cluster3_PAV_152strains.tiff",p, dpi = 600, device = "tiff", units = "in", width = 3.5, height = 3.6)




