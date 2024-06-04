#R code to perform Weighted Gene Co-expression Network Analysis using Rhizobia Genes

library(WGCNA)

library(gplots)


# setwd() #set working directory

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

#load averageLogCPM values (calculated using aveLogCPM function)
expData.Avelcpm <- readRDS("expData_AverageLogCPM.RDS")
#5993 genes after filtering for low count genes

datExpr = as.data.frame(t(expData.Avelcpm[, ]))

names(datExpr) = rownames(expData.Avelcpm)

rownames(datExpr) = names(expData.Avelcpm)


gsg = goodSamplesGenes(datExpr, verbose = 3)

gsg$allOK

#cluster the samples to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr), method = "average")

jpeg(
        "Figures/sampleClustering_Rhizobia.jpeg",
        width = 20,
        height = 10,
        units = "in",
        res = 300,
        quality = 100
)
par(cex = 1)

par(mar = c(0, 4, 2, 0))
plot(
        sampleTree,
        main = "Sample clustering to detect outliers",
        sub = "",
        xlab = "",
        cex.lab = 1.5,
        cex.axis = 1.5,
        cex.main = 2
)
dev.off()

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# # Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
# read output of softConnectivity
sft <- readRDS("sft_rhizobia.RDS")

# Plot the results:
jpeg(
        "Figures//softThreshold_Rhizobia.jpeg",
        width = 15,
        height = 10,
        units = "in",
        res = 300,
        quality = 100
)
par(mfrow = c(1, 2))

cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(
        sft$fitIndices[, 1],
        -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        xlab = "Soft Threshold (power)",
        ylab = "Scale Free Topology Model Fit,signed R^2",
        type = "n",
        main = paste("Scale independence")
)

text(
        sft$fitIndices[, 1],
        -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
        labels = powers,
        cex = cex1,
        col = "red"
)

# this line corresponds to using an R^2 cut-off of h
abline(h = 0.865, col = "blue")
# Mean connectivity as a function of the soft-thresholding power
plot(
        sft$fitIndices[, 1],
        sft$fitIndices[, 5],
        xlab = "Soft Threshold (power)",
        ylab = "Mean Connectivity",
        type = "n",
        main = paste("Mean connectivity")
)
text(
        sft$fitIndices[, 1],
        sft$fitIndices[, 5],
        labels = powers,
        cex = cex1,
        col = "red"
)
abline(h = 54, col = "blue")
dev.off()


# we chose B=4
softPower = 4

# Network Construction - BLockwise
# net <- blockwiseModules(datExpr,
#                         power =softPower,
#                         maxBlockSize = 6000, #should be larger than num. of genes
#                         corType = "bicor", #lecture for caveats
#                         maxPOutliers = 0.04,
#                         networkType = "signed hybrid",
#                         loadTOM = FALSE, #nothing to load the first time
#                         saveTOMs = FALSE, #save the TOMS to re-load at a later time
#                         deepSplit = 2, #intermediate level
#                         minModuleSize = 10, #a module has to have this many genes
#                         mergeCutHeight = 0.2, #level for merging; suggest 0.2 to start
#                         numericLabels = TRUE,
#                         verbose = 3 )

#read network object
net <- readRDS("network_Rhizobia.RDS")

table(net$colors)

table(net$unmergedColors)

mergedColors <- labels2colors(net$colors)

jpeg(
        "Figures//GeneDendrogram_ModuleColors_Rhizobia.jpeg",
        width = 20,
        height = 10,
        units = "in",
        res = 300,
        quality = 100
)
plotDendroAndColors(
        net$dendrograms[[1]],
        mergedColors[net$blockGenes[[1]]],
        "Module Colors",
        main = "Gene dendrogram and module colors",
        dendroLabels = FALSE,
        hang = 0.03,
        addGuide = FALSE,
        guideHang = 0.05
)
dev.off()


# Calculate eigengenes
#MEList = moduleEigengenes(datExpr, colors = net$colors)
MEList <- readRDS("MEList_rhizobia.RDS")
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1 - cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
jpeg(
        "Figures//ClusteringModules_Rhizobia.jpeg",
        width = 20,
        height = 10,
        units = "in",
        res = 300,
        quality = 100
)
plot(METree,
     main = "Clustering of module eigengenes (DZA - Ensifer)",
     xlab = "",
     sub = "")
dev.off()

############ Relate modules to traits ######

#read trait data
allTraits <- readRDS("traits_info.RDS")

# Form a data frame analogous to expression data that will hold the traits.
allSamples = rownames(datExpr)

traitRows = match(allSamples, allTraits$strainID)

datTraits = allTraits[traitRows, -1]

rownames(datTraits) = allTraits[traitRows, 1]

#re-order MEs such that similar ones are next to each other
MEs = orderMEs(MEs, greyLast = TRUE, greyName = "ME0")

#compute correlation
#moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitCor <- readRDS("moduleTraitCor_Rhizobia.RDS")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

colnames(moduleTraitPvalue) <-
        paste("p-value.", colnames(moduleTraitPvalue), sep = "")
#moduleDF <- cbind(moduleTraitCor,moduleTraitPvalue)
#write.csv(moduleDF, "module-trait-correlation_Rhizobia.csv")


#plot
library(rcartocolor)

my_color <- colorpanel(150, "#cceeff", "white", "#ffccee")

png(
        "Figures//ModuleTraitRelationship_Rhizobia.png",
        width = 8,
        height = 9,
        units = "in",
        res = 600
)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor[c(1:40), ], 2),
                    " (",
                    signif(moduleTraitPvalue[c(1:40), ], 1),
                    ")",
                    sep = "")

dim(textMatrix) = dim(moduleTraitCor[c(1:40), ])
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
p <- labeledHeatmap(
        Matrix = moduleTraitCor[c(1:40), ],
        xLabels = names(datTraits),
        yLabels = rownames(moduleTraitCor[c(1:40), ]),
        ySymbols = rownames(moduleTraitCor[c(1:40), ]),
        colorLabels = FALSE,
        colors = my_color,
        textMatrix = textMatrix,
        setStdMargins = TRUE,
        cex.text = 0.5,
        cex.lab.x = 0.5,
        cex.lab.y = 0.5,
        xLabelsAngle = 40,
        zlim = c(-1, 1),
        main = paste("Module-trait relationships")
)

dev.off()


###### Module Membership and Gene Significance

# modNames = substring(names(MEs), 3)
# # Calculate MM
# geneModuleMembership = as.data.frame(signedKME(datExpr, MEs, corFnc = "bicor"));
# MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
# names(geneModuleMembership) = paste("MM", modNames, sep="");
# names(MMPvalue) = paste("p-value.MM", modNames, sep="");
#
# # calculate GS for shoot biomass
# trait_shoot = as.data.frame(datTraits$ShootBiomass);
# names(trait_shoot) = "ShootBiomass"
# geneTraitSignificance_shoot = as.data.frame(cor(datExpr, trait_shoot, use = "p"));
# GSPvalue_shoot = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_shoot), nSamples));
# names(geneTraitSignificance_shoot) = paste("GS.", names(trait_shoot), sep="");
# names(GSPvalue_shoot) = paste("p-value.GS.", names(trait_shoot), sep="");
#
#
# MM_GS <- data.frame(matrix(ncol = 6, nrow = 0))
# colnames(MM_GS) <- c("GeneID","Module","ModuleMembership","pValue_MM","GS_shootBiomass","pvalue_GS_shootBiomass")
#
# for (mdl in c(1:40)){
#   moduleGenes = net$colors==mdl;
#   geneList <- names(moduleGenes[moduleGenes == TRUE])
#
#   moduleData<-data.frame(GeneID = geneList, Module = paste("r-M",mdl,sep=""), ModuleMembership = geneModuleMembership[geneList,paste("MM",mdl,sep="")], pValue_MM = MMPvalue[geneList,paste("p-value.MM",mdl,sep="")],GS_shootBiomass=geneTraitSignificance_shoot[geneList,"GS.ShootBiomass"], pvalue_GS_shootBiomass=GSPvalue_shoot[geneList,"p-value.GS.ShootBiomass"])
#
#   MM_GS<- rbind(MM_GS,moduleData)
# }

MM_GS <- readRDS("MM_GS_rhizobia.RDS")

write.csv(MM_GS,"MM_GS.csv")
