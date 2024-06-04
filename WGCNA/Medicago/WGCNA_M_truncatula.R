#Script to perform Weighted Gene Coexpression Network Analysis using Medicago truncatula Genes

library(WGCNA);
library(gplots);

############ Data Input #########

#set current working directory

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#load averageLogCPM values (calculated using aveLogCPM function)
expData.Avelcpm <- readRDS("expData_AverageLogCPM.RDS")
#33749 genes after filtering for low count genes

datExpr = as.data.frame(t(expData.Avelcpm[,]));
names(datExpr) = rownames(expData.Avelcpm);
rownames(datExpr) = names(expData.Avelcpm);

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

#Next we cluster the samples to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr), method = "average");
#Plot the sample tree
jpeg("Figures/sampleClustering_Mtruncatula.jpeg",width=20,height=10,units="in",res=300,quality=100)
par(cex = 1);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
    cex.axis = 1.5, cex.main = 2)
dev.off()

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


## # Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# #Call the network topology analysis function
#sft <- pickSoftThreshold(datExpr, powerVector = powers,
#                         verbose = 5, corFnc = bicor,
#                          networkType = "signed hybrid",
#                         corOptions = list(maxPOutliers = 0.04))


sft <- readRDS("sft_mtruncatula.RDS")

# # # Plot the results
jpeg("Figures/SoftPowerThreshold.jpeg",width=10,height=7,units="in",res=300,quality=100)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

dev.off()

# we chose B=6 based on scale free topology and mean connectivity
softPower = 6

# enableWGCNAThreads();
# 
# # Network Construction - BLockwise
# net <- blockwiseModules(datExpr,
#                         power =softPower,
#                         maxBlockSize = 35000, #should be larger than num. of genes
#                         #corType = "pearson", #lecture for caveats
#                         corType = "bicor", #lecture for caveats
#                         maxPOutliers = 0.04,
#                         networkType = "signed hybrid",
#                         loadTOM = FALSE, #nothing to load the first time
#                         saveTOMs = FALSE, #save the TOMS to re-load at a later time
#                         deepSplit = 2, #intermediate level
#                         minModuleSize = 30, #a module has to have this many genes
#                         mergeCutHeight = 0.2, #level for merging; suggest 0.2 to start
#                         numericLabels = TRUE,
#                         verbose = 3 )

net <- readRDS("network_mtruncatula.RDS")

table(net$colors)

table(net$unmergedColors)

mergedColors <- labels2colors(net$colors)


# # Show the dendrogram for all genes, with both merged modules colors:
jpeg("Figures//GeneDendrogram_ModuleColors.jpeg",width=20,height=10,units="in",res=300,quality=100)
plotDendroAndColors(net$dendrograms[[1]],mergedColors,
                     "Module colors", main = "Gene dendrogram and module colors",
                     dendroLabels = FALSE, hang = 0.03,
                     addGuide = TRUE, guideHang = 0.05)
 dev.off()


MEs = net$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
jpeg("Figures//ModuleClustering.jpeg",width=20,height=10,units="in",res=300,quality=100)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
dev.off()


# Rename to moduleColors
moduleColors = net$colors;


############ Relate models to experiments######

#read trait data
allTraits <- readRDS("traits_info.RDS")

# Form a data frame analogous to expression data that will hold the traits.
allSamples = rownames(datExpr);
traitRows = match(allSamples, allTraits$strainID);
datTraits = allTraits[traitRows,-1];
rownames(datTraits) = allTraits[traitRows, 1]

#re-order MEs such that similar ones are next to each other
MEs = orderMEs(MEs,greyLast = TRUE, greyName = "ME0")

#moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitCor <- readRDS("moduleTraitCor_Mtruncatula.RDS")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
colnames(moduleTraitPvalue) <- paste("p-value.",colnames(moduleTraitPvalue),sep="")
moduleDF <- cbind(moduleTraitCor,moduleTraitPvalue)
#write.csv(moduleDF, "module-trait-correlation_Mtruncatula.csv")


#plot
library(rcartocolor)

my_color <- colorpanel(150,"#cceeff", "white", "#ffccee")

png("Figures//ModuleTraitRelationship_Mtruncatula.png",width=8,height=13,units="in",res=600)
# Will display correlations and their p-values
textMatrix=  paste(signif(moduleTraitCor, 2), " (",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = paste("p-",names(MEs),sep=""),
               ySymbols = paste("p-",names(MEs),sep=""),
               colorLabels = FALSE,
               colors = my_color,
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.5,
               cex.lab.x = 0.5,
               cex.lab.y = 0.5,
               xLabelsAngle = 40,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))

dev.off()



###### Module Membership and Gene Significance
# 
#  modNames = substring(names(MEs), 3)
# #Calculate MM
#  geneModuleMembership = as.data.frame(signedKME(datExpr, MEs, corFnc = "bicor"));
#  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
#  names(geneModuleMembership) = paste("MM", modNames, sep="");
#  names(MMPvalue) = paste("p-value.MM", modNames, sep="");
#  
# #calculate GS for shoot biomass
#  trait_shoot = as.data.frame(datTraits$ShootBiomass);
#  names(trait_shoot) = "ShootBiomass"
#  geneTraitSignificance_shoot = as.data.frame(cor(datExpr, trait_shoot, use = "p"));
#  GSPvalue_shoot = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_shoot), nSamples));
#  names(geneTraitSignificance_shoot) = paste("GS.", names(trait_shoot), sep="");
#  names(GSPvalue_shoot) = paste("p-value.GS.", names(trait_shoot), sep="");
# 
# MM_GS <- data.frame(matrix(ncol = 6, nrow = 0))
# colnames(MM_GS) <- c("GeneID","Module","ModuleMembership","pValue_MM","GS_shootBiomass","pvalue_GS_shootBiomass")
# 
#  for (mdl in c(1:65)){
#     moduleGenes = net$colors==mdl;
#     geneList <- names(moduleGenes[moduleGenes == TRUE])
#  
#     moduleData<-data.frame(GeneID = geneList, Module = paste("p-M",mdl,sep=""), ModuleMembership = geneModuleMembership[geneList,paste("MM",mdl,sep="")], pValue_MM = MMPvalue[geneList,paste("p-value.MM",mdl,sep="")],GS_shootBiomass=geneTraitSignificance_shoot[geneList,"GS.ShootBiomass"], pvalue_GS_shootBiomass=GSPvalue_shoot[geneList,"p-value.GS.ShootBiomass"])
#  
#     MM_GS<- rbind(MM_GS,moduleData)
#   }

MM_GS <- readRDS("MM_GS_mtruncatula.RDS")

#write.csv(MM_GS,"MM_GS.csv")

