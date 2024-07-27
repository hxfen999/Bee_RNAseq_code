library(DESeq2)
library(edgeR)
library(ggplot2)
library(gplots)
library( "RColorBrewer" )
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()


setwd("C:\\Hxfen\\25.BeeTranscriptome\\github_hxfen")


Gene <- read.csv("gene_count_matrix.csv", header=T)

ID.info <- read.table("ID_info.txt", header=T)
ID.info
Gene.count <- Gene[, c(2:28)]
row.names(Gene.count) <- Gene[,1]
colnames(Gene.count) <- ID.info[,1]
#head(Gene.count)
#head(Gene)
##Filter out the low expression data
group11 <- factor(ID.info[,4])
y <- DGEList(counts=Gene.count, group=group11)
keep <- rowSums(cpm(y) > 0.1) > 5; 
y <- y[keep, ];
dim(y$counts)
Gene.count <- y$counts
head(Gene.count)

Gene.count["MSTRG.5907",1:15] <- 0 #setting the expression value of Toy gene as 0 in the knockout individuals. 


dds <- DESeqDataSetFromMatrix(countData = as.matrix(Gene.count), colData = ID.info, design =~ condition + day)
dds2 <- DESeq(dds)
head(assay(dds2))


rld <- rlog(dds2, blind=FALSE)
dataExpr <- assay(rld)

head(dataExpr)
dim(dataExpr)
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                 max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dataExpr["MSTRG.5907",]
dim(dataExpr)
dataExpr <- as.data.frame(t(dataExprVar))

nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
head(dataExpr)[,1:8]
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.85,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")

power = sft$powerEstimate
power

#if (is.na(power)){
#  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
#          ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
#          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
#          ifelse(type == "unsigned", 6, 12))))
#}


dataExpr[] <- lapply(dataExpr, as.numeric)
cor = WGCNA::cor
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)
cor = stats::cor
table(net$colors)


moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

MEs = net$MEs

MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)


plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

nGenes = ncol(dataExpr);
nSamples = nrow(dataExpr);

MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

dataTraits <- ID.info[,2:3]
dataTraits[,1] <- c(rep(0,15),rep(1,12))
dataTraits[,2] <- c(rep(4,4),rep(3,4),rep(2,3),rep(1,4),rep(4,3),rep(3,3),rep(2,3),rep(1,3))

moduleTraitCor = cor(MEs, dataTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(dataTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

condition = as.data.frame(dataTraits$condition);
names(condition) = "condition";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(dataExpr, condition, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(condition), sep="");
names(GSPvalue) = paste("p.GS.", names(condition), sep="");

module = "lightcyan"   #selecting the module containing target gene (Toy gene) and the related genes.
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

names(dataExpr)
names(dataExpr)[moduleColors=="lightcyan"]

dissTOM = 1-TOMsimilarityFromExpr(dataExpr, power = 12); ## time-consuming calculation

plotTOM = dissTOM^7;
diag(plotTOM) = NA;
TOM = 1 - dissTOM 

modules = c("lightcyan");
probes = names(dataExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

modGenes = probes[(moduleColors==module)];
modGenes

modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
# cyt = exportNetworkToCytoscape(modTOM,
#                                edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), "0.05.txt", sep=""),
#                                nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), "0.05.txt", sep=""),
#                                weighted = TRUE,
#                                threshold = 0.05,
#                                nodeNames = modProbes,
#                                altNodeNames = modGenes,
#                                nodeAttr = moduleColors[inModule])



# Export the network data in a format readable and displayable by the VisANT software.
exportNetworkToVisANT(modTOM, file = paste("VisANTInput-", paste(modules, collapse="-"), ".txt", sep=""),
                      weighted = TRUE, threshold = 0)
exportNetworkToVisANT(modTOM, file = paste("VisANTInput-", paste(modules, collapse="-"), "0.05.txt", sep=""),
                      weighted = TRUE, threshold = 0.05)









