library(DESeq2)
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(amap)
library("genefilter")

setwd("C:\\hxfen\\25.BeeTranscriptome\\github_hxfen")  #set the diroctery containing the expression data as working path. 

Gene <- read.csv("gene_count_matrix.csv", header=T)
#Gene <- read.csv("gene_count.txt", header=T)

ID.info <- read.table("ID_info.txt", header=T)
ID.info
Gene.count <- Gene[, c(2:28)]
row.names(Gene.count) <- Gene[,1]
head(Gene.count)

##Filter out the low expression data
group11 <- factor(ID.info[,4])
y <- DGEList(counts=Gene.count, group=group11)
keep <- rowSums(cpm(y) > 0.1) > 5; 
y <- y[keep, ];
dim(y$counts)
Gene.count <- y$counts
head(Gene.count)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(Gene.count), colData = ID.info, design =~ condition + day)
dds2 <- DESeq(dds)
head(assay(dds2))

#pdf(file = './Unnormalized_normalized_transcript_expression.pdf',width = 12, height = 6)
par(mfrow=c(2,1))
boxplot(log2(assay(dds2)+1))
normalized_counts <- counts(dds2, normalized=TRUE)
head(normalized_counts)
boxplot(log2(counts(dds2, normalized=TRUE)+1))
#dev.off()

vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$condition, vsd$sample_id, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = 'sample euclidean distance')

pcaData <- plotPCA(rld, intgroup = c( "condition", "day"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar
#pdf(file = './res/PCA_by_ggplot.pdf',width = 8,height = 4)
ggplot(pcaData, aes(x = PC1, y = PC2, color = day, shape = condition))  + 
  geom_point(size =3)  + 
  xlab(paste0("PC1: ", percentVar[1], "% variance"))  +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))  +
  coord_fixed()  +
  ggtitle("PCA with RLD data")

######For the 6th day larvae############
Y_Gene.count <- Gene.count[,c("L6E1","L6E2","L6E3","L6E4","L6C1","L6C2","L6C3")]
Y_ID.info <- ID.info[c(12:15,25:27),]
Y_ID.info
dim(Y_Gene.count)
head(Y_Gene.count)
Y_dds <- DESeqDataSetFromMatrix(countData = as.matrix(Y_Gene.count), colData = Y_ID.info, design =~ condition)
Y_dds2 <- DESeq(Y_dds)
Y_rld <- rlog( Y_dds2, blind=FALSE )
Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("condition", "Editted", "Control"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig
Y_resSig_lfc1 <- subset(Y_resSig, abs(log2FoldChange) > 1 )
Y_resSig_lfc1
Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR
cut_off_logFC = 1      #set fold change

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)

#pdf(file = './Youch_volcano.pdf',width = 4.8, height = 4.5)
p <- ggplot(
  # data, point, color
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.6, size=4) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # Auxiliary line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # Coordinate axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  ylim(0, 45)+
  xlim(-30, 30)+
  theme_bw()+
  # legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p


######For the 6th day larvae############
Y_Gene.count <- Gene.count[,c("A1C1","A1C2","A1C3","P3C1","P3C2","P3C3","PP2C1","PP2C2","PP2C3","L6C1","L6C2","L6C3")]
Y_ID.info <- ID.info[c(16:27),]
Y_ID.info
dim(Y_Gene.count)
head(Y_Gene.count)
Y_dds <- DESeqDataSetFromMatrix(countData = as.matrix(Y_Gene.count), colData = Y_ID.info, design =~ groups)
Y_dds2 <- DESeq(Y_dds)
Y_rld <- rlog( Y_dds2, blind=FALSE )

Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "LarvaC", "PrepupaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig

Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue, 
                          ifelse(dataset$log2FoldChange > 0, 'Up', 'Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Larva.vs.Prepupa.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Larva.vs.Prepupa.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)


Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "LarvaC", "PupaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig

Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR


#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue, 
                          ifelse(dataset$log2FoldChange > 0, 'Up', 'Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Larva.vs.Pupa.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Larva.vs.Pupa.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)



Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "LarvaC", "AdultC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig


Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR


#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue, 
                          ifelse(dataset$log2FoldChange > 0, 'Up', 'Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Larva.vs.Adult.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Larva.vs.Adult.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)



Up_genes_L <- c("MSTRG.8539","MSTRG.8371","MSTRG.3261","MSTRG.7826","MSTRG.12240","MSTRG.4120","MSTRG.6049","MSTRG.11659","MSTRG.13033","MSTRG.9485","MSTRG.3343","MSTRG.9661",
"MSTRG.2079","MSTRG.9706","MSTRG.1523","MSTRG.11921","MSTRG.9521","MSTRG.6246","MSTRG.587","MSTRG.11767","MSTRG.13839","MSTRG.160","MSTRG.12828","MSTRG.952",
"MSTRG.455","MSTRG.8319","MSTRG.12630","MSTRG.10511","MSTRG.10241","MSTRG.8166","MSTRG.11696","MSTRG.4960","MSTRG.7895","MSTRG.8451","MSTRG.4855","MSTRG.9537",
"MSTRG.2685","MSTRG.7636","MSTRG.11672","MSTRG.12383","MSTRG.1757","MSTRG.13387","MSTRG.3637","MSTRG.12395","MSTRG.5324","MSTRG.11141","MSTRG.12776","MSTRG.6851",
"MSTRG.4924","MSTRG.10501","MSTRG.3023","MSTRG.13969","MSTRG.5967","MSTRG.7548","MSTRG.13652","MSTRG.12755","MSTRG.359","MSTRG.1507","MSTRG.12272","MSTRG.8267",
"MSTRG.13666","MSTRG.1561","MSTRG.5531","MSTRG.11053","MSTRG.10185","MSTRG.8672","MSTRG.7835","MSTRG.5639","MSTRG.8086","MSTRG.13676","MSTRG.2752","MSTRG.252",
"MSTRG.10631","MSTRG.2086","MSTRG.8297","MSTRG.5666","MSTRG.6170","MSTRG.5820","MSTRG.3772","MSTRG.2979","MSTRG.11831","MSTRG.1961","MSTRG.7819","MSTRG.3006",
"MSTRG.3124","MSTRG.4669","MSTRG.185","MSTRG.12880","MSTRG.12924","MSTRG.4325","MSTRG.10239","MSTRG.6852","MSTRG.11448","MSTRG.7345","MSTRG.2646","MSTRG.405",
"MSTRG.11695","MSTRG.8921","MSTRG.6399","MSTRG.12535","MSTRG.6378","MSTRG.8613","MSTRG.13503","MSTRG.2185","MSTRG.5338","MSTRG.7231","MSTRG.4259","MSTRG.12855",
"MSTRG.8730","MSTRG.4878","MSTRG.8525","MSTRG.1857","MSTRG.161","MSTRG.4350","MSTRG.13624","MSTRG.5423","MSTRG.7655","MSTRG.3833","MSTRG.7904","MSTRG.11551",
"MSTRG.7905","MSTRG.8378","MSTRG.9085","MSTRG.10502","MSTRG.13863","MSTRG.10136","MSTRG.13425","MSTRG.2763","MSTRG.3447","MSTRG.2976","MSTRG.12523","MSTRG.12756",
"MSTRG.11951","MSTRG.5648","MSTRG.13025","MSTRG.1847","MSTRG.6563","MSTRG.11929","MSTRG.13738","MSTRG.1194","MSTRG.9271","MSTRG.11626","MSTRG.3303","MSTRG.3188",
"MSTRG.1747","MSTRG.953","MSTRG.12817","MSTRG.253","MSTRG.251","MSTRG.3917","MSTRG.8352","MSTRG.8553","MSTRG.3706","MSTRG.4876","MSTRG.9313","MSTRG.5532",
"MSTRG.5874","MSTRG.7672","MSTRG.3753","MSTRG.2789","MSTRG.11844","MSTRG.10729","MSTRG.3773","MSTRG.13631","MSTRG.8089","MSTRG.7709","MSTRG.9453","MSTRG.10339",
"MSTRG.4485","MSTRG.9431","MSTRG.13963","MSTRG.13951","MSTRG.13319","MSTRG.5146","MSTRG.2912","MSTRG.7160","MSTRG.9922","MSTRG.5083","MSTRG.8622","MSTRG.11948",
"MSTRG.1629","MSTRG.11206","MSTRG.2751","MSTRG.4699","MSTRG.1229","MSTRG.13016","MSTRG.9990","MSTRG.4349","MSTRG.7579","MSTRG.3658","MSTRG.9274","MSTRG.7343",
"MSTRG.4834","MSTRG.1960","MSTRG.10143","MSTRG.7861","MSTRG.6675","MSTRG.9488","MSTRG.9333","MSTRG.13383","MSTRG.10692","MSTRG.2831","MSTRG.10195","MSTRG.11942",
"MSTRG.7762","MSTRG.13674","MSTRG.6209","MSTRG.8087","MSTRG.5866","MSTRG.10925","MSTRG.13020","MSTRG.310","MSTRG.13315","MSTRG.10472","MSTRG.1291","MSTRG.11959",
"MSTRG.8266","MSTRG.7602","MSTRG.1636","MSTRG.10442","MSTRG.4283","MSTRG.8067","MSTRG.13295","MSTRG.10045","MSTRG.5192","MSTRG.10905","MSTRG.618","MSTRG.8811",
"MSTRG.10525","MSTRG.8197","MSTRG.8134","MSTRG.5087","MSTRG.5702","MSTRG.2212","MSTRG.3308","MSTRG.13687","MSTRG.2975","MSTRG.9851","MSTRG.11969","MSTRG.8082",
"MSTRG.7353","MSTRG.10487","MSTRG.2633","MSTRG.10628","MSTRG.3584","MSTRG.4126","MSTRG.3369","MSTRG.925","MSTRG.3702","MSTRG.5236","MSTRG.1548","MSTRG.2750",
"MSTRG.9416","MSTRG.13180","MSTRG.12241","MSTRG.13938","MSTRG.784","MSTRG.4544","MSTRG.300","MSTRG.8088","MSTRG.221","MSTRG.5705","MSTRG.13659","MSTRG.3240",
"gene-LOC102654520","MSTRG.3654","gene-LOC724425","gene-LOC100578668","MSTRG.1430","MSTRG.11429","MSTRG.11171","MSTRG.7814","MSTRG.6957","MSTRG.12026",
"MSTRG.10530","MSTRG.398","MSTRG.9175","MSTRG.7346","MSTRG.6917","MSTRG.8169","MSTRG.7102","MSTRG.11212","MSTRG.2097","MSTRG.1572","MSTRG.11832","MSTRG.758",
"MSTRG.3381","MSTRG.9209","MSTRG.9828","MSTRG.9070","MSTRG.11142","MSTRG.7583","MSTRG.2375","MSTRG.12789","MSTRG.12595","MSTRG.155","MSTRG.2790","MSTRG.10156",
"MSTRG.6179","MSTRG.8912","MSTRG.12770","MSTRG.2087","MSTRG.8772","MSTRG.1841","MSTRG.13482","MSTRG.13047","MSTRG.10471","MSTRG.3771","MSTRG.12367",
"MSTRG.12551","MSTRG.7550","MSTRG.12816","MSTRG.7334","MSTRG.7106","MSTRG.13084","MSTRG.8137","MSTRG.8464","MSTRG.8975","MSTRG.13667","MSTRG.12251",
"MSTRG.11673","MSTRG.2796","MSTRG.10519","MSTRG.11094","MSTRG.2374","MSTRG.9654","MSTRG.13701","MSTRG.10913","MSTRG.9883","MSTRG.7176","MSTRG.12239",
"MSTRG.11926","MSTRG.393","MSTRG.7226","MSTRG.2245","MSTRG.7315","MSTRG.1597","MSTRG.12110","MSTRG.11181","MSTRG.11616","MSTRG.13664","MSTRG.6866",
"MSTRG.3755","gene-LOC725106","MSTRG.5880","MSTRG.8671")



heatmap( assay(Y_rld)[Up_genes_L, ], scale="row",
trace="none", Rowv=NA, Colv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))





#######PP##########

Y_Gene.count <- Gene.count[,c("A1C1","A1C2","A1C3","P3C1","P3C2","P3C3","PP2C1","PP2C2","PP2C3","L6C1","L6C2","L6C3")]
Y_ID.info <- ID.info[c(16:27),]
Y_ID.info
dim(Y_Gene.count)
head(Y_Gene.count)
Y_dds <- DESeqDataSetFromMatrix(countData = as.matrix(Y_Gene.count), colData = Y_ID.info, design =~ groups)
Y_dds2 <- DESeq(Y_dds)
Y_rld <- rlog( Y_dds2, blind=FALSE )


Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "PrepupaC", "LarvaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig

Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue , 
                          ifelse(dataset$log2FoldChange> 0 ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)

write.table(file="Down_Prepupa.vs.Larva.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Prepupa.vs.Larva.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)



Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "PrepupaC", "PupaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig

Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue, 
                          ifelse(dataset$log2FoldChange > 0,'Up','Down'),
                          'Stable')

dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Prepupa.vs.Pupa.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Prepupa.vs.Pupa.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)



Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "PrepupaC", "AdultC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig

Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR


#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue, 
                          ifelse(dataset$log2FoldChange > 0, 'Up', 'Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Prepupa.vs.Adult.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Prepupa.vs.Adult.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)


Up_genes_PP <- c("MSTRG.229","MSTRG.7706","MSTRG.4870","MSTRG.10234","MSTRG.10857","MSTRG.782","MSTRG.2736","MSTRG.6858","MSTRG.3849","MSTRG.7157","MSTRG.2337","MSTRG.5712",
"MSTRG.2447","MSTRG.7214","MSTRG.2291","MSTRG.1408","MSTRG.3712","MSTRG.3098","MSTRG.11238","MSTRG.12075","MSTRG.8616","MSTRG.3677","MSTRG.12311","MSTRG.3983",
"MSTRG.6845","MSTRG.4224","MSTRG.1746","MSTRG.3139","MSTRG.7700","MSTRG.13472","MSTRG.9092","MSTRG.1293","MSTRG.12543","MSTRG.11082","MSTRG.5812","MSTRG.6639",
"MSTRG.6210","MSTRG.13173","MSTRG.3493","MSTRG.6388","MSTRG.5237","MSTRG.2162","MSTRG.11379","MSTRG.385","MSTRG.4590","MSTRG.10116","MSTRG.11609","MSTRG.12049",
"MSTRG.394","MSTRG.8976","MSTRG.13523","MSTRG.5447","MSTRG.13001","MSTRG.786","MSTRG.7273","MSTRG.8149","MSTRG.1147","MSTRG.608","MSTRG.53","MSTRG.11037",
"MSTRG.8351","MSTRG.6648","MSTRG.4809","MSTRG.9118","MSTRG.4027","MSTRG.10723","MSTRG.2468","MSTRG.6392","MSTRG.7998","MSTRG.10505","MSTRG.2135","MSTRG.9364",
"MSTRG.10053","MSTRG.1199","MSTRG.9659","MSTRG.3842","MSTRG.8108","MSTRG.4653","MSTRG.6133","MSTRG.10041","MSTRG.1162","MSTRG.5375","MSTRG.1116","MSTRG.348",
"MSTRG.13257","MSTRG.4336","MSTRG.6018","MSTRG.12868","MSTRG.11865","MSTRG.1295","MSTRG.10352","MSTRG.263","MSTRG.2559","MSTRG.11338","MSTRG.12634","MSTRG.11298",
"MSTRG.7991","MSTRG.13317","MSTRG.12810","MSTRG.10335","MSTRG.7355","MSTRG.13617","MSTRG.2961","MSTRG.6099","MSTRG.2288","MSTRG.13360","MSTRG.3128","MSTRG.3675",
"MSTRG.3080","MSTRG.11782","MSTRG.10640","MSTRG.5775","MSTRG.3508","MSTRG.6334","MSTRG.10878","MSTRG.7816","MSTRG.12379","MSTRG.4690","MSTRG.10475","MSTRG.4521",
"MSTRG.9316","MSTRG.1328","MSTRG.11336","MSTRG.6996","MSTRG.12475","MSTRG.2854","MSTRG.10290","MSTRG.10653","MSTRG.9619","MSTRG.12093","MSTRG.4083","MSTRG.2403",
"MSTRG.12303","MSTRG.5724","MSTRG.3652","MSTRG.10541","MSTRG.12357","MSTRG.6552","MSTRG.5393","MSTRG.3617","MSTRG.12442","MSTRG.8760","MSTRG.5920","MSTRG.12066",
"MSTRG.5100","MSTRG.12481","MSTRG.9066","MSTRG.4030","MSTRG.12531","MSTRG.11368","MSTRG.3230","MSTRG.6025","MSTRG.12134","MSTRG.10198","MSTRG.12384","MSTRG.2244",
"MSTRG.6993","MSTRG.4598","MSTRG.11567","MSTRG.7403","MSTRG.1646","MSTRG.7430","MSTRG.8521","MSTRG.3120","MSTRG.1379","MSTRG.8498","MSTRG.10959","MSTRG.11823",
"MSTRG.10722","MSTRG.1176","MSTRG.12527","MSTRG.5887","MSTRG.3659","MSTRG.2357","MSTRG.8945","MSTRG.12104","MSTRG.6848","MSTRG.6521","MSTRG.9226","MSTRG.11683",
"MSTRG.2988","MSTRG.11922","MSTRG.9406","MSTRG.7846","MSTRG.11915","MSTRG.11892","MSTRG.9566","MSTRG.2811","MSTRG.12228","MSTRG.2659","MSTRG.2848","MSTRG.3739",
"MSTRG.5870","MSTRG.12090","MSTRG.1172","MSTRG.816","MSTRG.10813","MSTRG.11591","MSTRG.3137","MSTRG.8528","MSTRG.6131","MSTRG.64","MSTRG.7977","MSTRG.1805",
"MSTRG.5380","MSTRG.215","MSTRG.4415","MSTRG.3564","MSTRG.3668","MSTRG.7201","MSTRG.151","MSTRG.10071","MSTRG.6098","MSTRG.9921","MSTRG.4969","MSTRG.6661",
"MSTRG.3345","MSTRG.10581","MSTRG.6916","MSTRG.2982","MSTRG.397","MSTRG.8537","MSTRG.12562","MSTRG.224","MSTRG.3623","MSTRG.8017","MSTRG.793","MSTRG.8536",
"MSTRG.3812","MSTRG.417","MSTRG.5177","MSTRG.4549","MSTRG.11740","MSTRG.13823","MSTRG.9344","MSTRG.7541","MSTRG.11884","MSTRG.10069","MSTRG.2551","MSTRG.4338",
"MSTRG.8421","MSTRG.8998","MSTRG.44","MSTRG.2892","MSTRG.2290","MSTRG.10851","MSTRG.4932","MSTRG.5962","MSTRG.3452","MSTRG.1333","MSTRG.4633","MSTRG.5122",
"MSTRG.3494","MSTRG.4982","MSTRG.9211","MSTRG.4869","MSTRG.3158","MSTRG.4884","MSTRG.10599","MSTRG.10171","MSTRG.867","MSTRG.10728","MSTRG.12443","MSTRG.10274",
"MSTRG.2667","MSTRG.9571","MSTRG.10681","MSTRG.4188","MSTRG.2287","MSTRG.13657","MSTRG.11074","MSTRG.11594","MSTRG.1177","MSTRG.11073","MSTRG.2024","MSTRG.3150",
"MSTRG.10200","MSTRG.5328","MSTRG.574","MSTRG.6654","MSTRG.9025","MSTRG.10366","MSTRG.12885","MSTRG.2519","MSTRG.10651","MSTRG.9318","MSTRG.4593","MSTRG.1513",
"MSTRG.8052","MSTRG.2539","MSTRG.2391","MSTRG.4039","MSTRG.7502","MSTRG.7534","MSTRG.7851","MSTRG.5391","MSTRG.807","MSTRG.8308","MSTRG.7710","MSTRG.6391",
"MSTRG.7664","MSTRG.5040","MSTRG.9275","MSTRG.12458","MSTRG.6247","MSTRG.3428","MSTRG.13778","MSTRG.3912","MSTRG.5182","MSTRG.1192","MSTRG.2196","MSTRG.6156",
"MSTRG.1371","MSTRG.2581","MSTRG.1662","MSTRG.12659","MSTRG.11761","MSTRG.299","MSTRG.7423","MSTRG.5256","MSTRG.7720","MSTRG.2404","MSTRG.13861","MSTRG.3045",
"MSTRG.8544","MSTRG.11717","MSTRG.5088","MSTRG.5755","MSTRG.11163","MSTRG.13288","MSTRG.13518","MSTRG.415","MSTRG.13686","MSTRG.3202","MSTRG.8362","MSTRG.5691",
"MSTRG.9930","MSTRG.11862","MSTRG.12754","MSTRG.3145","MSTRG.10991","MSTRG.5286","MSTRG.10291","MSTRG.12351","MSTRG.2289","MSTRG.8611","MSTRG.13658","MSTRG.11339",
"MSTRG.7711","MSTRG.13584","MSTRG.825","MSTRG.1955","MSTRG.5381","MSTRG.2345","MSTRG.349","MSTRG.5600","MSTRG.5044","MSTRG.8412","MSTRG.4380","MSTRG.7684",
"MSTRG.2847","MSTRG.2293","MSTRG.10236","MSTRG.10193","MSTRG.5386","MSTRG.2195","MSTRG.2147")


heatmap( assay(Y_rld)[Up_genes_PP, ], scale="row",
trace="none", Rowv=NA, Colv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))




#######Pupa##########

Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "PupaC", "LarvaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig

Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR
cut_off_logFC = 1      #set fold change

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Pupa.vs.Larva.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Pupa.vs.Larva.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)

#pdf(file = './Youch_volcano.pdf',width = 4.8, height = 4.5)
p <- ggplot(
  # data, point, color
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.6, size=4) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # Auxiliary line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # Coordinate axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  ylim(0, 45)+
  xlim(-30, 30)+
  theme_bw()+
  # legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p



Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "PupaC", "PrepupaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig
Y_resSig_lfc1 <- subset(Y_resSig, abs(log2FoldChange) > 1 )
Y_resSig_lfc1
Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR
cut_off_logFC = 1      #set fold change

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Pupa.vs.Prepupa.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Pupa.vs.Prepupa.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)

#pdf(file = './Youch_volcano.pdf',width = 4.8, height = 4.5)
p <- ggplot(
  # data, point, color
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.6, size=4) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # Auxiliary line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # Coordinate axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  ylim(0, 45)+
  xlim(-30, 30)+
  theme_bw()+
  # legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p




Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "PupaC", "AdultC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig
Y_resSig_lfc1 <- subset(Y_resSig, abs(log2FoldChange) > 1 )
Y_resSig_lfc1
Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR
cut_off_logFC = 1      #set fold change

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Pupa.vs.Adult.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Pupa.vs.Adult.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)

#pdf(file = './Youch_volcano.pdf',width = 4.8, height = 4.5)
p <- ggplot(
  # data, point, color
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.6, size=4) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # Auxiliary line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # Coordinate axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  ylim(0, 45)+
  xlim(-30, 30)+
  theme_bw()+
  # legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p



Up_genes_P <- c("MSTRG.12981","MSTRG.12982","MSTRG.13603","MSTRG.6046","MSTRG.9726","MSTRG.6113","MSTRG.8713","MSTRG.7891","MSTRG.1066","MSTRG.7751","MSTRG.549","MSTRG.5290","MSTRG.13357",
"MSTRG.11619","MSTRG.13619","MSTRG.9598","MSTRG.3298","MSTRG.9982","MSTRG.10659","MSTRG.8706","MSTRG.2062","MSTRG.3955","MSTRG.2944","MSTRG.3014","MSTRG.2830","MSTRG.10064",
"MSTRG.1290","MSTRG.11876","MSTRG.4434","MSTRG.8828","MSTRG.2884","MSTRG.11214","MSTRG.13184","MSTRG.8645","MSTRG.1485","MSTRG.2461","MSTRG.2201","MSTRG.10214","MSTRG.4840",
"MSTRG.8655","MSTRG.11540","MSTRG.11743","MSTRG.12987","MSTRG.12456","MSTRG.3318","MSTRG.11116","MSTRG.8677","MSTRG.12581","MSTRG.13371","MSTRG.4731","MSTRG.1094","MSTRG.5118",
"MSTRG.113","MSTRG.3219","MSTRG.145","MSTRG.10327","MSTRG.2527","MSTRG.3218","MSTRG.10997","MSTRG.11102","MSTRG.2018","MSTRG.1610","MSTRG.4271","MSTRG.12045","MSTRG.8311",
"MSTRG.7215","MSTRG.4159","MSTRG.1771","MSTRG.9004","MSTRG.12656","MSTRG.4214","MSTRG.5079","MSTRG.881","MSTRG.1266","MSTRG.6227","MSTRG.2945","MSTRG.8697","MSTRG.12968",
"MSTRG.11131","MSTRG.1303","MSTRG.8327","MSTRG.10736","MSTRG.13370","MSTRG.164","MSTRG.10007","MSTRG.13309","MSTRG.9373","MSTRG.5737","MSTRG.1812","MSTRG.173","MSTRG.8399",
"MSTRG.12532","MSTRG.8454","MSTRG.2526","MSTRG.623","MSTRG.13852","MSTRG.12977","MSTRG.7138","MSTRG.7219","MSTRG.13708","MSTRG.13008","MSTRG.12993","MSTRG.8673","MSTRG.6279",
"MSTRG.10473","MSTRG.1240","MSTRG.13620","MSTRG.11568","MSTRG.5076","MSTRG.13511","MSTRG.5839","MSTRG.6577","MSTRG.13348","MSTRG.11135","MSTRG.5248","MSTRG.11671","MSTRG.12563",
"MSTRG.8678","MSTRG.9041","MSTRG.10562","MSTRG.13407","MSTRG.8317","MSTRG.8453","MSTRG.5081","MSTRG.6469","MSTRG.4677","MSTRG.8044","MSTRG.12969","MSTRG.11316","MSTRG.1566",
"MSTRG.12929","MSTRG.2417","MSTRG.10225","MSTRG.12569","MSTRG.12978","gene-LOC100576122","MSTRG.9383","MSTRG.11667","MSTRG.12984","MSTRG.13605","MSTRG.5944","MSTRG.8789",
"MSTRG.2525","MSTRG.7542","MSTRG.8151","MSTRG.9636","MSTRG.7323","MSTRG.628","MSTRG.10354","MSTRG.7496","MSTRG.11744","MSTRG.12116","MSTRG.2063","MSTRG.3807","MSTRG.11511",
"MSTRG.718","MSTRG.12994","MSTRG.12990","MSTRG.8826","MSTRG.12628","MSTRG.8487","MSTRG.3387","MSTRG.1567","MSTRG.8289","MSTRG.12390","MSTRG.13042","MSTRG.4949","MSTRG.11596",
"MSTRG.4752","MSTRG.9456","MSTRG.6981","MSTRG.3016","MSTRG.8838","MSTRG.9529","MSTRG.3388","MSTRG.10733","MSTRG.7063","MSTRG.12561","MSTRG.2415","MSTRG.49","MSTRG.6576",
"MSTRG.2405","MSTRG.10734","MSTRG.2476","gene-LOC107965541","MSTRG.6437","MSTRG.10439","MSTRG.4163","MSTRG.13513","MSTRG.11734","MSTRG.1562","MSTRG.8452","MSTRG.8708","MSTRG.837",
"MSTRG.11317","MSTRG.7012","MSTRG.2641","MSTRG.5119","MSTRG.5745","MSTRG.11641","MSTRG.8656","MSTRG.2419","MSTRG.12985","MSTRG.9043","MSTRG.9665","MSTRG.12986","MSTRG.13012",
"MSTRG.5941","MSTRG.3993","MSTRG.13337","MSTRG.11127","MSTRG.10965","MSTRG.1364","MSTRG.13561","MSTRG.12126","MSTRG.10555","MSTRG.2200","MSTRG.178","MSTRG.8592","MSTRG.9582",
"MSTRG.1014","MSTRG.13181","MSTRG.3809","MSTRG.4164","MSTRG.4272","MSTRG.9570")


heatmap.2( assay(Y_rld)[Up_genes_P, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

Down_genes_P <- c("MSTRG.8913","MSTRG.11383","MSTRG.9642","MSTRG.9774","MSTRG.875","MSTRG.11767","MSTRG.4236","MSTRG.3751","MSTRG.11250","MSTRG.3752","gene-Kr","MSTRG.5418",
"MSTRG.3006","MSTRG.284","MSTRG.7934","MSTRG.8499","MSTRG.10462","MSTRG.9559","MSTRG.10321","MSTRG.12396","MSTRG.7904","MSTRG.7905","MSTRG.13304","MSTRG.12276",
"MSTRG.12523","MSTRG.10967","MSTRG.4862","MSTRG.12066","MSTRG.5221","MSTRG.11849","MSTRG.7883","MSTRG.8122","MSTRG.9534","MSTRG.7579","MSTRG.11689","MSTRG.12112",
"MSTRG.764","MSTRG.5222","MSTRG.39","MSTRG.6341","MSTRG.324","MSTRG.12895","MSTRG.11546","MSTRG.1688","MSTRG.1548","MSTRG.11820","MSTRG.3240","MSTRG.280",
"MSTRG.8964","MSTRG.7956","gene-LOC102656209","MSTRG.7603","MSTRG.10156","MSTRG.10929","MSTRG.3750","MSTRG.13575","MSTRG.10117","MSTRG.10526","MSTRG.11825",
"MSTRG.995","MSTRG.11024","MSTRG.11181","MSTRG.264")

heatmap.2( assay(Y_rld)[Down_genes_P, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


#######Adult##########

Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "AdultC", "LarvaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig
Y_resSig_lfc1 <- subset(Y_resSig, abs(log2FoldChange) > 1 )
Y_resSig_lfc1
Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR
cut_off_logFC = 1      #set fold change

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Adult.vs.Larva.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Adult.vs.Larva.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)

#pdf(file = './Youch_volcano.pdf',width = 4.8, height = 4.5)
p <- ggplot(
  # data, point, color
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.6, size=4) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # Auxiliary line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # Coordinate axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  ylim(0, 45)+
  xlim(-30, 30)+
  theme_bw()+
  # legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p



Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "AdultC", "PrepupaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig
Y_resSig_lfc1 <- subset(Y_resSig, abs(log2FoldChange) > 1 )
Y_resSig_lfc1
Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR
cut_off_logFC = 1      #set fold change

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Adult.vs.Prepupa.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Adult.vs.Prepupa.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)

#pdf(file = './Youch_volcano.pdf',width = 4.8, height = 4.5)
p <- ggplot(
  # data, point, color
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.6, size=4) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # Auxiliary line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # Coordinate axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  ylim(0, 45)+
  xlim(-30, 30)+
  theme_bw()+
  # legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p




Y_res <- results(Y_dds2, pAdjustMethod = "BH", contrast = c("groups", "AdultC", "PupaC"))  
Y_resSig <- subset(Y_res, padj < 0.05)
Y_resSig
Y_resSig_lfc1 <- subset(Y_resSig, abs(log2FoldChange) > 1 )
Y_resSig_lfc1
Y_resSig_order <- Y_resSig[order(Y_resSig$log2FoldChange),]
Y_resSig_orderlessthan0 <- Y_resSig_order[Y_resSig$log2FoldChange<0,]
dim(Y_resSig_orderlessthan0)
Y_resSig_names <- row.names(Y_resSig_order)

Y_rlogMat <- assay(Y_rld)[Y_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
Y_pearson_cor <- as.matrix(cor(Y_rlogMat, method="pearson"))
Y_hc <- hcluster(t(Y_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_Y_Pearsoncor.pdf",width=7,height=7)
heatmap.2(Y_pearson_cor, Rowv=as.dendrogram(Y_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- Y_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #set FDR
cut_off_logFC = 1      #set fold change

#According to the threshold parameters, set the up-regulated gene to 'Up', the down-regulated gene to 'Down', and the no difference gene to 'Stable', and save them in the "change" column.
dataset$change = ifelse(dataset$padj < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                          ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                          'Stable')
dataset.down <- subset(dataset, change == "Down")
dataset.up <- subset(dataset, change == "Up")
row.names(dataset.down)
row.names(dataset.up)
write.table(file="Down_Adult.vs.Pupa.txt", row.names(dataset.down),quote=F, row.names = F, col.names = F)
write.table(file="Up_Adult.vs.Pupa.txt", row.names(dataset.up),quote=F, row.names = F, col.names = F)

#pdf(file = './Youch_volcano.pdf',width = 4.8, height = 4.5)
p <- ggplot(
  # data, point, color
  dataset, aes(x = log2FoldChange, y = -log10(padj), colour=change)) +
  geom_point(alpha=0.6, size=4) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  # Auxiliary line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lty=4,col="black",lwd=0.8) +
  # Coordinate axis
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  ylim(0, 45)+
  xlim(-30, 30)+
  theme_bw()+
  # legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p



Up_genes_A <- c("MSTRG.5387","MSTRG.8871","MSTRG.5566","MSTRG.3615","MSTRG.9686","MSTRG.12284","MSTRG.6804","MSTRG.9642","MSTRG.9611","MSTRG.13067","MSTRG.8343","MSTRG.2035","MSTRG.2744","MSTRG.8590",
"MSTRG.1331","MSTRG.7383","MSTRG.9479","MSTRG.2495","MSTRG.1631","MSTRG.1935","MSTRG.7619","MSTRG.12368","MSTRG.7903","MSTRG.6685","MSTRG.5564","MSTRG.11872","MSTRG.11495","MSTRG.9697",
"MSTRG.1055","MSTRG.57","MSTRG.12149","MSTRG.2323","MSTRG.3804","MSTRG.3591","MSTRG.8771","MSTRG.4431","MSTRG.6422","MSTRG.1131","gene-5-ht7","MSTRG.2664","MSTRG.2386","MSTRG.3444",
"MSTRG.1793","MSTRG.13283","MSTRG.5811","MSTRG.13510","MSTRG.13509","MSTRG.9867","MSTRG.7955","MSTRG.9317","MSTRG.13272","MSTRG.763","MSTRG.9309","MSTRG.11912","MSTRG.6404",
"MSTRG.2429","MSTRG.7638","MSTRG.12161","MSTRG.12699","MSTRG.11931","MSTRG.8601","MSTRG.12356","MSTRG.11228","MSTRG.13076","MSTRG.7549","MSTRG.12288","MSTRG.4396","MSTRG.12437",
"MSTRG.10363","MSTRG.13750","MSTRG.6643","MSTRG.12808","MSTRG.6870","MSTRG.8460","MSTRG.7960","MSTRG.13848","MSTRG.1431","MSTRG.9563","MSTRG.13068","MSTRG.1964","MSTRG.3782",
"MSTRG.3883","MSTRG.13380","MSTRG.2056","MSTRG.4972","MSTRG.3577","MSTRG.12757","MSTRG.2219","MSTRG.5556","MSTRG.4985","MSTRG.8765","MSTRG.11397","MSTRG.4686","MSTRG.13037",
"MSTRG.9495","MSTRG.7979","MSTRG.12601","MSTRG.12463","MSTRG.11705","MSTRG.1238","MSTRG.1278","MSTRG.5557","MSTRG.526","MSTRG.13141","MSTRG.1078","MSTRG.10663","MSTRG.7785",
"MSTRG.6781","MSTRG.11130","MSTRG.7198","MSTRG.10388","MSTRG.11250","MSTRG.10927","MSTRG.1817","MSTRG.13207","MSTRG.4875","MSTRG.4124","MSTRG.3221","MSTRG.11333","MSTRG.7191",
"MSTRG.482","MSTRG.8664","MSTRG.4242","MSTRG.3963","MSTRG.4198","MSTRG.5973","MSTRG.5673","MSTRG.3153","MSTRG.11445","MSTRG.7142","MSTRG.6813","MSTRG.7305","MSTRG.9325",
"MSTRG.2009","MSTRG.1568","MSTRG.826","MSTRG.3791","MSTRG.2311","MSTRG.13046","MSTRG.11954","MSTRG.8892","MSTRG.6924","MSTRG.9040","MSTRG.1261","MSTRG.233","MSTRG.666",
"MSTRG.9100","MSTRG.4113","MSTRG.3616","MSTRG.11215","MSTRG.5538","MSTRG.1314","MSTRG.5601","MSTRG.2609","MSTRG.5193","MSTRG.8058","MSTRG.8125","MSTRG.5616","MSTRG.6562",
"MSTRG.4693","MSTRG.6144","MSTRG.4168","MSTRG.3987","MSTRG.4186","MSTRG.306","MSTRG.4701","MSTRG.5478","MSTRG.2720","gene-LOC409751","MSTRG.6283","MSTRG.1675","MSTRG.3360",
"MSTRG.10358","MSTRG.10867","MSTRG.8211","MSTRG.1766","MSTRG.12636","MSTRG.13024","MSTRG.13721","MSTRG.991","MSTRG.9925","MSTRG.10838","MSTRG.6183","MSTRG.5453","MSTRG.12199",
"MSTRG.3602","MSTRG.8244","MSTRG.313","MSTRG.2317","MSTRG.7958","MSTRG.3872","MSTRG.4477","MSTRG.1850","MSTRG.4288","MSTRG.5548","MSTRG.9433","MSTRG.11128","MSTRG.13261",
"MSTRG.2234","MSTRG.335","MSTRG.4019","MSTRG.1831","MSTRG.3880","MSTRG.6860","MSTRG.9825","MSTRG.1600","MSTRG.5169","MSTRG.7951","MSTRG.7856","MSTRG.10900","MSTRG.9761",
"MSTRG.7494","MSTRG.9747","MSTRG.9559","MSTRG.13124","MSTRG.1713","MSTRG.5606","MSTRG.10961","MSTRG.6931","MSTRG.6880","MSTRG.3727","MSTRG.5153","MSTRG.13981","MSTRG.5224",
"MSTRG.7892","MSTRG.3140","MSTRG.6704","MSTRG.9870","MSTRG.2319","MSTRG.9696","MSTRG.12121","MSTRG.1764","MSTRG.10975","MSTRG.746","MSTRG.1963","MSTRG.11930","MSTRG.12596",
"MSTRG.714","MSTRG.4600","MSTRG.3624","MSTRG.12470","MSTRG.4189","MSTRG.8250","MSTRG.8943","MSTRG.357","MSTRG.5510","MSTRG.5914","MSTRG.6909","MSTRG.3932","MSTRG.7121",
"MSTRG.8007","MSTRG.7601","MSTRG.8467","MSTRG.341","MSTRG.3097","MSTRG.5873","MSTRG.7019","MSTRG.3197","MSTRG.6424","MSTRG.12713","MSTRG.4262","MSTRG.6805","MSTRG.634",
"MSTRG.7556","MSTRG.3094","MSTRG.4495","MSTRG.1337","MSTRG.8213","MSTRG.958","MSTRG.6108","MSTRG.7949","MSTRG.11124","MSTRG.5565","MSTRG.2338","MSTRG.1630","MSTRG.4995",
"MSTRG.13059","MSTRG.1786","MSTRG.2554","MSTRG.8336","MSTRG.13159","MSTRG.7974","MSTRG.1352","MSTRG.13172","MSTRG.5211","MSTRG.3069","MSTRG.13304","MSTRG.2732","MSTRG.4489",
"MSTRG.1090","MSTRG.2644","MSTRG.6974","MSTRG.1054","MSTRG.6927","MSTRG.3894","MSTRG.13293","MSTRG.4617","MSTRG.8187","MSTRG.771","MSTRG.12488","MSTRG.13805","MSTRG.13453",
"MSTRG.8069","MSTRG.8837","MSTRG.10035","MSTRG.1452","MSTRG.9609","MSTRG.11933","MSTRG.9376","MSTRG.9625","MSTRG.1244","MSTRG.3895","MSTRG.7092","MSTRG.1101","MSTRG.2300",
"MSTRG.8926","MSTRG.2158","MSTRG.11863","MSTRG.9762","MSTRG.4679","MSTRG.7283","MSTRG.13149","MSTRG.5160","MSTRG.288","MSTRG.6150","MSTRG.8980","MSTRG.9322","MSTRG.13077",
"MSTRG.13273","MSTRG.4673","MSTRG.5682","MSTRG.7919","MSTRG.12074","MSTRG.6405","MSTRG.10967","MSTRG.13075","MSTRG.5774","MSTRG.5034","MSTRG.3865","MSTRG.3061","MSTRG.5208",
"MSTRG.12822","MSTRG.12162","MSTRG.9327","MSTRG.4087","MSTRG.9562","MSTRG.10744","MSTRG.11420","MSTRG.12296","MSTRG.2142","MSTRG.1965","MSTRG.6942","MSTRG.9530","MSTRG.4721",
"MSTRG.2474","MSTRG.11369","MSTRG.10385","MSTRG.7954","MSTRG.4861","MSTRG.6952","MSTRG.4233","MSTRG.13045","MSTRG.6799","MSTRG.5460","MSTRG.2568","MSTRG.5664","MSTRG.13754",
"MSTRG.7418","MSTRG.10159","MSTRG.5099","MSTRG.3890","MSTRG.13112","MSTRG.11545","MSTRG.1755","MSTRG.13160","MSTRG.8411","MSTRG.1237","MSTRG.5167","MSTRG.8359","MSTRG.40",
"MSTRG.10897","MSTRG.10398","MSTRG.7939","MSTRG.7250","MSTRG.6642","MSTRG.2900","MSTRG.11725","MSTRG.2401","MSTRG.8800","MSTRG.11356","MSTRG.2829","MSTRG.3272","MSTRG.6615",
"MSTRG.3438","MSTRG.3412","MSTRG.6746","MSTRG.2316","MSTRG.5292","MSTRG.11878","MSTRG.13051","MSTRG.2008","MSTRG.13019","MSTRG.3135","MSTRG.3038","MSTRG.10642","MSTRG.13055",
"MSTRG.8696","MSTRG.2870","MSTRG.262","MSTRG.410","MSTRG.13985","MSTRG.7566","MSTRG.11871","MSTRG.8887","MSTRG.11839","gene-LOC107964265","MSTRG.7921","MSTRG.2502",
"MSTRG.6873","MSTRG.7224","MSTRG.11643","MSTRG.2486","MSTRG.5904","MSTRG.12519","MSTRG.12893","MSTRG.5337","MSTRG.13880","MSTRG.12546","MSTRG.7945","MSTRG.8589","MSTRG.12621",
"MSTRG.5956","MSTRG.1465","MSTRG.9748","MSTRG.6068","MSTRG.9933","MSTRG.3057","MSTRG.990","MSTRG.1892","MSTRG.4076","MSTRG.8858","MSTRG.7347","MSTRG.1064","MSTRG.8448",
"MSTRG.4476","MSTRG.11446","MSTRG.8526","MSTRG.5683","MSTRG.3350","MSTRG.4777","MSTRG.11941","MSTRG.1210","MSTRG.8122","MSTRG.3715","MSTRG.6597","MSTRG.11245","MSTRG.9059",
"MSTRG.6117","MSTRG.13258","MSTRG.11934","MSTRG.3873","MSTRG.5200","MSTRG.12139","MSTRG.4988","MSTRG.13305","MSTRG.13276","MSTRG.8962","MSTRG.13061","MSTRG.639","MSTRG.9546",
"MSTRG.7957","MSTRG.4280","MSTRG.11688","MSTRG.3841","MSTRG.8323","MSTRG.6182","MSTRG.3348","MSTRG.7758","MSTRG.10026","MSTRG.10293","MSTRG.12147","MSTRG.1880","MSTRG.3062",
"MSTRG.1423","MSTRG.9036","MSTRG.1174","MSTRG.138","MSTRG.12297","MSTRG.4776","MSTRG.7633","MSTRG.7101","MSTRG.3874","MSTRG.5692","MSTRG.8123","MSTRG.668","MSTRG.1033",
"MSTRG.12603","MSTRG.4353","MSTRG.10688","MSTRG.13572","MSTRG.9805","MSTRG.13336","MSTRG.8888","MSTRG.9058","MSTRG.9576","MSTRG.142","MSTRG.6941","MSTRG.7429","MSTRG.13259",
"MSTRG.13679","MSTRG.13058","MSTRG.1065","MSTRG.12759","MSTRG.13846","MSTRG.9497","MSTRG.1091","MSTRG.13793","MSTRG.2007","MSTRG.1071","MSTRG.7546","MSTRG.3495","MSTRG.3238",
"MSTRG.1243","MSTRG.2029","MSTRG.12761","MSTRG.41","MSTRG.3417","MSTRG.11335","MSTRG.6423","MSTRG.8891","MSTRG.4316","MSTRG.3891","MSTRG.11240","MSTRG.10831","MSTRG.13274",
"MSTRG.9174","MSTRG.680","MSTRG.5554","MSTRG.11409","MSTRG.8288","MSTRG.11974","MSTRG.7689","MSTRG.194","MSTRG.6856","MSTRG.12883","MSTRG.5491","MSTRG.8799","MSTRG.3827",
"MSTRG.6598","MSTRG.1602","MSTRG.12286","MSTRG.4010","MSTRG.5519","MSTRG.13065","MSTRG.3077","MSTRG.1505","MSTRG.5222","MSTRG.13162","MSTRG.11610","MSTRG.966","MSTRG.9906",
"MSTRG.12223","MSTRG.2533","MSTRG.11426","MSTRG.13104","MSTRG.7275","MSTRG.11854","MSTRG.1451","MSTRG.7129","MSTRG.2406","MSTRG.9422","MSTRG.10520","MSTRG.13454","MSTRG.9088",
"MSTRG.3130","MSTRG.5155","MSTRG.2432","MSTRG.6401","MSTRG.12159","MSTRG.2731","MSTRG.4608","MSTRG.6044","MSTRG.1670","MSTRG.9786","MSTRG.2901","MSTRG.5659","MSTRG.3331",
"MSTRG.9124","MSTRG.4497","MSTRG.2610","MSTRG.8041","MSTRG.13570","MSTRG.802","MSTRG.11589","MSTRG.8153","MSTRG.7531","MSTRG.9357","MSTRG.7953","MSTRG.13692","MSTRG.10518",
"MSTRG.1966","MSTRG.13279","MSTRG.11447","MSTRG.6402","MSTRG.1280","MSTRG.12650","MSTRG.13376","MSTRG.13074","MSTRG.6869","MSTRG.3600","MSTRG.9393","MSTRG.1869","MSTRG.12255",
"MSTRG.11147","MSTRG.2444","MSTRG.9991","MSTRG.10747","MSTRG.12845","MSTRG.3724","MSTRG.2155","MSTRG.6759","MSTRG.12821","MSTRG.8679","MSTRG.9749","MSTRG.11973","MSTRG.1404",
"MSTRG.10956","MSTRG.3354","MSTRG.11838","MSTRG.8475","MSTRG.1810","MSTRG.5726","MSTRG.9934","MSTRG.2899","MSTRG.3196","MSTRG.13993","MSTRG.9966","MSTRG.1741","MSTRG.12466",
"MSTRG.2605","MSTRG.11427","MSTRG.6850","MSTRG.4268","MSTRG.10544","MSTRG.9557","MSTRG.11633","MSTRG.2535","MSTRG.4729","MSTRG.13544","MSTRG.2924","MSTRG.11309","MSTRG.7125",
"MSTRG.11584","MSTRG.427","MSTRG.5568","MSTRG.5480","MSTRG.6814","MSTRG.11774","MSTRG.8296","MSTRG.1767","MSTRG.4957","MSTRG.3459","MSTRG.10839","MSTRG.6843","MSTRG.12201",
"MSTRG.2497","MSTRG.3349","MSTRG.7893","MSTRG.4978","MSTRG.4504","MSTRG.10469","MSTRG.6271","MSTRG.1954","MSTRG.13837","MSTRG.9163","MSTRG.12787","MSTRG.5317","MSTRG.13191",
"MSTRG.8305","MSTRG.1591","MSTRG.2743","MSTRG.1785","MSTRG.152","MSTRG.11925","MSTRG.6706","MSTRG.9963","MSTRG.10662","MSTRG.973","MSTRG.6632","MSTRG.2998","MSTRG.9006",
"MSTRG.11956","MSTRG.11699","MSTRG.1570","MSTRG.10146","MSTRG.5950","MSTRG.10876","MSTRG.5917","MSTRG.11352","MSTRG.688","MSTRG.3714","MSTRG.6596","MSTRG.5366","MSTRG.12298",
"MSTRG.8006","MSTRG.9128","MSTRG.1358","MSTRG.118","MSTRG.12760","MSTRG.8299","MSTRG.12503","gene-LOC107965675","MSTRG.6947","MSTRG.4986","MSTRG.3892","MSTRG.1984","MSTRG.4756",
"MSTRG.13036","MSTRG.11243","MSTRG.9432","MSTRG.2156","MSTRG.1641","MSTRG.1565","MSTRG.11847","MSTRG.11248","MSTRG.10840","MSTRG.3403","MSTRG.1103","MSTRG.11821","MSTRG.6071",
"MSTRG.9746","MSTRG.3896","MSTRG.5455","MSTRG.4835","MSTRG.1959","MSTRG.3889","MSTRG.8851","MSTRG.546","MSTRG.974","MSTRG.10138","MSTRG.8468","MSTRG.11408","MSTRG.9910",
"MSTRG.7639","MSTRG.1787","MSTRG.3521","MSTRG.6633","MSTRG.3174","MSTRG.261","MSTRG.12487","MSTRG.2136","MSTRG.4674","MSTRG.7956","MSTRG.5282","MSTRG.6940","MSTRG.6859",
"MSTRG.5620","MSTRG.7650","MSTRG.10552","MSTRG.11927","MSTRG.6454","MSTRG.8900","MSTRG.4361","MSTRG.8447","MSTRG.5141","MSTRG.10357","MSTRG.7873","MSTRG.3875","MSTRG.6594",
"MSTRG.11455","MSTRG.13057","MSTRG.7942","MSTRG.13786","MSTRG.11216","MSTRG.6260","MSTRG.7593","MSTRG.3201","MSTRG.1264","MSTRG.7154","MSTRG.5210","MSTRG.11261","MSTRG.11496",
"MSTRG.4866","MSTRG.9200","MSTRG.5385","MSTRG.8737","MSTRG.1097","MSTRG.3562","MSTRG.9866","MSTRG.6293","MSTRG.3339","MSTRG.2718","MSTRG.6517","MSTRG.3013","MSTRG.676",
"MSTRG.5975","MSTRG.11302","MSTRG.8872","MSTRG.5663","MSTRG.11780","MSTRG.9384","MSTRG.6686","MSTRG.3876","MSTRG.5198","MSTRG.7564","MSTRG.5443","MSTRG.12945","MSTRG.10817",
"MSTRG.9379","MSTRG.10196","MSTRG.755","MSTRG.9528","MSTRG.2030","MSTRG.5618","MSTRG.10400","MSTRG.7715","MSTRG.1926","MSTRG.3078","MSTRG.7888","MSTRG.5761","MSTRG.14",
"MSTRG.6835","MSTRG.1756","MSTRG.534","MSTRG.8889","MSTRG.3611","MSTRG.6481","MSTRG.1721","MSTRG.5545","MSTRG.6925","MSTRG.4409","MSTRG.10409","MSTRG.10739","MSTRG.4856",
"MSTRG.9087","MSTRG.7254","MSTRG.4956","MSTRG.6483","MSTRG.12016","MSTRG.3647","MSTRG.13316","MSTRG.10164","MSTRG.4652","MSTRG.6943","MSTRG.9760","MSTRG.7818","MSTRG.9021",
"MSTRG.9356","MSTRG.6324","MSTRG.6902","MSTRG.11665","MSTRG.268","MSTRG.10834","MSTRG.12143","MSTRG.12692","MSTRG.4204","MSTRG.8009","MSTRG.4343","MSTRG.11825","MSTRG.6269",
"MSTRG.6459","MSTRG.2148","MSTRG.2558","MSTRG.12256","MSTRG.10160","MSTRG.643","MSTRG.3162","MSTRG.1276","MSTRG.2071","MSTRG.7959","MSTRG.5718","MSTRG.12940","MSTRG.4324",
"MSTRG.11231","MSTRG.12497","MSTRG.4187","MSTRG.11275","MSTRG.9873","MSTRG.3919","MSTRG.2027","MSTRG.1339","MSTRG.3095","MSTRG.4239","MSTRG.8474","MSTRG.13192","MSTRG.9387",
"MSTRG.12316","MSTRG.13533","MSTRG.12153","MSTRG.11022","MSTRG.687","MSTRG.4964","MSTRG.13307","MSTRG.8916","MSTRG.1546","MSTRG.10917","MSTRG.4758","MSTRG.4235","MSTRG.12155",
"MSTRG.9660","MSTRG.8523","MSTRG.11882","MSTRG.12510","MSTRG.10647","MSTRG.13343","MSTRG.5567","MSTRG.9640","MSTRG.10495","MSTRG.9511","MSTRG.3384","MSTRG.11727","MSTRG.2536",
"MSTRG.6682","MSTRG.13545","MSTRG.9725","MSTRG.11535","MSTRG.11752","MSTRG.5408","MSTRG.2566","MSTRG.4786","MSTRG.2336","MSTRG.12842","MSTRG.8809","MSTRG.11156","MSTRG.530",
"MSTRG.11386","MSTRG.5252","MSTRG.2006")


heatmap.2( assay(Y_rld)[Up_genes_A, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

Down_genes_A <- c("MSTRG.11532","MSTRG.13467","MSTRG.2390","MSTRG.229","MSTRG.7706","MSTRG.5879","MSTRG.4870","MSTRG.6203","MSTRG.10234","MSTRG.4668","MSTRG.6858","MSTRG.3621","MSTRG.9830","MSTRG.3849",
"MSTRG.4660","MSTRG.2337","MSTRG.12640","MSTRG.1388","MSTRG.5694","MSTRG.11777","MSTRG.316","MSTRG.9390","MSTRG.6622","gene-LOC113219228","MSTRG.5712","MSTRG.8659","MSTRG.2447","MSTRG.5765",
"MSTRG.1408","MSTRG.13044","MSTRG.1439","MSTRG.9302","MSTRG.10113","MSTRG.4675","MSTRG.13344","MSTRG.8616","MSTRG.2326","MSTRG.11619","MSTRG.7088","MSTRG.9263","MSTRG.6433","MSTRG.7440",
"MSTRG.5405","MSTRG.2753","MSTRG.1873","MSTRG.9056","MSTRG.10659","MSTRG.7066","MSTRG.1183","MSTRG.7700","MSTRG.13472","MSTRG.9846","MSTRG.2813","MSTRG.8717","MSTRG.5363","MSTRG.6348",
"MSTRG.7292","MSTRG.5812","MSTRG.1559","MSTRG.3156","MSTRG.6639","MSTRG.10748","MSTRG.3493","MSTRG.13669","MSTRG.10936","MSTRG.2830","MSTRG.11505","MSTRG.5237","MSTRG.2210","MSTRG.13652",
"MSTRG.7081","MSTRG.7351","MSTRG.1679","MSTRG.11609","MSTRG.12049","MSTRG.11811","MSTRG.13364","MSTRG.4290","MSTRG.11876","MSTRG.8976","MSTRG.4712","MSTRG.5447","MSTRG.10714","MSTRG.3673",
"MSTRG.786","MSTRG.2243","MSTRG.10280","MSTRG.11480","MSTRG.8091","MSTRG.7273","MSTRG.10958","MSTRG.10670","MSTRG.1136","MSTRG.1147","MSTRG.608","MSTRG.7724","MSTRG.2084","MSTRG.12701",
"MSTRG.3610","MSTRG.12764","MSTRG.8086","MSTRG.12657","MSTRG.12377","MSTRG.9682","MSTRG.12091","MSTRG.4664","MSTRG.12870","MSTRG.1023","MSTRG.4205","MSTRG.3492","MSTRG.7293","MSTRG.10660",
"MSTRG.8543","MSTRG.1843","MSTRG.4193","MSTRG.13492","MSTRG.1635","MSTRG.4422","MSTRG.10360","MSTRG.2135","MSTRG.107","MSTRG.5623","MSTRG.1199","MSTRG.3842","MSTRG.217","MSTRG.6400",
"MSTRG.6133","MSTRG.8270","MSTRG.1162","MSTRG.5375","MSTRG.10708","MSTRG.12564","MSTRG.4731","MSTRG.8050","MSTRG.4782","MSTRG.97","MSTRG.348","MSTRG.13257","MSTRG.4336","MSTRG.6018",
"MSTRG.5225","MSTRG.12868","MSTRG.7404","MSTRG.10209","MSTRG.9249","MSTRG.11865","MSTRG.1295","MSTRG.3419","MSTRG.12389","MSTRG.1594","MSTRG.10352","MSTRG.9001","MSTRG.9142","MSTRG.3194",
"MSTRG.3958","MSTRG.1284","MSTRG.10305","MSTRG.11298","MSTRG.13317","MSTRG.4682","MSTRG.416","MSTRG.7397","MSTRG.10456","MSTRG.4928","MSTRG.3543","MSTRG.1258","MSTRG.13367","MSTRG.9301",
"MSTRG.4281","MSTRG.6822","MSTRG.6416","MSTRG.988","MSTRG.13617","MSTRG.8788","MSTRG.4430","MSTRG.1440","MSTRG.7657","MSTRG.4213","MSTRG.868","MSTRG.13917","MSTRG.13360","MSTRG.7215",
"MSTRG.7705","MSTRG.11611","MSTRG.3675","MSTRG.4130","MSTRG.2120","MSTRG.9218","MSTRG.11782","MSTRG.5868","MSTRG.7767","MSTRG.2711","MSTRG.11977","MSTRG.13475","MSTRG.11906","MSTRG.915",
"MSTRG.13671","MSTRG.3815","MSTRG.13624","MSTRG.8472","MSTRG.3416","MSTRG.7391","MSTRG.5079","MSTRG.13298","MSTRG.6334","MSTRG.10453","MSTRG.4781","MSTRG.1328","MSTRG.4906","MSTRG.2945",
"MSTRG.4564","MSTRG.13668","MSTRG.7780","MSTRG.7109","MSTRG.10946","MSTRG.12475","MSTRG.6604","MSTRG.9800","MSTRG.10290","MSTRG.12093","MSTRG.9370","MSTRG.2403","MSTRG.658","MSTRG.5072",
"MSTRG.10904","MSTRG.1872","MSTRG.7504","MSTRG.5724","MSTRG.12722","MSTRG.10541","MSTRG.2700","MSTRG.8952","MSTRG.9271","MSTRG.1736","MSTRG.12357","MSTRG.10707","MSTRG.12958","MSTRG.12817",
"MSTRG.6718","MSTRG.5393","MSTRG.6134","MSTRG.5920","MSTRG.4665","MSTRG.2362","gene-LOC724515","MSTRG.11805","MSTRG.1875","MSTRG.12481","MSTRG.173","MSTRG.8434","MSTRG.3109","MSTRG.5522",
"MSTRG.11368","MSTRG.1393","MSTRG.3706","MSTRG.12134","MSTRG.3778","MSTRG.10488","MSTRG.7989","MSTRG.1469","MSTRG.3372","MSTRG.9265","MSTRG.13522","MSTRG.7399","MSTRG.4598","MSTRG.1381",
"MSTRG.4832","MSTRG.1560","MSTRG.3323","MSTRG.1646","MSTRG.2214","gene-LOC726193","MSTRG.848","MSTRG.1379","MSTRG.9727","MSTRG.12081","MSTRG.10571","MSTRG.10959","MSTRG.114","MSTRG.4158",
"MSTRG.4930","MSTRG.5450","MSTRG.1176","MSTRG.2449","MSTRG.11896","MSTRG.7065","MSTRG.2638","MSTRG.10990","MSTRG.8803","MSTRG.3634","MSTRG.2621","MSTRG.11663","MSTRG.2509","MSTRG.5083",
"MSTRG.12589","MSTRG.6577","MSTRG.11586","MSTRG.564","MSTRG.9262","MSTRG.11671","MSTRG.6521","MSTRG.4433","MSTRG.1032","MSTRG.11621","MSTRG.7846","MSTRG.11915","MSTRG.11892","MSTRG.9680",
"MSTRG.9566","MSTRG.3729","MSTRG.565","MSTRG.10704","MSTRG.11458","MSTRG.7770","MSTRG.2659","MSTRG.10083","MSTRG.2848","MSTRG.3739","MSTRG.13665","MSTRG.7582","MSTRG.8453","MSTRG.10813",
"MSTRG.6194","MSTRG.13189","MSTRG.10661","MSTRG.6675","MSTRG.10891","MSTRG.9488","MSTRG.64","MSTRG.10167","MSTRG.1141","MSTRG.10960","MSTRG.9259","MSTRG.4415","MSTRG.7163","MSTRG.6176",
"MSTRG.8233","MSTRG.6533","MSTRG.6637","gene-LOC113219051","MSTRG.11795","MSTRG.5649","MSTRG.1770","MSTRG.6098","MSTRG.2917","MSTRG.4969","MSTRG.1871","MSTRG.8061","MSTRG.6441",
"MSTRG.12350","MSTRG.6903","MSTRG.9283","MSTRG.9284","MSTRG.397","MSTRG.12562","MSTRG.11488","MSTRG.3623","MSTRG.12666","MSTRG.8931","MSTRG.6047","gene-LOC100576122","MSTRG.8738",
"MSTRG.1525","gene-LOC102655179","MSTRG.4806","MSTRG.10306","MSTRG.10598","MSTRG.3812","MSTRG.8034","MSTRG.3593","MSTRG.417","MSTRG.13212","MSTRG.1636","MSTRG.4549","MSTRG.6081",
"MSTRG.9261","MSTRG.1718","MSTRG.9344","MSTRG.9888","MSTRG.7091","MSTRG.12593","MSTRG.8790","MSTRG.5505","MSTRG.3092","MSTRG.12566","MSTRG.8998","MSTRG.9170","MSTRG.4342","MSTRG.6412",
"MSTRG.9368","MSTRG.1082","MSTRG.5850","MSTRG.845","MSTRG.11893","MSTRG.6456","MSTRG.8403","MSTRG.6678","MSTRG.1748","MSTRG.2379","MSTRG.1724","MSTRG.5778","MSTRG.10654","MSTRG.12236",
"MSTRG.1878","MSTRG.4633","MSTRG.12539","MSTRG.9858","MSTRG.9247","MSTRG.1394","MSTRG.5122","MSTRG.1282","MSTRG.4982","MSTRG.1874","MSTRG.12753","MSTRG.5772","MSTRG.2241","MSTRG.12857",
"MSTRG.9918","MSTRG.4602","MSTRG.8109","MSTRG.8826","MSTRG.6618","MSTRG.9400","MSTRG.3123","MSTRG.12085","MSTRG.3387","MSTRG.3294","MSTRG.940","MSTRG.7610","MSTRG.12315","MSTRG.9232",
"MSTRG.4188","MSTRG.12375","MSTRG.3065","MSTRG.11074","MSTRG.5313","MSTRG.10216","MSTRG.6330","MSTRG.12834","MSTRG.9106","MSTRG.8612","MSTRG.5106","MSTRG.13281","MSTRG.13209","MSTRG.2065",
"MSTRG.1155","MSTRG.10200","MSTRG.5840","MSTRG.574","MSTRG.13865","MSTRG.7990","MSTRG.2400","MSTRG.10537","MSTRG.2651","MSTRG.2704","MSTRG.9318","MSTRG.1870","MSTRG.6085","MSTRG.6093",
"MSTRG.10621","MSTRG.3475","MSTRG.2391","MSTRG.3722","MSTRG.10716","MSTRG.11768","MSTRG.9579","MSTRG.5266","MSTRG.3718","MSTRG.1461","MSTRG.960","MSTRG.2947","gene-LOC107965541",
"MSTRG.10713","MSTRG.1224","MSTRG.5854","MSTRG.1752","MSTRG.1207","MSTRG.1407","MSTRG.12595","MSTRG.6414","MSTRG.3033","MSTRG.3749","MSTRG.3428","MSTRG.5684","MSTRG.4057","MSTRG.3155",
"MSTRG.7181","MSTRG.52","MSTRG.9107","MSTRG.8452","MSTRG.1192","MSTRG.1371","MSTRG.10135","MSTRG.8617","MSTRG.7405","MSTRG.7222","MSTRG.5189","MSTRG.8912","MSTRG.2087","MSTRG.1299",
"MSTRG.1021","MSTRG.6841","MSTRG.3371","MSTRG.2055","MSTRG.6729","MSTRG.12052","MSTRG.7089","MSTRG.7458","MSTRG.11877","MSTRG.4320","MSTRG.7423","MSTRG.4663","MSTRG.13340","MSTRG.3909",
"MSTRG.1389","MSTRG.1620","MSTRG.729","MSTRG.8357","MSTRG.12084","MSTRG.2819","MSTRG.3045","MSTRG.8544","MSTRG.11717","MSTRG.4810","gene-LOC725058","MSTRG.11898","MSTRG.5755","MSTRG.11526",
"MSTRG.10910","MSTRG.956","MSTRG.4825","MSTRG.8137","MSTRG.10852","MSTRG.3202","MSTRG.2507","MSTRG.9475","MSTRG.2242","MSTRG.10304","MSTRG.2209","MSTRG.12754","MSTRG.9722","MSTRG.13201",
"MSTRG.11575","MSTRG.911","MSTRG.11137","MSTRG.10204","MSTRG.12420","MSTRG.5849","MSTRG.10291","MSTRG.12556","MSTRG.5359","MSTRG.12835","MSTRG.2091","MSTRG.7176","MSTRG.12472","MSTRG.2462",
"MSTRG.12092","MSTRG.11034","MSTRG.3476","MSTRG.5500","MSTRG.4935","MSTRG.11595","MSTRG.7226","MSTRG.13637","MSTRG.7827","MSTRG.2200","MSTRG.12540","MSTRG.13297","MSTRG.8412","MSTRG.1877",
"MSTRG.2952","MSTRG.2034","MSTRG.7769","MSTRG.4380","MSTRG.4667","MSTRG.277","MSTRG.10303","MSTRG.11632","MSTRG.5267","MSTRG.10957","MSTRG.10236","MSTRG.5390","MSTRG.6106","MSTRG.8064",
"MSTRG.9113","MSTRG.2472","MSTRG.13706")



heatmap.2( assay(Y_rld)[Down_genes_A, ], scale="row",
trace="none", dendrogram="column",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


heatmap.2( assay(Y_rld)[c(Up_genes_L,Up_genes_PP,Up_genes_P,Up_genes_A), ], scale="row",
trace="none", dendrogram="none",
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

length(Up_genes_L)
length(Up_genes_PP)
length(Up_genes_P)
length(Up_genes_A)


Up_genes_all <- rbind(assay(Y_rld)[c(Up_genes_L), ],assay(Y_rld)[c(Up_genes_PP), ],assay(Y_rld)[c(Up_genes_P), ],assay(Y_rld)[c(Up_genes_A), ])

pdf("Upgenes_all.pdf",width=6,height=10)
heatmap( Up_genes_all, scale="row",
keep.dendro = F,  Colv = NA, Rowv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

Down_genes_all <- rbind(assay(Y_rld)[c(Down_genes_L), ],assay(Y_rld)[c(Down_genes_PP), ],assay(Y_rld)[c(Down_genes_P), ],assay(Y_rld)[c(Down_genes_A), ])

pdf("Downgenes_all.pdf",width=6,height=10)
heatmap(Down_genes_all, scale="row",
keep.dendro = F,  Colv = NA, Rowv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()




L <- assay(Y_rld)[c(Up_genes_L), ]
dim(L)
L[1:275,1:12] <- 0
dim(L)
L[1:5,]


PP <- assay(Y_rld)[c(Up_genes_PP), ]
dim(PP)
PP[1:150,1:12] <- 0
dim(PP)
PP[1:5,]

P <- assay(Y_rld)[c(Up_genes_P), ]
dim(P)
P[1:226,1:12] <- 0
dim(P)
P[1:5,]


A <- assay(Y_rld)[c(Up_genes_A), ]
dim(A)
A[1:928,1:12] <- 0
dim(A)
A[1:5,]

pdf("UpgenesL.pdf",width=8,height=8)

heatmap( rbind(assay(Y_rld)[c(Up_genes_L), ], PP, P, A), scale="row",
keep.dendro = F,  Colv = NA, Rowv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

dev.off()


pdf("UpgenesPP.pdf",width=8,height=8)

heatmap( rbind(assay(Y_rld)[c(Up_genes_PP), ], L, P, A), scale="row",
keep.dendro = F,  Colv = NA, Rowv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

dev.off()

pdf("UpgenesP.pdf",width=8,height=8)

heatmap( rbind(assay(Y_rld)[c(Up_genes_P), ], L, PP, A), scale="row",
keep.dendro = F,  Colv = NA, Rowv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

dev.off()

pdf("UpgenesA.pdf",width=8,height=8)

heatmap( rbind(assay(Y_rld)[c(Up_genes_A), ], L, PP, P), scale="row",
keep.dendro = F,  Colv = NA, Rowv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

dev.off()




heatmap.2( assay(Y_rld)[c(Up_genes_PP), ], scale="row",
trace="none",    dendrogram = "none",   
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
heatmap.2( assay(Y_rld)[c(Up_genes_P), ], scale="row",
trace="none",    dendrogram = "none",   
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

heatmap.2( assay(Y_rld)[c(Up_genes_L,Up_genes_PP,Up_genes_P,Up_genes_A), ], scale="row",
trace="none",    dendrogram = "none",  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 1, 0)), 
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))




assay(Y_dds2)["MSTRG.13288",]




Larva_up <-  c("MSTRG.7089","MSTRG.618","MSTRG.1572","MSTRG.8671","MSTRG.1847","MSTRG.6866","MSTRG.7602","MSTRG.1960","MSTRG.9922","MSTRG.12251","MSTRG.13738","MSTRG.11921","MSTRG.1701",
"MSTRG.1507","MSTRG.12776","MSTRG.11767","MSTRG.4259","MSTRG.11386","MSTRG.12232","MSTRG.7904","MSTRG.1523","MSTRG.5192","MSTRG.5666","MSTRG.2975","MSTRG.2831","MSTRG.10631","MSTRG.13687",
"MSTRG.1526","MSTRG.13016","MSTRG.11852","MSTRG.6399","MSTRG.3702","MSTRG.7895","MSTRG.11626","MSTRG.3124","MSTRG.8166","MSTRG.953","MSTRG.3769","MSTRG.310","MSTRG.6563","MSTRG.9416",
"MSTRG.9431","MSTRG.11171","MSTRG.8086","MSTRG.1857","MSTRG.4334","MSTRG.10472","MSTRG.8464","MSTRG.11929","MSTRG.13652","MSTRG.9175","MSTRG.3261","MSTRG.5083","MSTRG.11844","MSTRG.13387",
"MSTRG.8067","MSTRG.252","MSTRG.2796","MSTRG.4325","MSTRG.13863","MSTRG.6246","MSTRG.5087","MSTRG.3773","MSTRG.13047","MSTRG.251","MSTRG.11931","MSTRG.12551","MSTRG.11959","MSTRG.13951",
"MSTRG.1841","MSTRG.7176","MSTRG.12535","MSTRG.12816","MSTRG.13631","MSTRG.5532","MSTRG.13659","MSTRG.11424","MSTRG.6209","MSTRG.2501","MSTRG.393","MSTRG.3602","MSTRG.12880","MSTRG.12383",
"MSTRG.13938","MSTRG.10692","MSTRG.3771","MSTRG.4485","MSTRG.404","MSTRG.8912","MSTRG.253","MSTRG.8525","MSTRG.2685","MSTRG.13084","MSTRG.3753","MSTRG.11926","MSTRG.3188","MSTRG.1629",
"MSTRG.7158","MSTRG.2185","MSTRG.3706","MSTRG.6376","MSTRG.9706","MSTRG.7106","MSTRG.9537","MSTRG.11206","MSTRG.11932","MSTRG.8192","MSTRG.10293","MSTRG.12089","MSTRG.7814","MSTRG.2897",
"MSTRG.4502","MSTRG.5146","MSTRG.5639","MSTRG.2245","MSTRG.10045","MSTRG.7334","MSTRG.7226","MSTRG.3833","MSTRG.5338","MSTRG.8451","MSTRG.13025","MSTRG.11648","MSTRG.9828","MSTRG.10185",
"MSTRG.9883","MSTRG.4834","MSTRG.7655","MSTRG.8811","MSTRG.10195","MSTRG.758","MSTRG.8371","MSTRG.10471","MSTRG.10502","MSTRG.1597","MSTRG.11094","MSTRG.7343","MSTRG.12140","MSTRG.10143",
"MSTRG.7579","MSTRG.1803","MSTRG.10487","MSTRG.2979","MSTRG.12497","MSTRG.1757","MSTRG.4960","MSTRG.11551","MSTRG.7849","MSTRG.7345","MSTRG.7762","MSTRG.10576","MSTRG.3343",
"MSTRG.2763","MSTRG.12523","MSTRG.12026","MSTRG.13562","MSTRG.10652","MSTRG.10628","MSTRG.2212","MSTRG.11695","MSTRG.8082","MSTRG.9038","MSTRG.2633","MSTRG.7583",
"MSTRG.12945","MSTRG.8378","MSTRG.4305","MSTRG.9085","MSTRG.9425","MSTRG.13319","MSTRG.6049","MSTRG.5423","MSTRG.4126","MSTRG.5704","MSTRG.12922","MSTRG.7315","MSTRG.10551","MSTRG.9654",
"MSTRG.11410","MSTRG.5874","MSTRG.5726","MSTRG.3339","MSTRG.9453","MSTRG.10511","MSTRG.13503","MSTRG.3764","MSTRG.944","MSTRG.8553","MSTRG.9661","MSTRG.13295","MSTRG.2841",
"MSTRG.8772","MSTRG.3739","MSTRG.8352","MSTRG.398","MSTRG.9488","MSTRG.10439","MSTRG.3643","MSTRG.12770","MSTRG.925","MSTRG.4744","MSTRG.2976","MSTRG.11311","MSTRG.4815","MSTRG.13020",
"MSTRG.9209","MSTRG.7636","MSTRG.6955","MSTRG.5236","MSTRG.587","MSTRG.9334","MSTRG.1220","MSTRG.4850","MSTRG.2079","MSTRG.10802","MSTRG.1194","MSTRG.7835","MSTRG.7548","MSTRG.7826",
"MSTRG.11863","MSTRG.7709","MSTRG.2746","MSTRG.12272","MSTRG.10729","MSTRG.13315","MSTRG.12367","MSTRG.10339","MSTRG.9313","MSTRG.11500","MSTRG.3637","MSTRG.11779","MSTRG.10196",
"MSTRG.3006","MSTRG.7819","MSTRG.1229","MSTRG.10850","MSTRG.13383")

Prepupa_up <- c("MSTRG.10069","MSTRG.7502","MSTRG.8527","MSTRG.793","MSTRG.53","MSTRG.10991","MSTRG.11823","MSTRG.10651","MSTRG.5386","MSTRG.12659","MSTRG.9364","MSTRG.2551","MSTRG.10366",
"MSTRG.5812","MSTRG.7991","MSTRG.10236","MSTRG.11037","MSTRG.7273","MSTRG.4090","MSTRG.3063","MSTRG.4039","MSTRG.6025","MSTRG.8760","MSTRG.13657","MSTRG.1405","MSTRG.2468","MSTRG.12475",
"MSTRG.10505","MSTRG.9316","MSTRG.10681","MSTRG.10813","MSTRG.11336","MSTRG.12384","MSTRG.8412","MSTRG.4653","MSTRG.1172","MSTRG.10653","MSTRG.13360","MSTRG.3493","MSTRG.9092","MSTRG.12228",
"MSTRG.4336","MSTRG.9519","MSTRG.8976","MSTRG.1295","MSTRG.12090","MSTRG.5057","MSTRG.11298","MSTRG.11911","MSTRG.7710","MSTRG.8945","MSTRG.9318","MSTRG.4590","MSTRG.11379",
"MSTRG.6334","MSTRG.11899","MSTRG.3508","MSTRG.786","MSTRG.11591","MSTRG.4521","MSTRG.5600","MSTRG.5962","MSTRG.4869","MSTRG.4982","MSTRG.10193","MSTRG.5182","MSTRG.816","MSTRG.7851",
"MSTRG.10541","MSTRG.299","MSTRG.11782","MSTRG.2404","MSTRG.11708","MSTRG.6648","MSTRG.6210","MSTRG.4188","MSTRG.12754","MSTRG.5691","MSTRG.10291","MSTRG.12379","MSTRG.11567","MSTRG.44",
"MSTRG.3803","MSTRG.12351","MSTRG.151","MSTRG.5920","MSTRG.3617","MSTRG.13861","MSTRG.6388","MSTRG.8611","MSTRG.7214","MSTRG.13867","MSTRG.12634","MSTRG.12602","MSTRG.1147","MSTRG.2345",
"MSTRG.8149","MSTRG.8521","MSTRG.3137","MSTRG.11717","MSTRG.10352","MSTRG.6996","MSTRG.3230","MSTRG.12885","MSTRG.7418","MSTRG.10581","MSTRG.7541","MSTRG.825","MSTRG.4030","MSTRG.5040",
"MSTRG.6661","MSTRG.3842","MSTRG.6392","MSTRG.10198","MSTRG.12352","MSTRG.2811","MSTRG.3158","MSTRG.2961","MSTRG.1662","MSTRG.1953","MSTRG.5177","MSTRG.807","MSTRG.11683","MSTRG.5375",
"MSTRG.6018","MSTRG.2162","MSTRG.1162","MSTRG.8052","MSTRG.4870","MSTRG.12134","MSTRG.11865","MSTRG.348","MSTRG.13262","MSTRG.4442","MSTRG.13173","MSTRG.10842","MSTRG.3150","MSTRG.3778",
"MSTRG.12017","MSTRG.8544","MSTRG.1746","MSTRG.3739","MSTRG.1328","MSTRG.3659","MSTRG.4027","MSTRG.11163","MSTRG.10053","MSTRG.10335","MSTRG.5122","MSTRG.1177","MSTRG.3428","MSTRG.12357",
"MSTRG.11224","MSTRG.3128","MSTRG.11311","MSTRG.12531","MSTRG.5237","MSTRG.13288","MSTRG.5044","MSTRG.5256","MSTRG.9619","MSTRG.10071","MSTRG.12543","MSTRG.3345","MSTRG.7157","MSTRG.9275",
"MSTRG.2736","MSTRG.10599","MSTRG.2519","MSTRG.10234","MSTRG.10290","MSTRG.5328","MSTRG.3202","MSTRG.8536","MSTRG.13472","MSTRG.12230","MSTRG.2854","MSTRG.5447","MSTRG.5100","MSTRG.867",
"MSTRG.7684","MSTRG.2244","MSTRG.13257","MSTRG.2135","MSTRG.13823","MSTRG.3452","MSTRG.6858","MSTRG.12481","MSTRG.12075","MSTRG.12066","MSTRG.2147","MSTRG.13269","MSTRG.13584","MSTRG.10200",
"MSTRG.7664","MSTRG.12104","MSTRG.7094","MSTRG.4635","MSTRG.5381","MSTRG.12303","MSTRG.9921","MSTRG.11704","MSTRG.3080","MSTRG.5393","MSTRG.1116","MSTRG.3094","MSTRG.1333","MSTRG.2982",
"MSTRG.13518","MSTRG.11609","MSTRG.2196","MSTRG.1805","MSTRG.565","MSTRG.10722","MSTRG.1379","MSTRG.4969","MSTRG.11238","MSTRG.6845","MSTRG.10959","MSTRG.10171","MSTRG.9344","MSTRG.12049",
"MSTRG.8308","MSTRG.11803","MSTRG.10274","MSTRG.3120","MSTRG.6916","MSTRG.3849","MSTRG.10109","MSTRG.2581","MSTRG.9226","MSTRG.3652","MSTRG.6164","MSTRG.10475","MSTRG.2447","MSTRG.1176",
"MSTRG.13591","MSTRG.7711","MSTRG.7977","MSTRG.12140","MSTRG.4224","MSTRG.1293","MSTRG.4549","MSTRG.7845","MSTRG.3145","MSTRG.10280","MSTRG.224","MSTRG.9118","MSTRG.9659","MSTRG.2337")

Pupa_up <- c("MSTRG.9636","MSTRG.12873","MSTRG.4840","MSTRG.3809","MSTRG.12554","MSTRG.3371","MSTRG.2267","MSTRG.9456","MSTRG.11088","MSTRG.2405","MSTRG.1240","MSTRG.10889","MSTRG.9383","MSTRG.1441",
"MSTRG.1364","MSTRG.12776","MSTRG.6393","MSTRG.3218","MSTRG.10214","MSTRG.13619","MSTRG.11286","MSTRG.3401","MSTRG.6569","MSTRG.1498","MSTRG.5598","MSTRG.10997","MSTRG.7971","MSTRG.7755",
"MSTRG.8452","MSTRG.2102","MSTRG.12990","MSTRG.10984","MSTRG.13513","MSTRG.2381","MSTRG.6469","MSTRG.13012","MSTRG.11568","MSTRG.12405","MSTRG.75","MSTRG.8044","MSTRG.12569","MSTRG.8678",
"MSTRG.9038","MSTRG.5515","MSTRG.1139","MSTRG.2056","MSTRG.6357","MSTRG.11929","MSTRG.12116","MSTRG.6829","MSTRG.10733","MSTRG.11630","MSTRG.4214","MSTRG.12984","MSTRG.8592","MSTRG.3993",
"MSTRG.5130","MSTRG.6227","MSTRG.3803","MSTRG.11404","MSTRG.4937","MSTRG.8327","MSTRG.11619","MSTRG.1574","MSTRG.8289","MSTRG.6463","MSTRG.11596","MSTRG.9725","MSTRG.5865","MSTRG.13181",
"MSTRG.13493","MSTRG.3284","MSTRG.8656","MSTRG.1290","MSTRG.12982","MSTRG.12045","MSTRG.7138","MSTRG.9582","MSTRG.3014","MSTRG.2945","MSTRG.8311","MSTRG.9780","MSTRG.1485","MSTRG.2752",
"MSTRG.7542","MSTRG.4752","MSTRG.10562","MSTRG.1415","MSTRG.10007","MSTRG.5366","MSTRG.3955","MSTRG.10965","MSTRG.12581","MSTRG.8091","MSTRG.12563","MSTRG.12929","MSTRG.12987",
"MSTRG.2501","MSTRG.10048","MSTRG.11489","MSTRG.4270","MSTRG.2884","MSTRG.13036","MSTRG.1162","MSTRG.2476","MSTRG.8677","MSTRG.11214","MSTRG.12985","MSTRG.10064","MSTRG.9373","MSTRG.11147",
"MSTRG.9982","MSTRG.11511","MSTRG.253","MSTRG.4677","MSTRG.2417","MSTRG.12017","MSTRG.13008","MSTRG.5076","MSTRG.4731","MSTRG.10170","MSTRG.1328","MSTRG.10555","MSTRG.10533","MSTRG.471",
"MSTRG.10439","MSTRG.2143","MSTRG.4434","MSTRG.11743","MSTRG.5788","MSTRG.10425","MSTRG.3621","MSTRG.2527","MSTRG.718","MSTRG.12969","MSTRG.13620","MSTRG.7215","MSTRG.3345",
"MSTRG.11667","MSTRG.622","MSTRG.13577","MSTRG.10807","MSTRG.2897","MSTRG.7557","MSTRG.12803","MSTRG.11744","MSTRG.5457","MSTRG.10802","MSTRG.7209","MSTRG.7301","MSTRG.12993","MSTRG.7496",
"MSTRG.13337","MSTRG.13400","MSTRG.10354","MSTRG.11123","MSTRG.13203","MSTRG.837","MSTRG.11671","MSTRG.1771","MSTRG.881","MSTRG.9635","MSTRG.9250","MSTRG.8260","MSTRG.3318","MSTRG.2415",
"MSTRG.13852","MSTRG.12683","MSTRG.4211","MSTRG.13722","MSTRG.13357","MSTRG.12456","MSTRG.11500","MSTRG.11803","MSTRG.13407","MSTRG.1812","MSTRG.13205","MSTRG.7323","MSTRG.2641","MSTRG.10196",
"MSTRG.623","MSTRG.1566","MSTRG.11131","MSTRG.13309","MSTRG.10473","MSTRG.10327","MSTRG.13310","MSTRG.2830","MSTRG.7726","#N/A")

Adult_up <- c("MSTRG.4124","MSTRG.7593","MSTRG.9562","MSTRG.1675","MSTRG.11634","MSTRG.137","MSTRG.5545","MSTRG.10469","MSTRG.12288","MSTRG.4204","MSTRG.3125","MSTRG.1966","MSTRG.8800","MSTRG.12497",
"MSTRG.427","MSTRG.8916","MSTRG.10552","MSTRG.1892","MSTRG.8927","MSTRG.12466","MSTRG.12392","MSTRG.2566","MSTRG.6880","MSTRG.5385","MSTRG.12801","MSTRG.10286","MSTRG.11369","MSTRG.9805",
"MSTRG.11780","MSTRG.2775","MSTRG.9748","MSTRG.8872","MSTRG.10747","MSTRG.11660","MSTRG.7394","MSTRG.11231","MSTRG.1960","MSTRG.6873","MSTRG.13122","MSTRG.3403","MSTRG.4948",
"MSTRG.8299","MSTRG.5455","MSTRG.13699","MSTRG.6087","MSTRG.2743","MSTRG.9625","MSTRG.7888","MSTRG.7856","MSTRG.13104","MSTRG.5160","MSTRG.13509","MSTRG.4735","MSTRG.13293","MSTRG.6271",
"MSTRG.12662","MSTRG.6974","MSTRG.6856","MSTRG.13587","MSTRG.7564","MSTRG.4268","MSTRG.8980","MSTRG.4826","MSTRG.9097","MSTRG.5917","MSTRG.1243","MSTRG.10158","MSTRG.4555","MSTRG.5110",
"MSTRG.3174","MSTRG.4693","MSTRG.11386","MSTRG.8633","MSTRG.12691","MSTRG.11705","MSTRG.9200","MSTRG.7254","MSTRG.1209","MSTRG.2767","MSTRG.7945","MSTRG.5664","MSTRG.5443",
"MSTRG.2605","MSTRG.1404","MSTRG.12074","MSTRG.11420","MSTRG.1954","MSTRG.3293","MSTRG.6459","MSTRG.6329","MSTRG.10827","MSTRG.973","MSTRG.194","MSTRG.9660","MSTRG.10400","MSTRG.1965",
"MSTRG.1817","MSTRG.9289","MSTRG.11111","MSTRG.12023","MSTRG.13484","MSTRG.13376","MSTRG.1276","MSTRG.10819","MSTRG.13380","MSTRG.1101","MSTRG.9305","MSTRG.1212","MSTRG.11261","MSTRG.6517",
"MSTRG.2338","MSTRG.9479","MSTRG.10363","MSTRG.10866","MSTRG.4187","MSTRG.6704","MSTRG.10967","MSTRG.9174","MSTRG.10026","MSTRG.5155","MSTRG.7546","MSTRG.11585","MSTRG.11838",
"MSTRG.7619","MSTRG.7639","MSTRG.9384","MSTRG.5057","MSTRG.826","MSTRG.6065","MSTRG.1926","MSTRG.6404","MSTRG.4431","MSTRG.9338","MSTRG.13176","MSTRG.8323","MSTRG.11156","MSTRG.12525",
"MSTRG.3412","MSTRG.12510","MSTRG.11847","MSTRG.13279","MSTRG.1933","MSTRG.11320","MSTRG.11852","MSTRG.11535","MSTRG.11821","MSTRG.8945","MSTRG.6383","MSTRG.7531","MSTRG.4396",
"MSTRG.10962","MSTRG.6978","MSTRG.8556","MSTRG.10642","MSTRG.12808","MSTRG.7121","MSTRG.6630","MSTRG.6632","MSTRG.3715","MSTRG.6068","MSTRG.9611","MSTRG.11118","MSTRG.13191","MSTRG.3067",
"MSTRG.1565","MSTRG.142","MSTRG.10161","MSTRG.9393","MSTRG.9356","MSTRG.6805","MSTRG.12296","MSTRG.8296","MSTRG.6044","MSTRG.6423","MSTRG.11907","MSTRG.3932","MSTRG.4285","MSTRG.6403",
"MSTRG.12286","MSTRG.11302","MSTRG.2600","MSTRG.4497","MSTRG.9510","MSTRG.1569","MSTRG.2300","MSTRG.13078","MSTRG.9038","MSTRG.6746","MSTRG.11375","MSTRG.5251","MSTRG.1339",
"MSTRG.1139","MSTRG.2056","MSTRG.12159","MSTRG.2155","MSTRG.13082","MSTRG.2019","MSTRG.13459","MSTRG.4679","MSTRG.12121","MSTRG.9420","MSTRG.3947","MSTRG.12583","MSTRG.746","MSTRG.4089",
"MSTRG.3444","MSTRG.9234","MSTRG.8250","MSTRG.10398","MSTRG.9863","MSTRG.8303","MSTRG.12418","MSTRG.335","MSTRG.13019","MSTRG.9528","MSTRG.5337","MSTRG.11929","MSTRG.6759","MSTRG.6683",
"MSTRG.10357","MSTRG.9057","MSTRG.7383","MSTRG.6422","MSTRG.2880","MSTRG.10544","MSTRG.5224","MSTRG.1956","MSTRG.7198","MSTRG.6947","MSTRG.10358","MSTRG.12713","MSTRG.4343",
"MSTRG.2432","MSTRG.6817","MSTRG.8887","MSTRG.624","MSTRG.1431","MSTRG.13544","MSTRG.3962","MSTRG.546","MSTRG.6813","MSTRG.11665","MSTRG.5538","MSTRG.3467","MSTRG.1640","MSTRG.2664",
"MSTRG.13745","MSTRG.4729","MSTRG.9594","MSTRG.1631","MSTRG.2444","MSTRG.12945","MSTRG.4978","MSTRG.9679","MSTRG.3993","MSTRG.11553","MSTRG.2406","MSTRG.4239","MSTRG.13316",
"MSTRG.2253","MSTRG.4476","MSTRG.13477","MSTRG.10160","MSTRG.8943","MSTRG.4649","MSTRG.13149","MSTRG.8153","MSTRG.5134","MSTRG.1959","MSTRG.2606","MSTRG.13364","MSTRG.6685",
"MSTRG.5619","MSTRG.375","MSTRG.9725","MSTRG.4721","MSTRG.9317","MSTRG.6108","MSTRG.11274","MSTRG.11084","MSTRG.11752","MSTRG.13272","MSTRG.10139","MSTRG.8753","MSTRG.10035",
"MSTRG.1810","MSTRG.13024","MSTRG.8448","MSTRG.3108","MSTRG.3360","MSTRG.11931","MSTRG.5614","MSTRG.5478","MSTRG.6869","MSTRG.12487","MSTRG.2071","MSTRG.5704","MSTRG.8920",
"MSTRG.7893","MSTRG.1747","MSTRG.6402","MSTRG.2900","MSTRG.5212","MSTRG.8244","MSTRG.1244","MSTRG.11091","MSTRG.8523","MSTRG.4736","MSTRG.13226","MSTRG.12636","MSTRG.3963",
"MSTRG.6182","MSTRG.3061","MSTRG.4355","MSTRG.8809","MSTRG.5601","MSTRG.4861","MSTRG.13169","MSTRG.4129","MSTRG.3616","MSTRG.482","MSTRG.13141","MSTRG.10917","MSTRG.2732",
"MSTRG.7418","MSTRG.9717","MSTRG.2998","MSTRG.9357","MSTRG.5378","MSTRG.12348","MSTRG.11882","MSTRG.10409","MSTRG.366","MSTRG.5673","MSTRG.2829","MSTRG.5460","MSTRG.9128",
"MSTRG.3914","MSTRG.1383","MSTRG.8213","MSTRG.9873","MSTRG.5366","MSTRG.9696","MSTRG.5185","MSTRG.4957","MSTRG.8444","MSTRG.10686","MSTRG.9563","MSTRG.226","MSTRG.958","MSTRG.5726",
"MSTRG.4651","MSTRG.2136","MSTRG.10880","MSTRG.10538","MSTRG.7920","MSTRG.3837","MSTRG.3339","MSTRG.6952","MSTRG.13264","MSTRG.5774","MSTRG.7818","MSTRG.1546","MSTRG.4477","MSTRG.8590",
"MSTRG.12650","MSTRG.11424","MSTRG.10831","MSTRG.1602","MSTRG.5153","MSTRG.3577","MSTRG.3602","MSTRG.1130","MSTRG.7785","MSTRG.12316","MSTRG.2336","MSTRG.3197","MSTRG.5099",
"MSTRG.1850","MSTRG.2142","MSTRG.10752","MSTRG.13570","MSTRG.5802","MSTRG.3106","MSTRG.11459","MSTRG.7549","MSTRG.13036","MSTRG.288","MSTRG.13632","MSTRG.4756","MSTRG.11200",
"MSTRG.3173","MSTRG.802","MSTRG.13336","MSTRG.2924","MSTRG.2035","MSTRG.1337","MSTRG.13761","MSTRG.2731","MSTRG.11214","MSTRG.1140","MSTRG.3069","MSTRG.10817","MSTRG.5761","MSTRG.3140",
"MSTRG.10739","MSTRG.11147","MSTRG.13343","MSTRG.6942","MSTRG.12223","MSTRG.7555","MSTRG.1884","MSTRG.4168","MSTRG.13767","MSTRG.233","MSTRG.1280","MSTRG.11701","MSTRG.1755",
"MSTRG.346","MSTRG.12200","MSTRG.13112","MSTRG.3987","MSTRG.1787","MSTRG.5896","MSTRG.9376","MSTRG.1003","MSTRG.9006","MSTRG.2234","MSTRG.13786","MSTRG.3791","MSTRG.12463","MSTRG.9433",
"MSTRG.2323","MSTRG.11243","MSTRG.13065","MSTRG.5956","MSTRG.11926","MSTRG.9059","MSTRG.4777","MSTRG.3077","MSTRG.11275","MSTRG.6369","MSTRG.7398","MSTRG.8367","MSTRG.11954",
"MSTRG.12619","MSTRG.5811","MSTRG.10822","MSTRG.12893","MSTRG.4817","MSTRG.11228","MSTRG.3057","MSTRG.11496","MSTRG.471","MSTRG.3647","MSTRG.3153","MSTRG.10948","MSTRG.13880",
"MSTRG.6200","MSTRG.2881","MSTRG.6850","MSTRG.12843","MSTRG.10495","MSTRG.8900","MSTRG.11798","MSTRG.13533","MSTRG.4744","MSTRG.5208","MSTRG.3641","MSTRG.668","MSTRG.4280",
"MSTRG.13837","MSTRG.11311","MSTRG.306","MSTRG.3203","MSTRG.5019","MSTRG.12161","MSTRG.10425","MSTRG.8601","MSTRG.4815","MSTRG.13207","MSTRG.341","MSTRG.11723","MSTRG.11932","MSTRG.7942",
"MSTRG.1925","MSTRG.3827","MSTRG.2009","MSTRG.3621","MSTRG.2527","MSTRG.11254","MSTRG.5548","MSTRG.10293","MSTRG.1103","MSTRG.9925","MSTRG.2474","MSTRG.2568","MSTRG.9088","MSTRG.1091",
"MSTRG.3345","MSTRG.1126","MSTRG.1483","MSTRG.6955","MSTRG.12621","MSTRG.4361","MSTRG.11445","MSTRG.9021","MSTRG.3354","MSTRG.1331","MSTRG.1984","MSTRG.4347","MSTRG.4186","MSTRG.5390",
"MSTRG.12786","MSTRG.5950","MSTRG.11267","MSTRG.12437","MSTRG.4010","MSTRG.11335","MSTRG.1097","MSTRG.12256","MSTRG.9087","MSTRG.11258","MSTRG.5403","MSTRG.4850","MSTRG.7557",
"MSTRG.1238","MSTRG.8207","MSTRG.7251","MSTRG.1210","MSTRG.7650","MSTRG.9991","MSTRG.9432","MSTRG.12505","MSTRG.2148","MSTRG.6424","MSTRG.6666","MSTRG.6930","MSTRG.9309","MSTRG.10950",
"MSTRG.4608","MSTRG.10802","MSTRG.530","MSTRG.11454","MSTRG.12489","MSTRG.5606","MSTRG.5193","MSTRG.13258","MSTRG.4324","MSTRG.5146","MSTRG.7359","MSTRG.9642","MSTRG.10829",
"MSTRG.11371","MSTRG.1231","MSTRG.7301","MSTRG.11871","MSTRG.6481","MSTRG.11245","MSTRG.152","MSTRG.7275","MSTRG.8474","MSTRG.1288","MSTRG.5481","MSTRG.5387","MSTRG.4166","MSTRG.13276",
"MSTRG.13678","MSTRG.5152","MSTRG.7466","MSTRG.1647","MSTRG.2533","MSTRG.13259","MSTRG.22","MSTRG.2535","MSTRG.11138","MSTRG.2554","MSTRG.12147","MSTRG.1227","MSTRG.11912",
"MSTRG.10867","MSTRG.13960","MSTRG.4489","MSTRG.11643","MSTRG.12462","MSTRG.11863","MSTRG.9557","MSTRG.2311","MSTRG.3998","MSTRG.7758","MSTRG.1830","MSTRG.5565","MSTRG.7154","MSTRG.771",
"MSTRG.5327","MSTRG.10834","MSTRG.9111","MSTRG.6732","MSTRG.1253","MSTRG.10876","MSTRG.7569","MSTRG.1314","MSTRG.4503","MSTRG.4961","MSTRG.13162","MSTRG.3094","MSTRG.2746","MSTRG.3804",
"MSTRG.4758","MSTRG.4086","MSTRG.1880","MSTRG.12574","MSTRG.2744","MSTRG.7731","MSTRG.8174","MSTRG.5663","MSTRG.1821","MSTRG.13524","MSTRG.11927","MSTRG.1451","MSTRG.7092",
"MSTRG.10359","MSTRG.11356","MSTRG.1579","MSTRG.4222","MSTRG.4471","MSTRG.431","MSTRG.13754","MSTRG.1957","MSTRG.5973","MSTRG.4686","MSTRG.5317","MSTRG.182","MSTRG.12738","MSTRG.4303",
"MSTRG.8468","MSTRG.13750","MSTRG.1704","MSTRG.12503","MSTRG.2429","MSTRG.4189","MSTRG.12883","MSTRG.8459","MSTRG.11500","MSTRG.11688","MSTRG.6966","gene-5-ht7","MSTRG.8515",
"MSTRG.11124","MSTRG.7556","MSTRG.3062","MSTRG.5169","MSTRG.12546","MSTRG.5566","MSTRG.12596","MSTRG.13067","MSTRG.5034","MSTRG.5210","MSTRG.242","MSTRG.12047","MSTRG.7601",
"MSTRG.12842","MSTRG.13692","MSTRG.8211","MSTRG.2027","MSTRG.10195","MSTRG.838","MSTRG.6620","MSTRG.2502","MSTRG.11854","MSTRG.4464","MSTRG.2757","MSTRG.3925","MSTRG.1882","MSTRG.12821",
"MSTRG.13469","MSTRG.12356","MSTRG.10196","MSTRG.12032","MSTRG.9497","MSTRG.10769","MSTRG.4289","MSTRG.5292","MSTRG.10115","MSTRG.6164","MSTRG.13307","MSTRG.4786","MSTRG.6835",
"MSTRG.7974","MSTRG.1077","MSTRG.2221","MSTRG.13454","MSTRG.13075","MSTRG.11447","MSTRG.13076","MSTRG.1237","MSTRG.10043","MSTRG.3004","MSTRG.9914","MSTRG.10897","MSTRG.1600",
"MSTRG.1741","MSTRG.12570","MSTRG.8851","MSTRG.9526","MSTRG.4734","MSTRG.11399","MSTRG.8336","MSTRG.9306","MSTRG.8187","MSTRG.2870","MSTRG.10850","MSTRG.8526","MSTRG.1793",
"MSTRG.2156","MSTRG.1131","MSTRG.10435","MSTRG.8865","MSTRG.8006","MSTRG.5659","MSTRG.1803","MSTRG.2616","MSTRG.7019","MSTRG.8549","MSTRG.2486","MSTRG.13045","MSTRG.9422","MSTRG.7129",
"MSTRG.11682","MSTRG.12199")

length(Larva_up)
heatmap( rbind(assay(Y_rld)[c(Larva_up[c(1:106,108:237)]), ]), scale="row",
keep.dendro = F,  Colv = NA, Rowv=NA,
col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

c(Prepupa_up[55:56])
