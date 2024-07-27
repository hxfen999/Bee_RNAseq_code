library(DESeq2)
library(edgeR)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(amap)

setwd("C:\\hxfen\\25.BeeTranscriptome\\github_hxfen")  #set the diroctery containing the expression data as working path. 

#Gene <- read.csv("gene_count_matrix.csv", header=T)
Gene <- read.csv("gene_count.txt", header=T)

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



######For the 2nd day prepupa############
M_Gene.count <- Gene.count[,c("PP2E1","PP2E2","PP2E3","PP2C1","PP2C2","PP2C3")]
M_ID.info <- ID.info[c(9:11,22:24),]
M_ID.info
dim(M_Gene.count)
head(M_Gene.count)
M_dds <- DESeqDataSetFromMatrix(countData = as.matrix(M_Gene.count), colData = M_ID.info, design =~ condition)
M_dds2 <- DESeq(M_dds)
M_rld <- rlog( M_dds2, blind=FALSE )
M_res <- results(M_dds2, pAdjustMethod = "BH", contrast = c("condition", "Editted", "Control"))  
M_resSig <- subset(M_res, padj < 0.05)
M_resSig
M_resSig_lfc1 <- subset(M_resSig, abs(log2FoldChange) > 1 )
M_resSig_lfc1
M_resSig_order <- M_resSig[order(M_resSig$log2FoldChange),]
M_resSig_orderlessthan0 <- M_resSig_order[M_resSig$log2FoldChange<0,]
dim(M_resSig_orderlessthan0)
M_resSig_names <- row.names(M_resSig_order)

M_rlogMat <- assay(M_rld)[M_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
M_pearson_cor <- as.matrix(cor(M_rlogMat, method="pearson"))
M_hc <- hcluster(t(M_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_M_Pearsoncor.pdf",width=7,height=7)
heatmap.2(M_pearson_cor, Rowv=as.dendrogram(M_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- M_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #setting FDR
cut_off_logFC = 1      #setting fold change

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
  # Data, point, color
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
  # Legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p



######For the 3rd day pupa############
P_Gene.count <- Gene.count[,c("P3E1","P3E2","P3E3","P3E4","P3C1","P3C2","P3C3")]
P_ID.info <- ID.info[c(5:8,19:21),]
P_ID.info
dim(P_Gene.count)
head(P_Gene.count)
P_dds <- DESeqDataSetFromMatrix(countData = as.matrix(P_Gene.count), colData = P_ID.info, design =~ condition)
P_dds2 <- DESeq(P_dds)
P_rld <- rlog( P_dds2, blind=FALSE )
P_res <- results(P_dds2, pAdjustMethod = "BH", contrast = c("condition", "Editted", "Control"))  
P_resSig <- subset(P_res, padj < 0.05)
P_resSig
P_resSig_lfc1 <- subset(P_resSig, abs(log2FoldChange) > 1 )
P_resSig_lfc1
P_resSig_order <- P_resSig[order(P_resSig$log2FoldChange),]
P_resSig_orderlessthan0 <- P_resSig_order[P_resSig$log2FoldChange<0,]
dim(P_resSig_orderlessthan0)
P_resSig_names <- row.names(P_resSig_order)

P_rlogMat <- assay(P_rld)[P_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
P_pearson_cor <- as.matrix(cor(P_rlogMat, method="pearson"))
P_hc <- hcluster(t(P_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_P_Pearsoncor.pdf",width=7,height=7)
heatmap.2(P_pearson_cor, Rowv=as.dendrogram(P_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- P_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #setting FDR
cut_off_logFC = 1      #setting fold change

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
  # Data, point, color
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
  # Legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p




######For the 1st day adult############
A_Gene.count <- Gene.count[,c("A1E1","A1E2","A1E3","A1E4","A1C1","A1C2","A1C3")]
A_ID.info <- ID.info[c(1:4,16:18),]
A_ID.info
dim(A_Gene.count)
head(A_Gene.count)
A_dds <- DESeqDataSetFromMatrix(countData = as.matrix(A_Gene.count), colData = A_ID.info, design =~ condition)
A_dds2 <- DESeq(A_dds)
A_rld <- rlog( A_dds2, blind=FALSE )
A_res <- results(A_dds2, pAdjustMethod = "BH", contrast = c("condition", "Editted", "Control"))  
A_resSig <- subset(A_res, padj < 0.05)
A_resSig
A_resSig_lfc1 <- subset(A_resSig, abs(log2FoldChange) > 1 )
A_resSig_lfc1
A_resSig_order <- A_resSig[order(A_resSig$log2FoldChange),]
A_resSig_orderlessthan0 <- A_resSig_order[A_resSig$log2FoldChange<0,]
dim(A_resSig_orderlessthan0)
A_resSig_names <- row.names(A_resSig_order)

A_rlogMat <- assay(A_rld)[A_resSig_names, ]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# calculating pearson correlation
A_pearson_cor <- as.matrix(cor(A_rlogMat, method="pearson"))
A_hc <- hcluster(t(A_rlogMat), method="pearson") 

#pdf("Heatmap_ResSig_A_Pearsoncor.pdf",width=7,height=7)
heatmap.2(A_pearson_cor, Rowv=as.dendrogram(A_hc), symm=T, trace="none", density.info="none",
col=hmcol, margins=c(5,5), main="The Pearson correlation 
of each 6-day larval sample")
#dev.off()

dataset <- A_res
head(dataset)
dataset <- as.data.frame(dataset)
resSig <- subset(dataset, padj < 0.05)
resSig[order(-log(resSig$padj)),]
regSig.dataframe <- as.data.frame(resSig[order(-log(resSig$padj)),])
tail(regSig.dataframe)

cut_off_pvalue = 0.05  #setting FDR
cut_off_logFC = 1      #setting fold change

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
  # Data, point, color
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
  # Legend
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text = element_text(size = 20))
p

