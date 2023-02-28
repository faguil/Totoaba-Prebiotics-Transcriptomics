## Totaba DEG analysis with DESeq2
## Felipe Aguilera, June 2020

## Import and pre-processing

# Set working directory
setwd("/Users/faguilera/Documents/")

# Install libraries
install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("readxl")
BiocManager::install("EnhancedVolcano")
install.packages("NMF")

# Load libraries
library(readxl)
library(DESeq2)
library(EnhancedVolcano)
library(NMF)

# Import data from feature-count Excel file
Totoaba.rawreadcounts <- read_excel("Totoaba_feature-counts.xlsx")

# Get read counts from excel file
Totoaba.readcounts <- Totoaba.rawreadcounts[ ,2:ncol(Totoaba.rawreadcounts)]

# Convert table of read counts to matrix
Totoaba.readcounts <- as.matrix(Totoaba.readcounts)

# Add Gene IDs to matrix of read counts
row.names(Totoaba.readcounts) <- Totoaba.rawreadcounts$`Gene-IDs`

# Assign conditions (first two are control, second two are bglucan treatment, third two are inulin treatment, fourth two are quitosan treatment)
Totoaba.sample_conditions <- factor(c(rep("Control", 2), rep("Bglucan", 2), rep("Inulin", 2), rep("Chitosan", 2)))

# Set control group as reference level-factor - control
Totoaba.sample_conditions <- relevel(Totoaba.sample_conditions, ref="Control")

# Create a colData frame for DESeq2
Totoaba.colData <- data.frame(row.names=colnames(Totoaba.readcounts), Totoaba.sample_conditions)

## Normalization and transformation of read counts

# Load input data for DESeq2
Totoaba.DESeq2inputdata <- DESeqDataSetFromMatrix(countData=Totoaba.readcounts, colData=Totoaba.colData, design=~Totoaba.sample_conditions)

# Remove genes without any counts
Totoaba.DESeq2inputdata <- Totoaba.DESeq2inputdata[ rowSums(counts(Totoaba.DESeq2inputdata)) > 0, ]

# Calculate the size factor and add it to the dataset - normalization
Totoaba.DESeq2inputdata <- estimateSizeFactors(Totoaba.DESeq2inputdata)

# Retrieve normalized read counts
Totoaba.deseq2counts.normalized <- counts(Totoaba.DESeq2inputdata, normalized=TRUE)

# Transform normalized read counts to log2 transformed read counts using a pseudocount of 1
Totoaba.log2norm.counts <- log2(Totoaba.deseq2counts.normalized +1)

# Plot normalized and log2-transformed read counts
pdf("Totoaba.normalized.log2-transformed.read-counts.DESeq2.pdf")
par(mfrow=c(2,1))
boxplot(Totoaba.deseq2counts.normalized, notch=TRUE, main="normalized read counts", ylab="size-factor normalized read counts")
boxplot(Totoaba.log2norm.counts, notch=TRUE, main="log2-transformed read counts", ylab="log2(read counts)")
dev.off()

# Plot pairwise comparison of Totoaba replicate samples
pdf("Totoaba.pairwise-comparison.control-replicates.DESEq2.pdf")
plot(Totoaba.log2norm.counts[,1:2], cex=0.2, main="log2-transformed read counts", ylab="Control replicate 1 log2(read counts +1)", xlab="Control replicate 2 log2(read counts +1)")
dev.off()
pdf("Totoaba.pairwise-comparison.bglucan-replicates.DESEq2.pdf")
plot(Totoaba.log2norm.counts[,3:4], cex=0.2, main="log2-transformed read counts", ylab="Bglucan replicate 1 log2(read counts +1)", xlab="Bglucan replicate 2 log2(read counts +1)")
dev.off()
pdf("Totoaba.pairwise-comparison.inulin-replicates.DESEq2.pdf")
plot(Totoaba.log2norm.counts[,5:6], cex=0.2, main="log2-transformed read counts", ylab="Inulin replicate 1 log2(read counts +1)", xlab="Inulin replicate 2 log2(read counts +1)")
dev.off()
pdf("Totoaba.pairwise-comparison.chitosan-replicates.DESEq2.pdf")
plot(Totoaba.log2norm.counts[,7:8], cex=0.2, main="log2-transformed read counts", ylab="Chitosan replicate 1 log2(read counts +1)", xlab="Chitosan replicate 2 log2(read counts +1)")
dev.off()

# Pairwise Pearson correlation between Totoaba samples
Totoaba.pearson.correlation.DESeq2 <- cor(Totoaba.log2norm.counts)
write.table(Totoaba.pearson.correlation.DESeq2, file="Totoaba.pairwise-Pearson-correlation.DESeq2.txt", quote=FALSE, sep="\t", col.names=NA)

# Transform input data for DESeq2 to regularized log-transformed values for clustering and PCA analyses
Totoaba.DESeq2inputdata.rlog <- rlogTransformation(Totoaba.DESeq2inputdata, blind=TRUE)
Totoaba.rlog.norm.counts <- assay(Totoaba.DESeq2inputdata.rlog)

# Hierarchical clustering of Hspinulosus samples using complete linkage function
Totoaba.distance.m_rlog <- as.dist(1 - cor(Totoaba.rlog.norm.counts, method="pearson"))
pdf("Totoaba.hierarchical-clustering-samples.DESEq2.pdf")
plot(hclust(Totoaba.distance.m_rlog), labels=colnames(Totoaba.rlog.norm.counts), main="Hierarchical clustering/Dendogram")
dev.off()

# Principal Components Analysis of Totoaba samples
pdf("Totoaba.PCA-samples.DESeq2.pdf")
plotPCA(Totoaba.DESeq2inputdata.rlog, intgroup="Totoaba.sample_conditions")
dev.off()

## Differential gene expression analysis with DESeq2

# Set control group as reference level-factor - control
colData(Totoaba.DESeq2inputdata)$Totoaba.sample_conditions <- relevel(colData(Totoaba.DESeq2inputdata)$Totoaba.sample_conditions, "Control")

# Run DESeq2 pipeline and processing results
Totoaba.DESeq2.pipeline <- DESeq(Totoaba.DESeq2inputdata)
Totoaba.DESeq2.Bglucan <- results(Totoaba.DESeq2.pipeline, contrast=c("Totoaba.sample_conditions", "Control", "Bglucan"), alpha=0.05)
Totoaba.DESeq2.Bglucan$threshold <- as.logical(Totoaba.DESeq2.Bglucan$padj < 0.05 & abs(Totoaba.DESeq2.Bglucan$log2FoldChange) > 2.0)
write.table(Totoaba.DESeq2.Bglucan, file="Totoaba.DESeq2-results.Bglucan.txt", quote=FALSE, sep="\t", col.names=NA)
Totoaba.DESeq2.Inulin <- results(Totoaba.DESeq2.pipeline, contrast=c("Totoaba.sample_conditions", "Control", "Inulin"), alpha=0.05)
Totoaba.DESeq2.Inulin$threshold <- as.logical(Totoaba.DESeq2.Inulin$padj < 0.05 & abs(Totoaba.DESeq2.Inulin$log2FoldChange) > 2.0)
write.table(Totoaba.DESeq2.Inulin, file="Totoaba.DESeq2-results.Inulin.txt", quote=FALSE, sep="\t", col.names=NA)
Totoaba.DESeq2.Chitosan <- results(Totoaba.DESeq2.pipeline, contrast=c("Totoaba.sample_conditions", "Control", "Chitosan"), alpha=0.05)
Totoaba.DESeq2.Chitosan$threshold <- as.logical(Totoaba.DESeq2.Chitosan$padj < 0.05 & abs(Totoaba.DESeq2.Chitosan$log2FoldChange) > 2.0)
write.table(Totoaba.DESeq2.Chitosan, file="Totoaba.DESeq2-results.Chitosan.txt", quote=FALSE, sep="\t", col.names=NA)
Totoaba.DEGs.DESeq2.Bglucan <- row.names(Totoaba.DESeq2.Bglucan)[which(Totoaba.DESeq2.Bglucan$threshold)]
write.table(Totoaba.DEGs.DESeq2.Bglucan, file="Totoaba.DEGs.DESeq2.Bglucan.txt", quote=FALSE, sep="\t", col.names=NA)
Totoaba.DEGs.DESeq2.Inulin <- row.names(Totoaba.DESeq2.Inulin)[which(Totoaba.DESeq2.Inulin$threshold)]
write.table(Totoaba.DEGs.DESeq2.Inulin, file="Totoaba.DEGs.DESeq2.Inulin.txt", quote=FALSE, sep="\t", col.names=NA)
Totoaba.DEGs.DESeq2.Chitosan <- row.names(Totoaba.DESeq2.Chitosan)[which(Totoaba.DESeq2.Chitosan$threshold)]
write.table(Totoaba.DEGs.DESeq2.Chitosan, file="Totoaba.DEGs.DESeq2.Chitosan.txt", quote=FALSE, sep="\t", col.names=NA)

# Volcano plot of Totoaba between control group and treatments in a independent manner
Totoaba.DESeq2.Bglucan.volcanoplot <- as.data.frame(Totoaba.DESeq2.Bglucan)
Totoaba.DESeq2.Inulin.volcanoplot <- as.data.frame(Totoaba.DESeq2.Inulin)
Totoaba.DESeq2.Chitosan.volcanoplot <- as.data.frame(Totoaba.DESeq2.Chitosan)
pdf("Totoaba.Volcano-plot.DESeq2.Bglucan.pdf")
EnhancedVolcano(Totoaba.DESeq2.Bglucan.volcanoplot, x='log2FoldChange', y='padj', lab=rownames(Totoaba.DESeq2.Bglucan.volcanoplot), col=c('black', 'black', 'black', 'red'), FCcutoff=2.0, pCutoff=0.05, pointSize=0.5, cutoffLineCol='blue', gridlines.major=FALSE, gridlines.minor=FALSE, border='full', title='Control v/s Bglucan', legendPosition='none', selectLab="none")
dev.off()
pdf("Totoaba.Volcano-plot.DESeq2.Inulin.pdf")
EnhancedVolcano(Totoaba.DESeq2.Inulin.volcanoplot, x='log2FoldChange', y='padj', lab=rownames(Totoaba.DESeq2.Inulin.volcanoplot), col=c('black', 'black', 'black', 'red'), FCcutoff=2.0, pCutoff=0.05, pointSize=0.5, cutoffLineCol='blue', gridlines.major=FALSE, gridlines.minor=FALSE, border='full', title='Control v/s Inulin', legendPosition='none', selectLab="none")
dev.off()
pdf("Totoaba.Volcano-plot.DESeq2.Chitosan.pdf")
EnhancedVolcano(Totoaba.DESeq2.Chitosan.volcanoplot, x='log2FoldChange', y='padj', lab=rownames(Totoaba.DESeq2.Chitosan.volcanoplot), col=c('black', 'black', 'black', 'red'), FCcutoff=2.0, pCutoff=0.05, pointSize=0.5, cutoffLineCol='blue', gridlines.major=FALSE, gridlines.minor=FALSE, border='full', title='Control v/s Chitosan', legendPosition='none', selectLab="none")
dev.off()

# Heatmap of Totoaba DEGs
Totoaba.DEGenes.Bglucan <- as.data.frame(Totoaba.DEGs.DESeq2.Bglucan)
Totoaba.DEGenes.Inulin <- as.data.frame(Totoaba.DEGs.DESeq2.Inulin)
Totoaba.DEGenes.Chitosan <- as.data.frame(Totoaba.DEGs.DESeq2.Chitosan)
colnames(Totoaba.DEGenes.Bglucan) <- c("Gene-IDs")
colnames(Totoaba.DEGenes.Inulin) <- c("Gene-IDs")
colnames(Totoaba.DEGenes.Chitosan) <- c("Gene-IDs")
Totoaba.DEGenes.heatmap.IDs <- rbind(Totoaba.DEGenes.Bglucan, Totoaba.DEGenes.Inulin, Totoaba.DEGenes.Chitosan)
Totoaba.DEGenes.heatmap.IDs <- Totoaba.DEGenes.heatmap.IDs [!duplicated(Totoaba.DEGenes.heatmap.IDs[c(1)]),]
Totoaba.DEGenes.input.heatmap <- Totoaba.log2norm.counts[Totoaba.DEGenes.heatmap.IDs, ]
pdf("Totoaba.heatmap-DEGs.DESeq2.alltreatments.pdf")
aheatmap(Totoaba.DEGenes.input.heatmap, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()

# Heatmap of Totoaba DEGs individual treatments
Totoaba.DEGenes.heatmap.Bglucan <- Totoaba.log2norm.counts[Totoaba.DEGs.DESeq2.Bglucan, ]
pdf("Totoaba.heatmap-DEGs.DESeq2.Bglucan.pdf")
aheatmap(Totoaba.DEGenes.heatmap.Bglucan, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()
Totoaba.DEGenes.heatmap.Inulin <- Totoaba.log2norm.counts[Totoaba.DEGs.DESeq2.Inulin, ]
pdf("Totoaba.heatmap-DEGs.DESeq2.Inulin.pdf")
aheatmap(Totoaba.DEGenes.heatmap.Inulin, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()
Totoaba.DEGenes.heatmap.Chitosan <- Totoaba.log2norm.counts[Totoaba.DEGs.DESeq2.Chitosan, ]
pdf("Totoaba.heatmap-DEGs.DESeq2.Chitosan.pdf")
aheatmap(Totoaba.DEGenes.heatmap.Chitosan, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()

## END OF THE SCRIPT!