## Totaba DEG analysis with edgeR
## Felipe Aguilera, June 2020

## Import and pre-processing

# Set working directory
setwd("/Users/faguilera/Documents/")

# Install libraries
install.packages("BiocManager")
BiocManager::install("edgeR")
install.packages("readxl")
install.packages("ggplot2")
install.packages("NMF")

# Load libraries
library(readxl)
library(edgeR)
library(ggplot2)
library(NMF)

# Import data from feature-count Excel file
Totoaba.rawreadcounts <- read_excel("Totoaba_feature-counts.xlsx")

# Get read counts from excel file
Totoaba.readcounts <- Totoaba.rawreadcounts[ ,2:ncol(Totoaba.rawreadcounts)]

# Convert table of read counts to matrix
Totoaba.readcounts <- as.matrix(Totoaba.readcounts)

# Add Gene ID to matrix of read counts
row.names(Totoaba.readcounts) <- Totoaba.rawreadcounts$`Gene-IDs`

# Assign conditions (first two are control, second two are bglucan treatment, third two are inulin treatment, fourth two are quitosan treatment)
Totoaba.sample_conditions <- factor(c(rep("Control", 2), rep("Bglucan", 2), rep("Inulin", 2), rep("Chitosan", 2)))

## Normalization and transformation of read counts

# Create input file for edgeR
Totoaba.edgeRinputdata <- DGEList(counts=Totoaba.readcounts, group=Totoaba.sample_conditions)

# Remove genes without any counts at all samples
Totoaba.edgeRinputdata <- DGEList(counts=Totoaba.readcounts, group=Totoaba.sample_conditions, remove.zeros=TRUE)

# Recompute sequencing library sizes after filtering
Totoaba.edgeRinputdata$samples$lib.size <- colSums(Totoaba.edgeRinputdata$counts)

# Calculate normalization factor for the library sizes - normalization
Totoaba.edgeRinputdata <- calcNormFactors(Totoaba.edgeRinputdata, method="TMM")

# Retrieve normalized read counts
Totoaba.cpmcounts.normalized <- cpm(Totoaba.edgeRinputdata)

# Transform normalized read counts to log2 transformed read counts using a pseudocount of 1
Totoaba.log2norm.counts <- log2(Totoaba.cpmcounts.normalized +1)

# Plot normalized and log2-transformed read counts
pdf("Totoaba.normalized.log2-transformed.read-counts.edgeR.pdf")
par(mfrow=c(2,1))
boxplot(Totoaba.cpmcounts.normalized, nocth=TRUE, main="normalized read counts", ylab="TMM normalized read counts")
boxplot(Totoaba.log2norm.counts, notch=TRUE, main="log2-transformed read counts", ylab="log2(read counts)")
dev.off()

# Plot pairwise comparison of Totaba replicate samples
pdf("Totoaba.pairwise-comparison.control-replicates.edgeR.pdf")
plot(Totoaba.log2norm.counts[,1:2], cex=0.2, main="log-2 transformed read counts", ylab="Control replicate 1 log2(read counts +1)", xlab="Control replicate 2 log2(read counts +1)")
dev.off()
pdf("Totoaba.pairwise-comparison.bglucan-replicates.edgeR.pdf")
plot(Totoaba.log2norm.counts[,3:4], cex=0.2, main="log-2 transformed read counts", ylab="Bglucan replicate 1 log2(read counts +1)", xlab="Bglucan replicate 2 log2(read counts +1)")
dev.off()
pdf("Totoaba.pairwise-comparison.inulin-replicates.edgeR.pdf")
plot(Totoaba.log2norm.counts[,5:6], cex=0.2, main="log-2 transformed read counts", ylab="Inulin replicate 1 log2(read counts +1)", xlab="Inulin replicate 2 log2(read counts +1)")
dev.off()
pdf("Totoaba.pairwise-comparison.chitosan-replicates.edgeR.pdf")
plot(Totoaba.log2norm.counts[,7:8], cex=0.2, main="log-2 transformed read counts", ylab="Chitosan replicate 1 log2(read counts +1)", xlab="Chitosan replicate 2 log2(read counts +1)")
dev.off()

# Pairwise Pearson correlation between Totoaba samples
Totoaba.pearson.correlation.edgeR <- cor(Totoaba.log2norm.counts)
write.table(Totoaba.pearson.correlation.edgeR, file="Totoaba.pairwise-Pearson-correlation.edgeR.txt", quote=FALSE, sep="\t", col.names=NA)

# Hierarchical clustering of Totoaba samples using complete linkage function
Totoaba.distance.mlog2 <- as.dist(1 - cor(Totoaba.log2norm.counts, method="pearson"))
pdf("Totoaba.hierarchical-clustering-samples.edgeR.pdf")
plot(hclust(Totoaba.distance.mlog2), labels=colnames(Totoaba.distance.mlog2), main="Hierarchical clustering/Dendogram")
dev.off()

# Principal Components Analysis of Totoaba samples
pdf("Totoaba.PCA-samples.edgeR.pdf")
plotMDS(Totoaba.edgeRinputdata, method="bcv", col=as.numeric(Totoaba.edgeRinputdata$samples$group))
dev.off()

## Differential gene expression analysis with edgeR

# Set control group as reference level-factor - control
Totoaba.sample_conditions <- relevel(Totoaba.sample_conditions, ref="Control")

# Run edgeR pipeline with four treatment (control group vs three treatments) and processing results
Totoaba.design <- model.matrix(object=~Totoaba.sample_conditions)
Totoaba.edgeRinputdata <- estimateDisp(Totoaba.edgeRinputdata, Totoaba.design)
Totoaba.edgeR.glmFit <- glmFit(Totoaba.edgeRinputdata, Totoaba.design)
Totoaba.edgeR.gmlLRT <- glmLRT(Totoaba.edgeR.glmFit, coef=2:4)
Totoaba.edgeR.gmlLRT.FDR <- topTags(Totoaba.edgeR.gmlLRT, n=Inf, sort.by="none", adjust.method="BH")
Totoaba.edgeR.gmlLRT.FDR$table$thresholdBglucan <- as.logical(Totoaba.edgeR.gmlLRT.FDR$table$FDR < 0.05 & abs(Totoaba.edgeR.gmlLRT.FDR$table$logFC.Totoaba.sample_conditionsBglucan) > 2.0)
Totoaba.edgeR.gmlLRT.FDR$table$thresholdInulin <- as.logical(Totoaba.edgeR.gmlLRT.FDR$table$FDR < 0.05 & abs(Totoaba.edgeR.gmlLRT.FDR$table$logFC.Totoaba.sample_conditionsInulin) > 2.0)
Totoaba.edgeR.gmlLRT.FDR$table$thresholdChitosan <- as.logical(Totoaba.edgeR.gmlLRT.FDR$table$FDR < 0.05 & abs(Totoaba.edgeR.gmlLRT.FDR$table$logFC.Totoaba.sample_conditionsChitosan) > 2.0)
write.table(Totoaba.edgeR.gmlLRT.FDR$table, file="Totoaba.edgeR-results.txt", quote=FALSE, sep="\t", col.names=NA)
Totoaba.DEGs.edgeR.Bglucan <- row.names(Totoaba.edgeR.gmlLRT.FDR$table)[which(Totoaba.edgeR.gmlLRT.FDR$table$thresholdBglucan)]
Totoaba.DEGs.edgeR.Inulin <- row.names(Totoaba.edgeR.gmlLRT.FDR$table)[which(Totoaba.edgeR.gmlLRT.FDR$table$thresholdInulin)]
Totoaba.DEGs.edgeR.Chitosan <- row.names(Totoaba.edgeR.gmlLRT.FDR$table)[which(Totoaba.edgeR.gmlLRT.FDR$table$thresholdChitosan)]
write.table(Totoaba.DEGs.edgeR.Bglucan, file="Totoaba.DEGs.edgeR.Bglucan.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(Totoaba.DEGs.edgeR.Inulin, file="Totoaba.DEGs.edgeR.Inulin.txt", quote=FALSE, sep="\t", col.names=NA)
write.table(Totoaba.DEGs.edgeR.Chitosan, file="Totoaba.DEGs.edgeR.Chitosan.txt", quote=FALSE, sep="\t", col.names=NA)

# MA plot of Totoaba between control group and treatments in a independent manner
Totoaba.edgeR.MAplot <- as.data.frame(Totoaba.edgeR.gmlLRT.FDR)
pdf("Totoaba.MA-plot.edgeR.Bglucan.pdf")
ggplot(data=Totoaba.edgeR.MAplot, aes(x=logCPM, y=logFC.Totoaba.sample_conditionsBglucan, color=thresholdBglucan)) + geom_point(size=0.25) + scale_color_manual(values=c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept=2, linetype="dashed", color="blue", size=0.5) + geom_hline(yintercept=-2, linetype="dashed", color="blue", size=0.5) + guides(color=guide_legend("Significance"))
dev.off()
pdf("Totoaba.MA-plot.edgeR.Inulin.pdf")
ggplot(data=Totoaba.edgeR.MAplot, aes(x=logCPM, y=logFC.Totoaba.sample_conditionsInulin, color=thresholdInulin)) + geom_point(size=0.25) + scale_color_manual(values=c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept=2, linetype="dashed", color="blue", size=0.5) + geom_hline(yintercept=-2, linetype="dashed", color="blue", size=0.5) + guides(color=guide_legend("Significance"))
dev.off()
pdf("Totoaba.MA-plot.edgeR.Chitosan.pdf")
ggplot(data=Totoaba.edgeR.MAplot, aes(x=logCPM, y=logFC.Totoaba.sample_conditionsChitosan, color=thresholdChitosan)) + geom_point(size=0.25) + scale_color_manual(values=c("black", "red")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(yintercept=2, linetype="dashed", color="blue", size=0.5) + geom_hline(yintercept=-2, linetype="dashed", color="blue", size=0.5) + guides(color=guide_legend("Significance"))
dev.off()

# Heatmap of Totoaba DEGs
Totoaba.edgeR.heatmap <- Totoaba.edgeR.gmlLRT.FDR$table[order(Totoaba.edgeR.gmlLRT.FDR$table$FDR), ]
Totoaba.DEGenes.Bglucan <- rownames(subset(Totoaba.edgeR.heatmap, FDR < 0.05 & abs(Totoaba.edgeR.heatmap$logFC.Totoaba.sample_conditionsBglucan) > 2.0))
Totoaba.DEGenes.Inulin <- rownames(subset(Totoaba.edgeR.heatmap, FDR < 0.05 & abs(Totoaba.edgeR.heatmap$logFC.Totoaba.sample_conditionsInulin) > 2.0))
Totoaba.DEGenes.Chitosan <- rownames(subset(Totoaba.edgeR.heatmap, FDR < 0.05 & abs(Totoaba.edgeR.heatmap$logFC.Totoaba.sample_conditionsChitosan) > 2.0))
Totoaba.DEGenes.Bglucan <- as.data.frame(Totoaba.DEGenes.Bglucan)
Totoaba.DEGenes.Inulin <- as.data.frame(Totoaba.DEGenes.Inulin)
Totoaba.DEGenes.Chitosan <- as.data.frame(Totoaba.DEGenes.Chitosan)
colnames(Totoaba.DEGenes.Bglucan) <- c("Gene-IDs")
colnames(Totoaba.DEGenes.Inulin) <- c("Gene-IDs")
colnames(Totoaba.DEGenes.Chitosan) <- c("Gene-IDs")
Totoaba.DEGenes.heatmap.IDs <- rbind(Totoaba.DEGenes.Bglucan, Totoaba.DEGenes.Inulin, Totoaba.DEGenes.Chitosan)
Totoaba.DEGenes.heatmap.IDs <- Totoaba.DEGenes.heatmap.IDs [!duplicated(Totoaba.DEGenes.heatmap.IDs[c(1)]),]
Totoaba.DEGenes.input.heatmap <- Totoaba.log2norm.counts[Totoaba.DEGenes.heatmap.IDs, ]
pdf("Totoaba.heatmap-DEGs.edgeR.alltreatments.pdf")
aheatmap(Totoaba.DEGenes.input.heatmap, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()

# Heatmap of Totoaba DEGs individual treatments
Totoaba.DEGenes.heatmap.Bglucan <- Totoaba.log2norm.counts[Totoaba.DEGs.edgeR.Bglucan, ]
pdf("Totoaba.heatmap-DEGs.edgeR.Bglucan.pdf")
aheatmap(Totoaba.DEGenes.heatmap.Bglucan, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()
Totoaba.DEGenes.heatmap.Inulin <- Totoaba.log2norm.counts[Totoaba.DEGs.edgeR.Inulin, ]
pdf("Totoaba.heatmap-DEGs.edgeR.Inulin.pdf")
aheatmap(Totoaba.DEGenes.heatmap.Inulin, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()
Totoaba.DEGenes.heatmap.Chitosan <- Totoaba.log2norm.counts[Totoaba.DEGs.edgeR.Chitosan, ]
pdf("Totoaba.heatmap-DEGs.edgeR.Chitosan.pdf")
aheatmap(Totoaba.DEGenes.heatmap.Chitosan, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale="row")
dev.off()
## END OF THE SCRIPT!