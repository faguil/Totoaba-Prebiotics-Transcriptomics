## Totoaba functional enrichment analysis on DEGs from DESeq2 with topGO
## Felipe Aguilera, July 2020

## Import data and setting environments

# Set working directory
setwd("/Users/faguilera/Documents/")

# Install libraries
install.packages("BiocManager")
BiocManager::install("topGO")

# Load libraries
library(topGO)

# FUNCTIONAL GO ENRICHMENT ANALYSIS FOR BGLUCAN TREATMENT #

# Import GO annotations for the genes
Totoaba.GOterms <- readMappings(file="Tmacdonaldi.GO-annotations.txt")

# Define list of genes of interest
Totoaba.GeneIDs <- names(Totoaba.GOterms)
Totoaba.DEGs_Bglucan <- read.table("Totoaba.DEGs.DESeq2.Bglucan.txt", header=FALSE)
Totoaba.DEGs_Bglucan <- as.character(Totoaba.DEGs_Bglucan$V1)
DEGs.Bglucan <- factor(as.integer(Totoaba.GeneIDs%in%Totoaba.DEGs_Bglucan))
names(DEGs.Bglucan) <- Totoaba.GeneIDs

# Put data together into an R object #
Totoaba.GO_BP.Bglucan <- new("topGOdata", description="Totoaba-Bglucan (BP)", ontology="BP", allGenes=DEGs.Bglucan, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_MF.Bglucan <- new("topGOdata", description="Totoaba-Bglucan (MF)", ontology="MF", allGenes=DEGs.Bglucan, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_CC.Bglucan <- new("topGOdata", description="Totoaba-Bglucan (CC)", ontology="CC", allGenes=DEGs.Bglucan, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)

# Perform of GO enrichment test #
Totoaba.Bglucan.BP.Fisher <- runTest(Totoaba.GO_BP.Bglucan, algorithm="weight01", statistic="fisher")
Totoaba.Bglucan.MF.Fisher <- runTest(Totoaba.GO_MF.Bglucan, algorithm="weight01", statistic="fisher")
Totoaba.Bglucan.CC.Fisher <- runTest(Totoaba.GO_CC.Bglucan, algorithm="weight01", statistic="fisher")
Totoaba.Bglucan.allGO.BP <- usedGO(Totoaba.GO_BP.Bglucan)
Totoaba.Bglucan.topGO.BP.results <- GenTable(Totoaba.GO_BP.Bglucan, weightFisher=Totoaba.Bglucan.BP.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Bglucan.allGO.BP))
names(Totoaba.Bglucan.topGO.BP.results)[names(Totoaba.Bglucan.topGO.BP.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.Bglucan.allGO.MF <- usedGO(Totoaba.GO_MF.Bglucan)
Totoaba.Bglucan.topGO.MF.results <- GenTable(Totoaba.GO_MF.Bglucan, weightFisher=Totoaba.Bglucan.MF.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Bglucan.allGO.MF))
names(Totoaba.Bglucan.topGO.MF.results)[names(Totoaba.Bglucan.topGO.MF.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.Bglucan.allGO.CC <- usedGO(Totoaba.GO_CC.Bglucan)
Totoaba.Bglucan.topGO.CC.results <- GenTable(Totoaba.GO_CC.Bglucan, weightFisher=Totoaba.Bglucan.CC.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Bglucan.allGO.CC))
names(Totoaba.Bglucan.topGO.CC.results)[names(Totoaba.Bglucan.topGO.CC.results) == "weightFisher"] <- "weightFisher.pvalue"

# Correct for multiple testing - adjusted p-values #
Totoaba.Bglucan.topGO.BP.padj <- round(p.adjust(Totoaba.Bglucan.topGO.BP.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.BP.Bglucan <- cbind(Totoaba.Bglucan.topGO.BP.results,Totoaba.Bglucan.topGO.BP.padj)
Totoaba.GOenrichment.BP.Bglucan <- Totoaba.GOenrichment.BP.Bglucan[order(Totoaba.GOenrichment.BP.Bglucan$Totoaba.Bglucan.topGO.BP.padj),]
Totoaba.GOenrichment.BP.Bglucan.BH <- Totoaba.GOenrichment.BP.Bglucan[which(Totoaba.GOenrichment.BP.Bglucan$Totoaba.Bglucan.topGO.BP.padj<=0.05),]
names(Totoaba.GOenrichment.BP.Bglucan.BH)[names(Totoaba.GOenrichment.BP.Bglucan.BH) == "Totoaba.Bglucan.topGO.BP.padj"] <- "BH.adjustedpvalue"
Totoaba.Bglucan.topGO.MF.padj <- round(p.adjust(Totoaba.Bglucan.topGO.MF.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.MF.Bglucan <- cbind(Totoaba.Bglucan.topGO.MF.results,Totoaba.Bglucan.topGO.MF.padj)
Totoaba.GOenrichment.MF.Bglucan <- Totoaba.GOenrichment.MF.Bglucan[order(Totoaba.GOenrichment.MF.Bglucan$Totoaba.Bglucan.topGO.MF.padj),]
Totoaba.GOenrichment.MF.Bglucan.BH <- Totoaba.GOenrichment.MF.Bglucan[which(Totoaba.GOenrichment.MF.Bglucan$Totoaba.Bglucan.topGO.MF.padj<=0.05),]
names(Totoaba.GOenrichment.MF.Bglucan.BH)[names(Totoaba.GOenrichment.MF.Bglucan.BH) == "Totoaba.Bglucan.topGO.MF.padj"] <- "BH.adjustedpvalue"
Totoaba.Bglucan.topGO.CC.padj <- round(p.adjust(Totoaba.Bglucan.topGO.CC.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.CC.Bglucan <- cbind(Totoaba.Bglucan.topGO.CC.results,Totoaba.Bglucan.topGO.CC.padj)
Totoaba.GOenrichment.CC.Bglucan <- Totoaba.GOenrichment.CC.Bglucan[order(Totoaba.GOenrichment.CC.Bglucan$Totoaba.Bglucan.topGO.CC.padj),]
Totoaba.GOenrichment.CC.Bglucan.BH <- Totoaba.GOenrichment.CC.Bglucan[which(Totoaba.GOenrichment.CC.Bglucan$Totoaba.Bglucan.topGO.CC.padj<=0.05),]
names(Totoaba.GOenrichment.CC.Bglucan.BH)[names(Totoaba.GOenrichment.CC.Bglucan.BH) == "Totoaba.Bglucan.topGO.CC.padj"] <- "BH.adjustedpvalue"
write.table(Totoaba.GOenrichment.BP.Bglucan, file="Totoaba.GOenrichment.BP.Bglucan.txt", sep="\t")
write.table(Totoaba.GOenrichment.BP.Bglucan.BH, file="Totoaba.GOenrichment.BP.Bglucan.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.Bglucan, file="Totoaba.GOenrichment.MF.Bglucan.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.Bglucan.BH, file="Totoaba.GOenrichment.MF.Bglucan.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.Bglucan, file="Totoaba.GOenrichment.CC.Bglucan.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.Bglucan.BH, file="Totoaba.GOenrichment.CC.Bglucan.BH.txt", sep="\t")



# FUNCTIONAL GO ENRICHMENT ANALYSIS FOR INULIN TREATMENT #

# Import GO annotations for the genes
Totoaba.GOterms <- readMappings(file="Tmacdonaldi.GO-annotations.txt")

# Define list of genes of interest
Totoaba.GeneIDs <- names(Totoaba.GOterms)
Totoaba.DEGs_Inulin <- read.table("Totoaba.DEGs.DESeq2.Inulin.txt", header=FALSE)
Totoaba.DEGs_Inulin <- as.character(Totoaba.DEGs_Inulin$V1)
DEGs.Inulin <- factor(as.integer(Totoaba.GeneIDs%in%Totoaba.DEGs_Inulin))
names(DEGs.Inulin) <- Totoaba.GeneIDs

# Put data together into an R object #
Totoaba.GO_BP.Inulin <- new("topGOdata", description="Totoaba-Inulin (BP)", ontology="BP", allGenes=DEGs.Inulin, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_MF.Inulin <- new("topGOdata", description="Totoaba-Inulin (MF)", ontology="MF", allGenes=DEGs.Inulin, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_CC.Inulin <- new("topGOdata", description="Totoaba-Inulin (CC)", ontology="CC", allGenes=DEGs.Inulin, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)

# Perform of GO enrichment test #
Totoaba.Inulin.BP.Fisher <- runTest(Totoaba.GO_BP.Inulin, algorithm="weight01", statistic="fisher")
Totoaba.Inulin.MF.Fisher <- runTest(Totoaba.GO_MF.Inulin, algorithm="weight01", statistic="fisher")
Totoaba.Inulin.CC.Fisher <- runTest(Totoaba.GO_CC.Inulin, algorithm="weight01", statistic="fisher")
Totoaba.Inulin.allGO.BP <- usedGO(Totoaba.GO_BP.Inulin)
Totoaba.Inulin.topGO.BP.results <- GenTable(Totoaba.GO_BP.Inulin, weightFisher=Totoaba.Inulin.BP.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Inulin.allGO.BP))
names(Totoaba.Inulin.topGO.BP.results)[names(Totoaba.Inulin.topGO.BP.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.Inulin.allGO.MF <- usedGO(Totoaba.GO_MF.Inulin)
Totoaba.Inulin.topGO.MF.results <- GenTable(Totoaba.GO_MF.Inulin, weightFisher=Totoaba.Inulin.MF.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Inulin.allGO.MF))
names(Totoaba.Inulin.topGO.MF.results)[names(Totoaba.Inulin.topGO.MF.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.Inulin.allGO.CC <- usedGO(Totoaba.GO_CC.Inulin)
Totoaba.Inulin.topGO.CC.results <- GenTable(Totoaba.GO_CC.Inulin, weightFisher=Totoaba.Inulin.CC.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Inulin.allGO.CC))
names(Totoaba.Inulin.topGO.CC.results)[names(Totoaba.Inulin.topGO.CC.results) == "weightFisher"] <- "weightFisher.pvalue"

# Correct for multiple testing - adjusted p-values #
Totoaba.Inulin.topGO.BP.padj <- round(p.adjust(Totoaba.Inulin.topGO.BP.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.BP.Inulin <- cbind(Totoaba.Inulin.topGO.BP.results,Totoaba.Inulin.topGO.BP.padj)
Totoaba.GOenrichment.BP.Inulin <- Totoaba.GOenrichment.BP.Inulin[order(Totoaba.GOenrichment.BP.Inulin$Totoaba.Inulin.topGO.BP.padj),]
Totoaba.GOenrichment.BP.Inulin.BH <- Totoaba.GOenrichment.BP.Inulin[which(Totoaba.GOenrichment.BP.Inulin$Totoaba.Inulin.topGO.BP.padj<=0.05),]
names(Totoaba.GOenrichment.BP.Inulin.BH)[names(Totoaba.GOenrichment.BP.Inulin.BH) == "Totoaba.Inulin.topGO.BP.padj"] <- "BH.adjustedpvalue"
Totoaba.Inulin.topGO.MF.padj <- round(p.adjust(Totoaba.Inulin.topGO.MF.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.MF.Inulin <- cbind(Totoaba.Inulin.topGO.MF.results,Totoaba.Inulin.topGO.MF.padj)
Totoaba.GOenrichment.MF.Inulin <- Totoaba.GOenrichment.MF.Inulin[order(Totoaba.GOenrichment.MF.Inulin$Totoaba.Inulin.topGO.MF.padj),]
Totoaba.GOenrichment.MF.Inulin.BH <- Totoaba.GOenrichment.MF.Inulin[which(Totoaba.GOenrichment.MF.Inulin$Totoaba.Inulin.topGO.MF.padj<=0.05),]
names(Totoaba.GOenrichment.MF.Inulin.BH)[names(Totoaba.GOenrichment.MF.Inulin.BH) == "Totoaba.Inulin.topGO.MF.padj"] <- "BH.adjustedpvalue"
Totoaba.Inulin.topGO.CC.padj <- round(p.adjust(Totoaba.Inulin.topGO.CC.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.CC.Inulin <- cbind(Totoaba.Inulin.topGO.CC.results,Totoaba.Inulin.topGO.CC.padj)
Totoaba.GOenrichment.CC.Inulin <- Totoaba.GOenrichment.CC.Inulin[order(Totoaba.GOenrichment.CC.Inulin$Totoaba.Inulin.topGO.CC.padj),]
Totoaba.GOenrichment.CC.Inulin.BH <- Totoaba.GOenrichment.CC.Inulin[which(Totoaba.GOenrichment.CC.Inulin$Totoaba.Inulin.topGO.CC.padj<=0.05),]
names(Totoaba.GOenrichment.CC.Inulin.BH)[names(Totoaba.GOenrichment.CC.Inulin.BH) == "Totoaba.Inulin.topGO.CC.padj"] <- "BH.adjustedpvalue"
write.table(Totoaba.GOenrichment.BP.Inulin, file="Totoaba.GOenrichment.BP.Inulin.txt", sep="\t")
write.table(Totoaba.GOenrichment.BP.Inulin.BH, file="Totoaba.GOenrichment.BP.Inulin.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.Inulin, file="Totoaba.GOenrichment.MF.Inulin.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.Inulin.BH, file="Totoaba.GOenrichment.MF.Inulin.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.Inulin, file="Totoaba.GOenrichment.CC.Inulin.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.Inulin.BH, file="Totoaba.GOenrichment.CC.Inulin.BH.txt", sep="\t")



# FUNCTIONAL GO ENRICHMENT ANALYSIS FOR CHITOSAN TREATMENT #

# Import GO annotations for the genes
Totoaba.GOterms <- readMappings(file="Tmacdonaldi.GO-annotations.txt")

# Define list of genes of interest
Totoaba.GeneIDs <- names(Totoaba.GOterms)
Totoaba.DEGs_Chitosan <- read.table("Totoaba.DEGs.DESeq2.Chitosan.txt", header=FALSE)
Totoaba.DEGs_Chitosan <- as.character(Totoaba.DEGs_Chitosan$V1)
DEGs.Chitosan <- factor(as.integer(Totoaba.GeneIDs%in%Totoaba.DEGs_Chitosan))
names(DEGs.Chitosan) <- Totoaba.GeneIDs

# Put data together into an R object #
Totoaba.GO_BP.Chitosan <- new("topGOdata", description="Totoaba-Chitosan (BP)", ontology="BP", allGenes=DEGs.Chitosan, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_MF.Chitosan <- new("topGOdata", description="Totoaba-Chitosan (MF)", ontology="MF", allGenes=DEGs.Chitosan, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_CC.Chitosan <- new("topGOdata", description="Totoaba-Chitosan (CC)", ontology="CC", allGenes=DEGs.Chitosan, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)

# Perform of GO enrichment test #
Totoaba.Chitosan.BP.Fisher <- runTest(Totoaba.GO_BP.Chitosan, algorithm="weight01", statistic="fisher")
Totoaba.Chitosan.MF.Fisher <- runTest(Totoaba.GO_MF.Chitosan, algorithm="weight01", statistic="fisher")
Totoaba.Chitosan.CC.Fisher <- runTest(Totoaba.GO_CC.Chitosan, algorithm="weight01", statistic="fisher")
Totoaba.Chitosan.allGO.BP <- usedGO(Totoaba.GO_BP.Chitosan)
Totoaba.Chitosan.topGO.BP.results <- GenTable(Totoaba.GO_BP.Chitosan, weightFisher=Totoaba.Chitosan.BP.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Chitosan.allGO.BP))
names(Totoaba.Chitosan.topGO.BP.results)[names(Totoaba.Chitosan.topGO.BP.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.Chitosan.allGO.MF <- usedGO(Totoaba.GO_MF.Chitosan)
Totoaba.Chitosan.topGO.MF.results <- GenTable(Totoaba.GO_MF.Chitosan, weightFisher=Totoaba.Chitosan.MF.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Chitosan.allGO.MF))
names(Totoaba.Chitosan.topGO.MF.results)[names(Totoaba.Chitosan.topGO.MF.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.Chitosan.allGO.CC <- usedGO(Totoaba.GO_CC.Chitosan)
Totoaba.Chitosan.topGO.CC.results <- GenTable(Totoaba.GO_CC.Chitosan, weightFisher=Totoaba.Chitosan.CC.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.Chitosan.allGO.CC))
names(Totoaba.Chitosan.topGO.CC.results)[names(Totoaba.Chitosan.topGO.CC.results) == "weightFisher"] <- "weightFisher.pvalue"

# Correct for multiple testing - adjusted p-values #
Totoaba.Chitosan.topGO.BP.padj <- round(p.adjust(Totoaba.Chitosan.topGO.BP.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.BP.Chitosan <- cbind(Totoaba.Chitosan.topGO.BP.results,Totoaba.Chitosan.topGO.BP.padj)
Totoaba.GOenrichment.BP.Chitosan <- Totoaba.GOenrichment.BP.Chitosan[order(Totoaba.GOenrichment.BP.Chitosan$Totoaba.Chitosan.topGO.BP.padj),]
Totoaba.GOenrichment.BP.Chitosan.BH <- Totoaba.GOenrichment.BP.Chitosan[which(Totoaba.GOenrichment.BP.Chitosan$Totoaba.Chitosan.topGO.BP.padj<=0.05),]
names(Totoaba.GOenrichment.BP.Chitosan.BH)[names(Totoaba.GOenrichment.BP.Chitosan.BH) == "Totoaba.Chitosan.topGO.BP.padj"] <- "BH.adjustedpvalue"
Totoaba.Chitosan.topGO.MF.padj <- round(p.adjust(Totoaba.Chitosan.topGO.MF.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.MF.Chitosan <- cbind(Totoaba.Chitosan.topGO.MF.results,Totoaba.Chitosan.topGO.MF.padj)
Totoaba.GOenrichment.MF.Chitosan <- Totoaba.GOenrichment.MF.Chitosan[order(Totoaba.GOenrichment.MF.Chitosan$Totoaba.Chitosan.topGO.MF.padj),]
Totoaba.GOenrichment.MF.Chitosan.BH <- Totoaba.GOenrichment.MF.Chitosan[which(Totoaba.GOenrichment.MF.Chitosan$Totoaba.Chitosan.topGO.MF.padj<=0.05),]
names(Totoaba.GOenrichment.MF.Chitosan.BH)[names(Totoaba.GOenrichment.MF.Chitosan.BH) == "Totoaba.Chitosan.topGO.MF.padj"] <- "BH.adjustedpvalue"
Totoaba.Chitosan.topGO.CC.padj <- round(p.adjust(Totoaba.Chitosan.topGO.CC.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.CC.Chitosan <- cbind(Totoaba.Chitosan.topGO.CC.results,Totoaba.Chitosan.topGO.CC.padj)
Totoaba.GOenrichment.CC.Chitosan <- Totoaba.GOenrichment.CC.Chitosan[order(Totoaba.GOenrichment.CC.Chitosan$Totoaba.Chitosan.topGO.CC.padj),]
Totoaba.GOenrichment.CC.Chitosan.BH <- Totoaba.GOenrichment.CC.Chitosan[which(Totoaba.GOenrichment.CC.Chitosan$Totoaba.Chitosan.topGO.CC.padj<=0.05),]
names(Totoaba.GOenrichment.CC.Chitosan.BH)[names(Totoaba.GOenrichment.CC.Chitosan.BH) == "Totoaba.Chitosan.topGO.CC.padj"] <- "BH.adjustedpvalue"
write.table(Totoaba.GOenrichment.BP.Chitosan, file="Totoaba.GOenrichment.BP.Chitosan.txt", sep="\t")
write.table(Totoaba.GOenrichment.BP.Chitosan.BH, file="Totoaba.GOenrichment.BP.Chitosan.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.Chitosan, file="Totoaba.GOenrichment.MF.Chitosan.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.Chitosan.BH, file="Totoaba.GOenrichment.MF.Chitosan.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.Chitosan, file="Totoaba.GOenrichment.CC.Chitosan.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.Chitosan.BH, file="Totoaba.GOenrichment.CC.Chitosan.BH.txt", sep="\t")



# FUNCTIONAL GO ENRICHMENT ANALYSIS FOR ALL TREATMENTS TREATMENT #

# Import GO annotations for the genes
Totoaba.GOterms <- readMappings(file="Tmacdonaldi.GO-annotations.txt")

# Define list of genes of interest
Totoaba.GeneIDs <- names(Totoaba.GOterms)
Totoaba.DEGs_alltreatments <- read.table("Totoaba.DEGs.DESeq2.alltreatments.txt", header=FALSE)
Totoaba.DEGs_alltreatments <- as.character(Totoaba.DEGs_alltreatments$V1)
DEGs.alltreatments <- factor(as.integer(Totoaba.GeneIDs%in%Totoaba.DEGs_alltreatments))
names(DEGs.alltreatments) <- Totoaba.GeneIDs

# Put data together into an R object #
Totoaba.GO_BP.alltreatments <- new("topGOdata", description="Totoaba-alltreatments (BP)", ontology="BP", allGenes=DEGs.alltreatments, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_MF.alltreatments <- new("topGOdata", description="Totoaba-alltreatments (MF)", ontology="MF", allGenes=DEGs.alltreatments, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)
Totoaba.GO_CC.alltreatments <- new("topGOdata", description="Totoaba-alltreatments (CC)", ontology="CC", allGenes=DEGs.alltreatments, annot=annFUN.gene2GO, gene2GO=Totoaba.GOterms)

# Perform of GO enrichment test #
Totoaba.alltreatments.BP.Fisher <- runTest(Totoaba.GO_BP.alltreatments, algorithm="weight01", statistic="fisher")
Totoaba.alltreatments.MF.Fisher <- runTest(Totoaba.GO_MF.alltreatments, algorithm="weight01", statistic="fisher")
Totoaba.alltreatments.CC.Fisher <- runTest(Totoaba.GO_CC.alltreatments, algorithm="weight01", statistic="fisher")
Totoaba.alltreatments.allGO.BP <- usedGO(Totoaba.GO_BP.alltreatments)
Totoaba.alltreatments.topGO.BP.results <- GenTable(Totoaba.GO_BP.alltreatments, weightFisher=Totoaba.alltreatments.BP.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.alltreatments.allGO.BP))
names(Totoaba.alltreatments.topGO.BP.results)[names(Totoaba.alltreatments.topGO.BP.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.alltreatments.allGO.MF <- usedGO(Totoaba.GO_MF.alltreatments)
Totoaba.alltreatments.topGO.MF.results <- GenTable(Totoaba.GO_MF.alltreatments, weightFisher=Totoaba.alltreatments.MF.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.alltreatments.allGO.MF))
names(Totoaba.alltreatments.topGO.MF.results)[names(Totoaba.alltreatments.topGO.MF.results) == "weightFisher"] <- "weightFisher.pvalue"
Totoaba.alltreatments.allGO.CC <- usedGO(Totoaba.GO_CC.alltreatments)
Totoaba.alltreatments.topGO.CC.results <- GenTable(Totoaba.GO_CC.alltreatments, weightFisher=Totoaba.alltreatments.CC.Fisher, orderBy="weightFisher", topNodes=length(Totoaba.alltreatments.allGO.CC))
names(Totoaba.alltreatments.topGO.CC.results)[names(Totoaba.alltreatments.topGO.CC.results) == "weightFisher"] <- "weightFisher.pvalue"

# Correct for multiple testing - adjusted p-values #
Totoaba.alltreatments.topGO.BP.padj <- round(p.adjust(Totoaba.alltreatments.topGO.BP.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.BP.alltreatments <- cbind(Totoaba.alltreatments.topGO.BP.results,Totoaba.alltreatments.topGO.BP.padj)
Totoaba.GOenrichment.BP.alltreatments <- Totoaba.GOenrichment.BP.alltreatments[order(Totoaba.GOenrichment.BP.alltreatments$Totoaba.alltreatments.topGO.BP.padj),]
Totoaba.GOenrichment.BP.alltreatments.BH <- Totoaba.GOenrichment.BP.alltreatments[which(Totoaba.GOenrichment.BP.alltreatments$Totoaba.alltreatments.topGO.BP.padj<=0.05),]
names(Totoaba.GOenrichment.BP.alltreatments.BH)[names(Totoaba.GOenrichment.BP.alltreatments.BH) == "Totoaba.alltreatments.topGO.BP.padj"] <- "BH.adjustedpvalue"
Totoaba.alltreatments.topGO.MF.padj <- round(p.adjust(Totoaba.alltreatments.topGO.MF.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.MF.alltreatments <- cbind(Totoaba.alltreatments.topGO.MF.results,Totoaba.alltreatments.topGO.MF.padj)
Totoaba.GOenrichment.MF.alltreatments <- Totoaba.GOenrichment.MF.alltreatments[order(Totoaba.GOenrichment.MF.alltreatments$Totoaba.alltreatments.topGO.MF.padj),]
Totoaba.GOenrichment.MF.alltreatments.BH <- Totoaba.GOenrichment.MF.alltreatments[which(Totoaba.GOenrichment.MF.alltreatments$Totoaba.alltreatments.topGO.MF.padj<=0.05),]
names(Totoaba.GOenrichment.MF.alltreatments.BH)[names(Totoaba.GOenrichment.MF.alltreatments.BH) == "Totoaba.alltreatments.topGO.MF.padj"] <- "BH.adjustedpvalue"
Totoaba.alltreatments.topGO.CC.padj <- round(p.adjust(Totoaba.alltreatments.topGO.CC.results$weightFisher.pvalue, method="BH"), digits=6)
Totoaba.GOenrichment.CC.alltreatments <- cbind(Totoaba.alltreatments.topGO.CC.results,Totoaba.alltreatments.topGO.CC.padj)
Totoaba.GOenrichment.CC.alltreatments <- Totoaba.GOenrichment.CC.alltreatments[order(Totoaba.GOenrichment.CC.alltreatments$Totoaba.alltreatments.topGO.CC.padj),]
Totoaba.GOenrichment.CC.alltreatments.BH <- Totoaba.GOenrichment.CC.alltreatments[which(Totoaba.GOenrichment.CC.alltreatments$Totoaba.alltreatments.topGO.CC.padj<=0.05),]
names(Totoaba.GOenrichment.CC.alltreatments.BH)[names(Totoaba.GOenrichment.CC.alltreatments.BH) == "Totoaba.alltreatments.topGO.CC.padj"] <- "BH.adjustedpvalue"
write.table(Totoaba.GOenrichment.BP.alltreatments, file="Totoaba.GOenrichment.BP.alltreatments.txt", sep="\t")
write.table(Totoaba.GOenrichment.BP.alltreatments.BH, file="Totoaba.GOenrichment.BP.alltreatments.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.alltreatments, file="Totoaba.GOenrichment.MF.alltreatments.txt", sep="\t")
write.table(Totoaba.GOenrichment.MF.alltreatments.BH, file="Totoaba.GOenrichment.MF.alltreatments.BH.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.alltreatments, file="Totoaba.GOenrichment.CC.alltreatments.txt", sep="\t")
write.table(Totoaba.GOenrichment.CC.alltreatments.BH, file="Totoaba.GOenrichment.CC.alltreatments.BH.txt", sep="\t")

## END OF THE SCRIPT!