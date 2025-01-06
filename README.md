# Physiological and transcriptomic effects of formulated diets including the prebiotics inulin, β-glucan, and chitosan on juveniles of *Totoaba macdonaldi* - Supplementary Data

This GitHub repository provides the supplementary data referenced in the publication cited below.

## How to use and cite these files 

All files are publicly available and can be used for further research or other applications. However, if you utilize these resources in your work, we kindly request that you cite our original publication.

**Physiological and transcriptomic effects of formulated diets including the prebiotics inulin, β-glucan, and chitosan on juveniles of *Totoaba macdonaldi.*** Oscar E. Juárez, Clara E. Galindo-Sánchez, Fabiola Lafarga-De la Cruz, Sara Enciso, Edgar A. López-Landavery, Camilo Muñoz, Felipe Aguilera, Juan Pablo Lazo. *Aquaculture International* **32**:61-85 (2024). https://doi.org/10.1007/s10499-023-01144-1

**Abstract**

In this study, we evaluated the effects of three prebiotics (inulin, β-glucan, and chitosan) on the physiological performance of Totoaba macdonaldi juveniles under culture conditions. The respiratory burst and the leucocyte content were measured in the blood to assess innate immune responses. The intestinal digestive capacity was evaluated by analyzing trypsin, amylase, and lipase activities, whereas the effects of such prebiotics at the transcriptomic level were assessed by implementing the RNA-Seq of liver tissue. After 60 days, fish fed with 0.5% chitosan diets showed the highest respiratory burst, the lowest lipase activity, and the highest number of differentially expressed genes (DEGs), where biological processes related to proteolysis, digestion, and lipid hydroxylation were the most affected. In addition, fish from the chitosan diet showed the highest expression of immunoglobulin genes. In contrast, fish fed with the 1% inulin diet presented the highest diet digestibility and trypsin and lipase activities. These physiological effects align with the highest expression of trypsin-like and chymotrypsin-like genes in the liver of fish from this diet. On the other hand, fish fed the 0.1% β-glucan diets showed the lowest amount of DEGs compared to the control group, most of which were associated with immune response, with an up-regulation of genes related to the complement system and a downregulation of immunoglobulin genes. Based on our results, we propose the inclusion of 1% dietary inulin to improve the digestibility of experimental diets and the addition of 0.5% chitosan to stimulate the immune system of T. macdonaldi juveniles.

## Author contact

- [Fabiola Lafarga-de la Cruz](mailto:flafarga@cicese.mx) (corresponding author)
- [Felipe Aguilera](mailto:faguilera@udec.cl) (Assistant professor)
- [Juan Pablo Lazo](mailto:jplazo@cicese.mx) (corresponding author)

## How to download data files

The files in this repository are ready for use. Simply click on the file, and you will be directed to the GitHub webpage to save it to your device.

## Index of data file contents

[01-Totoaba-transcriptome-BUSCO.docx](https://github.com/faguil/Totoaba_transcriptomics/blob/main/01-Totoaba-transcriptome-BUSCO.docx) contains the BUSCO results of the quality assessment of the Totoaba reference transcriptome.

[02-Trinotate-report-DEGs.docx](https://github.com/faguil/Totoaba_transcriptomics/blob/main/02-Trinotate-report-DEGs.docx) contains the Trinotate results of the differentially expressed genes used to perform enriched Gene Ontology heatmaps.

[03-Phylogenetic-tree-ETS4_PDEF.pdf](https://github.com/faguil/Pearl-Sac-Gene-Expression/blob/main/03-Phylogenetic-tree-ETS4_PDEF.pdf) contains the Maximum Likelihood phylogenetic analysis of the *Pinctada maxima* ETS4/PDEF protein.

[04-Supplementary-Tables.xlsx](https://github.com/faguil/Pearl-Sac-Gene-Expression/blob/main/04-Supplementary-Tables.xlsx) contains supplementary tables used in this study.

- **Table S1.** Transcripts differentially expressed between sacs producing high vs low lustre pearls. 
- **Table S2.** Transcripts differentially expressed between sacs producing pearls with 'calcification' and those without.
- **Table S3.** Gene Ontology 'Biological Process' enrichments within calcification-associated differentially expressed transcripts.
- **Table S4.** Transcripts differentially expressed between sacs producing pearls with 'underskin' and those without.
- **Table S5.** Transcripts differentially expressed between sacs producing pearls with different weights.



# Totoaba_transcriptomics
Identification of differentially expressed genes in Totoaba juveniles exposed to different prebiotics

Public repository for bioinformatics analysis conducted in the evaluation of different prebiotics (i.e., chitosan, inulin, beta-glucan) on juveniles of Totoaba macdonaldi

This repository contains:

* Bash script to make a fake annotation file of the de novo reference transcriptome of Totoaba macdonaldi (filename = Totoaba.fake-annotation-transcriptome.sh)
* R script to conduct the identification of differentially expressed genes (DEGs) using DESeq2 (filename = Totoaba.DESeq2.R)
* R script to conduct the identification of differentially expressed genes (DEGs) using edgeR (filename = Totoaba.edgeR.R
* R script to conduct the gene ontology (GO) enrichment of DEGs using topGO (filename = topGO.R)
