# Physiological and transcriptomic effects of formulated diets including the prebiotics inulin, β-glucan, and chitosan on juveniles of *Totoaba macdonaldi* - Supplementary Data

This GitHub repository provides the supplementary data referenced in the publication cited below.

## How to use and cite these files 

All files are publicly available and can be used for further research or other applications. However, if you utilize these resources in your work, we kindly request that you cite our original publication.

**Physiological and transcriptomic effects of formulated diets including the prebiotics inulin, β-glucan, and chitosan on juveniles of *Totoaba macdonaldi.*** Oscar E. Juárez, Clara E. Galindo-Sánchez, Fabiola Lafarga-De la Cruz, Sara Enciso, Edgar A. López-Landavery, Camilo Muñoz, Felipe Aguilera, Juan Pablo Lazo. *Aquaculture International* **32**:61-85 (2024). https://doi.org/10.1007/s10499-023-01144-1

**Abstract**

Pearls are highly prized biomineralized gemstones produced by molluscs. The appearance and mineralogy of cultured pearls can vary markedly, greatly affecting their commercial value. To begin to understand the role of pearl sacs—organs that form in host oysters from explanted mantle tissues that surround and synthesize pearls—we undertook transcriptomic analyses to identify genes that are differentially expressed in sacs producing pearls with different surface and structural characteristics. Our results indicate that gene expression profiles correlate with different pearl defects, suggesting that gene regulation in the pearl sac contributes to pearl appearance and quality. For instance, pearl sacs that produced pearls with surface non-lustrous calcification significantly down-regulate genes associated with cilia and microtubule function compared to pearl sacs giving rise to lustrous pearls. These results suggest that gene expression profiling can advance our understanding of processes that control biomineralization, which may be of direct value to the pearl industry, particularly in relation to defects that result in low value pearls.

## Author contact

- [Carmel McDougall](mailto:c.mcdougall@uq.edu.au) (senior author - PostDoc)
- [Felipe Aguilera](mailto:f.aguilera@uq.edu.au) (first author - PhD student)
- [Bernie Degnan](b.degnan@uq.edu.au) (corresponding author)

## How to download data files

The files in this repository are ready for use. Simply click on the file, and you will be directed to the GitHub webpage to save it to your device.

## Index of data file contents

[01-Photographs-pearls.pdf](https://github.com/faguil/Pearl-Sac-Gene-Expression/blob/main/01-Photographs-pearls.pdf) contains the photographs of pearls extracted from the pearl sacs used in this study.

[02-Phylogenetic-tree-FoxJ1.pdf](https://github.com/faguil/Pearl-Sac-Gene-Expression/blob/main/02-Phylogenetic-tree-FoxJ1.pdf) contains the Maximum Likelihood phylogenetic analysis of the *Pinctada maxima* FoxJ1 protein.

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
