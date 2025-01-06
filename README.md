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

## How to download code files

The code in this repository is stored in the [00-Code folder]( and can be downloaded as outlined below. The ***protein_motif_searching.pl*** script, written in Perl, scans one or more protein motif patterns in FASTA files. For detailed instructions on how to use the script, please refer to the README file in the 00-Code folder. While we believe the script and accompanying README make it largely self-explanatory, do not hesitate to reach out if you have any questions or concerns.

## How to download data files

The files in this repository are ready for use. Simply click on the file, and you will be directed to the GitHub webpage to save it to your device.

## Index of data file contents

[01-Totoaba-transcriptome-BUSCO.docx](https://github.com/faguil/Totoaba_transcriptomics/blob/main/01-Totoaba-transcriptome-BUSCO.docx) contains the BUSCO results of the quality assessment of the Totoaba reference transcriptome.

[02-Trinotate-report-DEGs.docx](https://github.com/faguil/Totoaba_transcriptomics/blob/main/02-Trinotate-report-DEGs.docx) contains the Trinotate results of the differentially expressed genes used to perform enriched Gene Ontology heatmaps.

[03-Gene-ontology-enrichment-DEGs.docx](https://github.com/faguil/Totoaba_transcriptomics/blob/main/03-Gene-ontology-enrichment-DEGs.docx) contains the complete Gene Ontology (GO) enrichment report for differentially expressed genes.

[04-Supplementary-Material.docx](https://github.com/faguil/Totoaba_transcriptomics/blob/main/04-Supplementary-Material.docx) contains supplementary material used in this study.

- **Figure 1.** Evaluation of independent sequence similarity search results across different protein databases. 
- **Figure 2.** Evaluation of GO and KEGG assignment.

[05-Totoaba-DESeq2.R](https://github.com/faguil/Totoaba_transcriptomics/blob/main/05-Totoaba.DESeq2.R) contains an R script to conduct the identification of differentially expressed genes (DEGs) using DESeq2.

[06-Totoaba-edgeR.R](https://github.com/faguil/Totoaba_transcriptomics/blob/main/06-Totoaba.edgeR.R) contains an R script to conduct the identification of differentially expressed genes (DEGs) using edgeR.

[07-Totoaba-fake-annotation-transcriptome.sh](https://github.com/faguil/Totoaba_transcriptomics/blob/main/07-Totoaba.fake-annotation-transcriptome.sh) contains a Bash script to make a fake annotation file of the de novo reference transcriptome of *Totoaba macdonaldi*.

[08-Totoaba-topGO.R](https://github.com/faguil/Totoaba_transcriptomics/blob/main/08-Totoaba.topGO.R) contains an R script to conduct the gene ontology (GO) enrichment of DEGs using topGO.
