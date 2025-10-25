# New insights into the molecular actions of Grosheimin, Costunolide, α- and β-Cyclocostunolide on primary cilia structure and Hedgehog signaling

This repository contains the data analysis workflow for the manuscript titled: "New insights into the molecular actions of Grosheimin, Costunolide, α- and β- Cyclocostunolide on primary cilia structure and Hedgehog signaling".

This study investigates the effects of four sesquiterpene lactones (SLs)—Grosheimin, Costunolide, α-Cyclocostunolide, and β-Cyclocostunolide on the structure of primary cilia and the Hedgehog signaling pathway in human primary fibroblasts and retinal pigment epithelial cells (RPE).

---

## Analysis

The `Analysis/` folder contains the R scripts and code used to perform the differential gene expression analysis of hTERT-RPE-1 cells treated with either Grosheimin or α-Cyclocostunolide, compared to a DMSO control.

The analysis workflow includes the following steps:
1.  **Quality Control (QC):** Initial quality assessment of the raw FASTQ files.
2.  **Read Alignment:** Alignment of reads to the human reference genome (GRCh38/hg38) using STAR.
3.  **Quantification:** Gene-level quantification of first-strand paired-end fragments was performed using featureCounts.
4.  **Differential Expression:** Analysis was conducted using DESeq2 in R to identify significantly up- and down-regulated genes.
5.  **Gene Ontology (GO) Enrichment Analysis:** Over-representation analysis was performed using ClusterProfiler to identify enriched biological processes and molecular functions.
6.  **Figure Generation:** The code to generate the volcano plots and GO analysis figures presented in the manuscript is also included.

---

## Data

The raw sequencing data (FASTQ files), raw gene counts, and Transcripts Per Million (TPMs) files used for this analysis are publicly available at the **ArrayExpress** database.

*   **Accession Number:** E-MTAB-15939 (currently in curation).

---

