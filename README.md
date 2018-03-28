# TCGA-CPTAC-analysis
Multi-omics analysis of breast and ovarian cancer data from National Cancer Institute's TCGA (The Cancer Genome Atlas) and CPTAC (Clinical Proteomic Tumor Analysis Consortium) repositories.

## Motivation
Set up in-house scripts for analyzing molecular data from patient tumor samples. Data on breast and ovarian tumor samples characterized by phospho-proteomics and transcriptomics were extracted from the CPTAC and TCGA public repositories. Kinase activities were predicted from phospho-proteomics data using the Kinase Set Enrichment Analysis (KSEA) proposed by Casado et al. (*Science Signaling* 2013).

## File Description
+ *BRCA_networKIN_input.tsv*: Formatted files for NetworKIN kinase-substrate predictions (SCRIPT NOT YET INCORPORATES PREDICTIONS).
+ *BRCA_networKIN_input_mod.tsv*: Formatted files for NetworKIN kinase-substrate predictions (SCRIPT NOT YET INCORPORATES PREDICTIONS).
+ *PhosphoSitePlus_Kinase_Substrate_Dataset*: Kinase-substrate relationships from PhosphoSitePlus online database.
+ *TCGA_CPTAC_analysis.R*: R script for data anlysis.
+ *breast_logcpm.csv*: RNAseq data for breast cancer samples.
+ *ovarian_logcpm.csv*: RNAseq data for ovarian cancer samples.
+ *pam50_annotation.txt*: List of 50 breast cancer gene signatures for stratifying patient samples.
