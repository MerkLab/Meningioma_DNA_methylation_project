# Meningioma_DNA_methylation_project
Our study highlights the role of a DNA hypermethylation signature and its downstream effects on signaling pathways involving clustered protocadherin genes during progression of meningiomas.

This study includes DNA methylation array data from a cross-sectional meningioma cohort (n=231), and a longitudinal cohort including primary and recurrent meningioma samples from a total of 18 patients.

All meningiomas and meningioma cell models were profiled using the MethylationEPIC v2.0 array. Raw data are archived at GEO (GSE304097), and meta data is available from the corresponding publication listed below. The following datasets were used for further validation: meningioma RNAseq data from UCSF (GSE183656), meningioma DNA methylation arrays from UCSF (GSE183656), and meningioma DNA methylation arrays from Toronto (GSE180061). DNA methylation arrays from normal meninges used in this study are subject to legal restrictions and access decision on a case-by-case basis, and can be retrieved from Dr. Tiit Mathiesen from the Rigshospitalet in Copenhagen, Denmark (tiit.illimar.mathiesen@regionh.dk). Additional intermediate data files necessary or hepful for reproduction of this study can be found at Figshare (https://figshare.com/s/ea7852e5723c1972329d).

This repository provides essential codes used: 

I) to identify a DNA hypermethylation signature associated with disease progression in meninigioma
II) to infer copy number variations in benign and malignant meningiomas
III) to perform survival analyses and assess prediction performances of Cox regression models
IV) to investigate differential methylation in between normal meninges and distinct meningioma settings
V) to perform trajectory inference in order to analyze DNA methylation profiles in the context of microevolutionary development from normal meninges to malignant meningioma
VI) to perform differential gene expression analysis of benign and malignant meningiomas

The code of this study has been subdivided into 12 sub-sections that individually fulfill distinct aims within the course of the study. Following, the main tasks of each code sub-section are highlighted:

01_QC-preselect: Several measures are employed to screen a total of ??? DNA methylation arrays from human meningioma as initial quality control, including probe signal intensities and signal-to-noise ratios. Output is a total of 231 human meningioma samples with high-quality DNA methylation arrays included in the discovery cohort.

02_discovery_benchmarking: Samples from the discovery cohort were validated to recapitulate know molecular and clinical features of human meningioma, including immune cell infiltration and survival stratification using molecular groups. 

03_discovery_clustering: This code identifies a DNA hypermethylation signature associated with risk of progression in human meningioma using principal component analysis for probe selection. Further, relation of this hypermethylation signature to previously established classification schemes is outlined. Unsupervised hierarchical clustering is used to stratify the discovery cohort into METHlow an METHhigh human meningioma, were METHlow mainly encompass benign meningioma groups (CNS WHO grade 1, benign, NF2-wildtype, immuneenriched), while METHhigh tumors show a strong enrichment of malignant meningiomas (CNS WHO grade 3, malignant, proliferative). Main output is the stratification of the discovery cohort into DNA methylation clusters.

04_machine_learning_prediction: This code provides all machine learning approaches to predict molecular groups defined by distinct groups and methylation clustering. Mainstay of this work is done by the caret R package.

05_differential_methylation: This part includes all code used to define DNA methylation changes between any given contrast (e.g. CNS WHO grade 1 vs 3) by using linear modeling, and their statistical analysis. Furthermore, for methylation clusters, a detailed gene-centric analysis (overenrichment) of differentially regulated CpG sites is performed, implicating several gene clusters (HOX genes, protocadherins) during progression of meningioma. This sub-section provides differentially methylated CpG sites defining distinct subgroups of meningioma (CNS WHO grade, molecular groups, methylation clusters).

06_suvival: This code accordingly provides detailed analyses of progression-free survival as stratified by the same contrasts used for differential DNA methylation analyses. This includes calculation of standard Kaplan-Meier curves, modeling of survival using Cox proportional hazards which accounts for competing variables of clinical meningioma behavior, and several measurements of model performance including but not limited to c-index, net reclassification improvement, and Brier scores. While these findings show high predictive power for our methylation cluster signature, they also provide evidence for a nonlinear, continuous risk gradient in meningioma that is directly mirrored by our signature. 

07_trajectory_analysis: This code investigates in detail the transition of human meningiomas, both with respect to benign vs malignant as well as primary vs recurrent, along the risk gradient, and calculates a disease trajectory that is strongly mirrored by the methylation signature.

08_cell_lines: Several cell models stemming from both benign and malignant meningiomas were investigated with regard to their ability to faithfully recapitulate the findings of our DNA methylation clusters in primary tumors. Furthermore, code encompasses data analyses from FACS analyses investigating cellular localization of beta catenin in meningioma cell models, and their graphical representation.

09_discovery_CNV: Copy number variations were analyzed based on DNA methylation array data using the conumee2 R package. This includes generation of genomplots and their corresponding .tsv files indicating bin-wise log2 ratios of signal intensities. 

10_longitudinal_CNV: Similar data analyses to 09_discovery_CNV. Furthermore, for a more gene-centric view of the data, code includes analyses using the CNV.focal function from conumee2 in order to define gene-level CNVs.


11_RNAseq: Standard differential gene expression analyses using DESeq2. These were performed for in-house meningiomas, stratified by DNA methylation cluster assignment (n = 10). Furthermore, this was further evaluated using the UCSF gene expression cohort comprising a total of 185 human meningiomas. Output is differential gene expression in between METHhigh and METHlow meningioma clusters.

12_classifier: Using the caret R package to generate a random forest model used as a binary classifier for METHhigh and METHlow meningioma clusters based in DNA methylation data from the Discovery cohort. The code also details specifics for the shiny web app hosted on the Hertie Institute server. 

Last, we provide an easy-to-use web application that predicts hypermethylation class association and genome-wide copy number analysis at https://mmcc.neurologie.uni-tuebingen.de. The code to this web app can be found within the 12_classifier code file.

Please find the version numbers of the most relevant R packages used within our code below. All packages are available at Bioconductor, except conumee2, which was cloned from Github (https://github.com/hovestadtlab/conumee2). All calculations were performed on a MacBook Pro, Apple M1 Chip, MacOS Sequoia 15.6.1.
R (v4.4.2)
SeSAMe (v1.24.0)
sesameData (v1.24.0)
IlluminaHumanMethylationEPICv2anno.20a1.hg38 (v1.0.0)
minfi (v1.52.1)
caret (v7.0-1)
factoextra (v1.0.7)
conumee2 (v2.1.2)
SummarizedExperiment (v1.26.0)
slingshot (v2.14.0)
ConsensusClusterPlus (v1.70.0)
survival (v3.8-3)
survminer (v0.5.0)
DESeq2 (v1.46.0)
limma (v3.62.2)
pheatmap (v1.0.13)
ggplot2 (v4.0.1)
ggpubr (v0.6.2)
networkD3 (v0.4.1)
lme4 (v1.1-37)
lmerTest (v3.1-3)
survIDINRI (1.1-2)

Please find the original publication here:





