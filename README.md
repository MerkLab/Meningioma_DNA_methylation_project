# Meningioma_DNA_methylation_project
Our study highlights the role of a DNA hypermethylation signature and its downstream effects on signaling pathways involving clustered protocadherin genes during progression of meningiomas.

This study includes DNA methylation array data from a cross-sectional meningioma cohort (n=231), and a longitudinal cohort including primary and recurrent meningioma samples from a total of 18 patients.

All meningiomas were profiled using the MethylationEPIC v2.0 array. Raw data are archived at GEO, and meta data is available from the corresponding publication listed below. This repository provides essential codes used: 

I) to identify a DNA hypermethylation signature associated with disease progression in meninigioma
II) to infer copy number variations in benign and malignant meningiomas
III) to perform survival analyses and assess prediction performances of Cox regression models
IV) to investigate differential methylation in between normal meninges and distinct meningioma settings
V) to perform trajectory inference in order to analyze DNA methylation profiles in the context of microevolutionary development from normal meninges to malignant meningioma
VI) to perform differential gene expression analysis of benign and malignant meningiomas

Last, we provide an easy-to-use web application that predicts hypermethylation class association and genome-wide copy number analysis at https://mmcc.neurologie.uni-tuebingen.de.

Please find the version numbers of the most relevant R packages used within our code below. All calculations were peformed on a MacBook Pro, Apple M1 Chip, MacOS Sequoia 15.6.1.
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




