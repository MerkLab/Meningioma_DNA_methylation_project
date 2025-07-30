library(sesame)
library(minfi)
library(caret)
library(ggstatsplot)
library(ggpubr)
library(gplots)
library(networkD3)
library(factoextra)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(slingshot)
library(mclust, quietly = TRUE)
library(BiocNeighbors)


#data with EPICv2 probes names
idat_dir = "directory"
targets_validation <- read.metharray.sheet(idat_dir, pattern="targets.csv")

rgSet <- read.metharray.exp(base = idat_dir, targets=targets_validation)
rgSet
sampleNames(rgSet) <- targets_validation$ID
rgSet
annotation(rgSet) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
rgSet

###Normalization
mSetSq <- preprocessNoob(rgSet) 
mSetSq

#filtering
pvalspOOBAH = openSesame(idat_dir, func = pOOBAH, return.pval = TRUE)
detP = as.data.frame(pvalspOOBAH)
detP <- detP[,match(targets_validation$Basename, colnames(detP))]
colnames(detP) = targets_validation$ID
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
keep <- rowSums(detP < 0.05) == ncol(mSetSq) 
table(keep)
#based on detP, 350,857 probes are excluded, 586,133 kept

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

#make MethylSet a GenomicRatioSet
mSetSqFlt = mapToGenome(mSetSqFlt)
mSetSqFlt = ratioConvert(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
#probes after SNP drop (n=571,977)

# tag sex chromosome probes for removal
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
data("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

keep <- !(featureNames(mSetSqFlt) %in% annoEPICv2$Name[annoEPICv2$chr %in% 
                                                         c("chrX","chrY")])
table(keep)
#10,053 probes in sex chromosomes

mSetSqFlt
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt
#final probe number n=561,924


#get B values
bVals <- getBeta(mSetSqFlt)

####meninges are on EPIC array, collapse validation beta probes first
bVals_mat = as.matrix(bVals)
bVals_mat = betasCollapseToPfx(bVals_mat)
bVals_mat = as.data.frame(bVals_mat)

#read in data for normal meninges, using the same pipeline as above for the validation cohort
#data with EPIC probes names
idat_dir = "directory"
targets_meninges <- read.metharray.sheet(idat_dir, pattern="targets_meninges.csv")

rgSet <- read.metharray.exp(base = idat_dir, targets=targets_meninges)
rgSet
sampleNames(rgSet) <- targets_meninges$ID
rgSet

###Normalization
mSetSq <- preprocessNoob(rgSet) 
mSetSq

#filtering
pvalspOOBAH = openSesame(idat_dir, func = pOOBAH, return.pval = TRUE)
detP = as.data.frame(pvalspOOBAH)
detP <- detP[,match(targets_meninges$Basename, colnames(detP))]
colnames(detP) = targets_meninges$ID
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
detP = na.omit(detP)
keep <- rowSums(detP < 0.05) == ncol(mSetSq) 
table(keep)
#based on detP, 83,080 probes are excluded, 783,137 kept

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

#make MethylSet a GenomicRatioSet
mSetSqFlt = mapToGenome(mSetSqFlt)
mSetSqFlt = ratioConvert(mSetSqFlt)
mSetSqFlt

#get B values
bVals_meninges <- getBeta(mSetSqFlt)
dim(bVals_meninges)
bVals_meninges = as.data.frame(bVals_meninges)

save.image()

#reduce validation and meninges beta values to same rownames, bring to same order, and combine
common_rows = intersect(rownames(bVals_mat), rownames(bVals_meninges))
bVals_mat_common = bVals_mat[common_rows,,drop=FALSE]
bVals_menninges_common = bVals_meninges[common_rows,,drop=FALSE]
all(rownames(bVals_mat_common) %in% rownames(bVals_menninges_common))
all(rownames(bVals_mat_common) == rownames(bVals_menninges_common))

bVals_combine = cbind(bVals_menninges_common,bVals_mat_common)

###########################check average METHhigh signature activity in meninges and MNG, avg per sample
#get METH high probe IDs for EPIC
hyper_probes_EPIC = read.csv(file="c2_hyper_probes_EPIC.csv", header = T)
#subset combined bVal data (is only 4941 probes left in combined bVal data)
bVals_combine_hyper_probes = bVals_combine %>% dplyr::filter(rownames(bVals_combine) %in% hyper_probes_EPIC$ID)
#calculate avg per sample
sample_avg_hyper_probes = apply(bVals_combine_hyper_probes,2,mean)
#make dataframe for plot
df_sample_signature_hyper = as.data.frame(sample_avg_hyper_probes)
df_sample_signature_hyper$setting_detail = targets_combine$setting_detail
write.csv(df_sample_signature_hyper, file="avg_hyper_signature.csv")

cairo_pdf(filename = "Boxplot_hyper_probe_avg_combined_data.pdf", width = 8, height = 8)
ggplot(df_sample_signature_hyper, 
       aes(setting_detail, sample_avg_hyper_probes))+
  geom_point(position = position_jitter(width=0.1), alpha=0.95, aes(color=sample_avg_hyper_probes), size=6)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  scale_y_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.2,0.82))+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red3"), limits=c(0,1))+
  theme_classic()+
  theme(axis.line=element_line(size=1))+
  theme(legend.position = "none")
dev.off()


##########check development of hyper signature in several samples from recurrence patients
targets_rec = targets_MNG[-c(1:13),]
write.csv(targets_rec, file="targets_rec.csv")
targets_rec = read.csv(file="targets_rec.csv")

cairo_pdf(filename = "Lineplot_hyper_signature_development_patient.pdf", width = 5, height = 8)
ggplot(targets_rec, aes(x=setting_detail, y = hyper_avg, group = patient))+
  geom_line(alpha=0.8)+
  geom_point(aes(color=hyper_avg),size=5)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red3"), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

###########################check average METH cluster signature activity in meninges and MNG, avg per sample
#re-do as above, but use cluster signature (top2000) from here on, not sig hyper probes
#get EPICv1 probes names for top2000 probes
EPICv1_v2 = data.frame(EPICv2 = rep(NA_character_, 930075))
EPICv1_v2$EPICv2 = rownames(annoEPICv2)
EPICv1_v2$EPICv1 = annoEPICv2$EPICv1_Loci
write.csv(EPICv1_v2, file="EPICv1_v2.csv")

#get probes from top2000 cluster on EPICv1
METH_cluster_top2000_EPICv1 = read.csv(file="METH_cluster_top2000_EPICv1.csv", header = T)

#subset combined bVal data by top2000 probes
bVals_combine_top2000 = bVals_combine %>% dplyr::filter(rownames(bVals_combine) %in% METH_cluster_top2000_EPICv1$ID)
#calculate avg per sample
sample_avg_top2000 = apply(bVals_combine_top2000,2,mean)
#make dataframe for plot
df_sample_signature_cluster = as.data.frame(sample_avg_top2000)
df_sample_signature_cluster$setting_detail = targets_combine$setting_detail
write.csv(df_sample_signature_cluster, file="avg_cluster_signature.csv")

library(ggplot2)
cairo_pdf(filename = "Boxplot_cluster_probes_avg_combined_data.pdf", width = 8, height = 6)
ggplot(df_sample_signature_cluster, 
       aes(setting_detail, sample_avg_top2000))+
  geom_point(position = position_jitter(width=0.1), alpha=0.95, aes(color=sample_avg_top2000), size=6)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  scale_y_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.85))+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red3"), limits=c(0,1))+
  theme_classic()+
  theme(axis.line=element_line(size=0.5))+
  theme(legend.position = "none")
dev.off()


##########check development of cluster signature in several samples from recurrence patients
targets_rec = read.csv(file="targets_rec.csv")
str(targets_rec)
targets_rec$WHO_2021 = as.numeric(targets_rec$WHO_2021)

cairo_pdf(filename = "Lineplot_cluster_signature_development_patient.pdf", width = 5, height = 6)
ggplot(targets_rec, aes(x=setting_detail, y = cluster_avg, group = patient))+
  geom_line(alpha=0.8)+
  geom_point(aes(color=cluster_avg),size=6)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red3"), limits=c(0,1))+
  theme_classic()+
  theme(axis.line=element_line(size=0.5))+
  theme(legend.position = "none")
dev.off()

###########check differential methylation, only include probes that are setting (meninges, Pri_no, Pri_yes, Recurrence) specific, then do PCA
###first make new meta file
colnames(targets_validation)[colnames(targets_validation) == "MCConsensus"] = "MCconsensus"
colnames(targets_meninges)[colnames(targets_meninges) == "Subtype"] = "MCsubtype"
colnames(targets_meninges)[colnames(targets_meninges) == "Mcconsensus"] = "MCconsensus"
targets_combine = rbind(targets_meninges, targets_validation)

#make summarized experiment
library(SummarizedExperiment)
se.combine <- SummarizedExperiment(as.matrix(bVals_combine), colData = targets_combine)
str(se.combine)
dim(se.combine)
colData(se.combine)$setting
assay(se.combine)
colData(se.combine)$setting = relevel(factor(colData(se.combine)$setting),"Meninges")
colData(se.combine)$setting
colData(se.combine)$patient = relevel(factor(colData(se.combine)$patient),"1")
colData(se.combine)$patient

#check for NAs
se_ok = (checkLevels(assay(se.combine), colData(se.combine)$setting) &
           checkLevels(assay(se.combine), colData(se.combine)$patient))
sum(se_ok)
#no NAs in the data

#differential methylation for setting contrast all
smry_combine = DML(se.combine, ~setting)
smry_combine
test_result_combine = summaryExtractTest(smry_combine)
colnames(test_result_combine)

#differential methylation for only contrast primary_no_rec, primary_rec, and recurrences (MNG probes)
#generate probes specific for only MNG, without meninges
#make subsets of bVals and targets first
bVals_MNG = bVals_combine[,-c(1:20)]
targets_MNG = targets_combine[-c(1:20),]
#make summarized experiment
se.MNG <- SummarizedExperiment(as.matrix(bVals_MNG), colData = targets_MNG)
str(se.MNG)
colData(se.MNG)$setting = relevel(factor(colData(se.MNG)$setting),"Primary_no_rec")
se_MNG_ok = (checkLevels(assay(se.MNG), colData(se.MNG)$setting))
sum(se_MNG_ok)
smry_MNG = DML(se.MNG, ~setting)
test_result_MNG = summaryExtractTest(smry_MNG)
colnames(test_result_MNG)


#check setting specific probes (112,864)
test_result_setting = test_result_combine %>%
  mutate(setting_specific =
           ifelse(FPval_setting < 0.01 & Eff_setting > 0.1, TRUE, FALSE)) 
setting_probes_combine = test_result_setting %>% dplyr::filter(setting_specific=="TRUE")
dim(setting_probes_combine)


#check MNG specific probes (17,161)
colnames(test_result_MNG)
test_result_MNG = test_result_MNG %>%
  mutate(setting_specific_MNG =
           ifelse(FPval_setting < 0.01 & Eff_setting > 0.1, TRUE, FALSE))
setting_MNG_probes = test_result_MNG %>% dplyr::filter(setting_specific_MNG=="TRUE")
dim(setting_MNG_probes)

#####subset
#make beta values subset of setting and setting_detail specific probes, and do PCA
bVals_combine_setting_probes = bVals_combine %>% filter(rownames(bVals_combine) %in% (setting_probes_combine$Probe_ID))
dim(bVals_combine_setting_probes)

bVals_combine_setting_probes_MNG = bVals_combine %>% filter(rownames(bVals_combine) %in% (setting_MNG_probes$Probe_ID))
dim(bVals_combine_setting_probes_MNG)

#make PCA
bVals_combine_setting_PCA = prcomp(t(bVals_combine_setting_probes), scale. = T)
bVals_combine_setting_MNG_PCA = prcomp(t(bVals_combine_setting_probes_MNG), scale. = T)

df_PCA_setting_all <- cbind(targets_combine, bVals_combine_setting_PCA$x[,1:3])
df_PCA_setting_MNG <- cbind(targets_combine, bVals_combine_setting_MNG_PCA$x[,1:3])

save.image()

#scree plots
library(factoextra)

tiff("Scree_PCs_setting_all_probes_combined_cohort_labels.tiff", res=300,width=1260,height=2048)
fviz_screeplot(bVals_combine_setting_PCA, addlabels=T)
dev.off()

tiff("Scree_PCs_setting_all_probes_combined_cohort_clean.tiff", res=300,width=1260,height=2048)
fviz_screeplot(bVals_combine_setting_PCA, addlabels=F)
dev.off()

tiff("Scree_PCs_setting_MNG_probes_combined_cohort_labels.tiff", res=300,width=1260,height=2048)
fviz_screeplot(bVals_combine_setting_MNG_PCA, addlabels=T)
dev.off()

tiff("Scree_PCs_setting_MNG_probes_combined_cohort_clean.tiff", res=300,width=1260,height=2048)
fviz_screeplot(bVals_combine_setting_MNG_PCA, addlabels=F)
dev.off()

#perform k-means clustering to add clusters in PCA plot
save.image()
#optimal number of clusters k for all setting probes
tiff("Elbow_clusters_combine_probes.tiff", res=300,width=1548,height=1548)
fviz_nbclust(bVals_combine_setting_probes, kmeans,  method = "wss")+
  geom_vline(xintercept = 2, linetype=2)
dev.off()

#perform kmeans clustering for all setting probes
k=2
km_res = kmeans(t(bVals_combine_setting_probes), centers = k, nstart = 25)
df_PCA_setting_all$kmeans = as.factor(km_res$cluster)

#make PCA plot with kmeans for all setting probes
cairo_pdf(filename = "PCA_combine_setting_detail_kmeans_cluster_all_probes_2.pdf",width = 8, height = 6)
ggplot(df_PCA_setting_all, aes(x=PC1, y=PC2, linetype=kmeans))+
  geom_point(size=6, aes(color=setting_detail))+
  scale_color_manual(values=c(Dura ="royalblue3",
                              Leptomeninges="mediumpurple4",
                              Primary_no_rec = "#1CFDB2",
                              Primary_rec = "darkgreen",
                              Rec_1="#FFCB1B",
                              Rec_2="#D95F02"))+
  stat_ellipse(type="t", level = 0.95)+
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.8))+
  theme(legend.position = "none")
dev.off()

theme(axis.line = element_line(linewidth=0.5))+
  #optimal number of clusters k for MNG setting probes
  tiff("Elbow_clusters_combine_MNG_probes.tiff", res=300,width=1548,height=1548)
fviz_nbclust(bVals_combine_setting_probes_MNG, kmeans,  method = "wss")+
  geom_vline(xintercept = 2, linetype=2)
dev.off()

#perform kmeans clustering for MNG probes
k=2
km_res = kmeans(t(bVals_combine_setting_probes_MNG), centers = k, nstart = 25)
df_PCA_setting_MNG$kmeans = as.factor(km_res$cluster)

#make PCA plot with kmeans for MNG probes
cairo_pdf(filename = "PCA_combine_setting_detail_kmeans_cluster_MNG_probes_2.pdf", width = 8, height = 6)
ggplot(df_PCA_setting_MNG, aes(x=PC1, y=PC2, linetype=kmeans))+
  geom_point(size=6, aes(color=setting_detail))+
  scale_color_manual(values=c(Dura ="royalblue3",
                              Leptomeninges="mediumpurple4",
                              Primary_no_rec = "#1CFDB2",
                              Primary_rec = "darkgreen",
                              Rec_1="#FFCB1B",
                              Rec_2="#D95F02"))+
  stat_ellipse(type="t", level = 0.95)+
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.8))+
  theme(legend.position = "none")
dev.off()

########make differential methylation to test differences between pri_no_rec, pri_rec, and rec as compared to meninges
###overlay METH cluster signature and common MNG probes

#differential methylation for kmeans cluster for all probes including pvalue adjustment
smry_meninges_contrast = DML(se_kmeans, ~setting)
test_result_meninges_contrast = summaryExtractTest(smry_meninges_contrast)
colnames(test_result_meninges_contrast)

padj_Pri_no_rec = p.adjust(test_result_meninges_contrast$Pval_settingPrimary_no_rec, method = "BH", n = length(test_result_meninges_contrast$Pval_settingPrimary_no_rec))
padj_Pri_rec = p.adjust(test_result_meninges_contrast$Pval_settingPrimary_rec, method = "BH", n = length(test_result_meninges_contrast$Pval_settingPrimary_rec))
padj_Rec = p.adjust(test_result_meninges_contrast$Pval_settingRecurrence, method = "BH", n = length(test_result_meninges_contrast$Pval_settingRecurrence))
test_result_meninges_contrast = cbind(test_result_meninges_contrast, padj_Pri_no_rec, padj_Pri_rec, padj_Rec)

save.image()

#determine which probes are shared across all 3 contrasts
test_result_meninges_contrast = test_result_meninges_contrast %>%
  mutate(Pri_no_rec_specific =
           ifelse(padj_Pri_no_rec < 0.05 & abs(Est_settingPrimary_no_rec) > 0.2, TRUE, FALSE),
         Pri_rec_specific =
           ifelse(padj_Pri_rec < 0.05 & abs(Est_settingPrimary_rec) > 0.2, TRUE, FALSE),
         Rec_specific =
           ifelse(padj_Rec < 0.05 & abs(Est_settingRecurrence) > 0.2, TRUE, FALSE)
  ) 

Pri_no_rec_probes = test_result_meninges_contrast %>% dplyr::filter(Pri_no_rec_specific=="TRUE")
Pri_rec_probes = test_result_meninges_contrast %>% dplyr::filter(Pri_rec_specific=="TRUE")
Rec_probes = test_result_meninges_contrast %>% dplyr::filter(Rec_specific=="TRUE")
str(Pri_no_rec_probes)
meninges_contrast_common = Reduce(intersect, list(Pri_no_rec_probes$Probe_ID, Pri_rec_probes$Probe_ID, Rec_probes$Probe_ID))


#make plots
#show contrast Primary_no_rec vs meninges
cairo_pdf(filename = "Volcano_meth_Pri_no_rec_meninges_common_cluster.pdf", width = 4, height = 6)
test_result_meninges_contrast %>%
  ggplot(aes(x=Est_settingPrimary_no_rec, y=-log10(padj_Pri_no_rec), size=abs(Est_settingPrimary_no_rec))) + 
  geom_point(alpha=0.5, shape=16)+
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% meninges_contrast_common), 
             aes(x=Est_settingPrimary_no_rec, y=-log10(padj_Pri_no_rec)), color="slateblue3", alpha=0.6, shape=16)+
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% METH_cluster_top2000_EPICv1$ID), 
             aes(x=Est_settingPrimary_no_rec, y=-log10(padj_Pri_no_rec)), color="red3", alpha=0.7, shape=16)+
  scale_radius(range = c(0.05,5))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.6,0.25), limits=c(-0.5,0.6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.8))+
  theme(legend.position = "none")
dev.off()

#show contrast Primary_rec vs meninges
cairo_pdf(filename = "Volcano_meth_Pri_rec_meninges_common_cluster.pdf", width = 4, height = 6)
test_result_meninges_contrast %>%
  ggplot(aes(x=Est_settingPrimary_rec, y=-log10(padj_Pri_rec), size=abs(Est_settingPrimary_rec))) + 
  geom_point(alpha=0.5, shape=16)+
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% meninges_contrast_common), 
             aes(x=Est_settingPrimary_rec, y=-log10(padj_Pri_rec)), color="slateblue3", alpha=0.6, shape=16)+
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% METH_cluster_top2000_EPICv1$ID), 
             aes(x=Est_settingPrimary_rec, y=-log10(padj_Pri_rec)), color="red3", alpha=0.7, shape=16)+
  scale_radius(range = c(0.05,5))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.6,0.25), limits=c(-0.5,0.7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.8))+
  theme(legend.position = "none")
dev.off()

#show contrast Recurrence vs meninges
cairo_pdf(filename = "Volcano_meth_Recurrence_meninges_common_cluster.pdf", width = 4, height = 6)
test_result_meninges_contrast %>%
  ggplot(aes(x=Est_settingRecurrence, y=-log10(padj_Rec), size=abs(Est_settingRecurrence))) + 
  geom_point(alpha=0.5, shape=16)+
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% meninges_contrast_common), 
             aes(x=Est_settingRecurrence, y=-log10(padj_Rec)), color="slateblue3", alpha=0.6, shape=16)+
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% METH_cluster_top2000_EPICv1$ID), 
             aes(x=Est_settingRecurrence, y=-log10(padj_Rec)), color="red3", alpha=0.7, shape=16)+
  scale_radius(range = c(0.05,5))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.6,0.25), limits=c(-0.5,0.7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.8))+
  theme(legend.position = "none")
dev.off()


##############perform slingshot for probes associated with conditions primary_no, primary_yes, and recurrence
####bVals_setting_probes
#####and read in as normal counts and norm counts
##subset all data by MNG probes and make matrix
bVals_combine_MNG_probes = bVals_combine %>% filter(rownames(bVals_combine) %in% (setting_MNG_probes$Probe_ID))
bVals_combine_sling = as.matrix(bVals_combine_MNG_probes)

#make singlecellexperiment
sce_combine <- SingleCellExperiment(assays = List(counts = bVals_combine_sling))
assays(sce_combine)$norm <- bVals_combine_sling

#perform PCA and add values to singlecellexperiment
pca_combine_sling <- prcomp(t(assays(sce_combine)$norm), scale. = TRUE)
rd1_combine <- pca_combine_sling$x[,1:2]

plot(rd1_combine, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(sce_combine) <- SimpleList(PCA = rd1_combine)

#add information on clusters
targets_combine_sling = read.csv(file="targets_combine_sling.csv", header = T)
all(targets_combine_sling$ID == sce_combine@colData@rownames)

cl1_combine = targets_combine_sling$cluster_development
names(cl1_combine) = targets_combine_sling$ID
cl1_combine
all(colData(sce_combine)$rownames == names(cl1_combine))
colData(sce_combine)$GMM <- cl1_combine

#check clusters on PCA
library(RColorBrewer)
plot(rd1_combine, col = brewer.pal(9,"Set1")[cl1_combine], pch=16, asp = 1)

#perform sling
sce_combine <- slingshot(sce_combine, clusterLabels = 'GMM', reducedDim = 'PCA')

summary(sce_combine$slingPseudotime_1)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00   69.02  318.03  250.67  377.53  467.46      16 

#show trajectory in plot
library(grDevices)

plot(reducedDims(sce_combine)$PCA, pch=16, asp = 1)
lines(SlingshotDataSet(sce_combine), lwd=2, col='black')


#define linages by start and end point
lin1 <- getLineages(rd1_combine, cl1_combine, start.clus="1", end.clus="3")
lin1
#check minimum span
plot(rd1_combine, col = brewer.pal(9,"Set1")[cl1_combine], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

#get ciurve
crv1 <- getCurves(lin1,extend="pc1")
crv1

plot(rd1_combine, col = brewer.pal(9,"Set1")[cl1_combine], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')

slingPseudotime(crv1)

save.image()

###use ggplot to visualize the trajectory inference
#make dataframe
df_combine = data.frame(rd1_combine, "cl1_combine" = as.character(cl1_combine))
all(rownames(df_combine) == targets_combine_sling$ID)
all(names(crv1)==rownames(df_combine))

#make curve data a dataframe
crv1_plot = SlingshotDataSet(crv1)
crv1_data = as.data.frame(crv1_plot@curves$Lineage1$s)
crv1_data$order = crv1_plot@curves$Lineage1$ord
crv1_data$lineage = as.integer(crv1_plot@curves$Lineage1$w)
crv1_data = crv1_data[order(crv1_data$order, decreasing=FALSE),]

all(rownames(df_combine)==rownames(crv1_data))
df_combine= df_combine[match(rownames(crv1_data),rownames(df_combine)),]
all(rownames(df_combine)==rownames(crv1_data))
df_combine$curve_x = crv1_data$PC1
df_combine$curve_y = crv1_data$PC2
df_combine$order = crv1_data$order
df_combine$lineage = crv1_data$lineage

#add methylation signature and setting information
all(rownames(df_combine)==targets_combine_sling$ID)
df_combine= df_combine[match(targets_combine_sling$ID,rownames(df_combine)),]
all(rownames(df_combine)==targets_combine_sling$ID)
avg_top2000_signature_combine = read.csv(file="avg_cluster_signature.csv", header = T)
rownames(avg_top2000_signature_combine) = avg_top2000_signature_combine$X
all(rownames(df_combine) %in% rownames(avg_top2000_signature_combine))
all(rownames(df_combine) == rownames(avg_top2000_signature_combine))
df_combine= df_combine[match(rownames(avg_top2000_signature_combine),rownames(df_combine)),]
all(rownames(df_combine) == rownames(avg_top2000_signature_combine))
df_combine$avg_signature = avg_top2000_signature_combine$sample_avg_top2000
df_combine$MCconsensus = avg_top2000_signature_combine$MCconsensus
df_combine$setting_detail = avg_top2000_signature_combine$setting_detail
df_combine$grade = avg_top2000_signature_combine$grade
df_combine$grade = factor(df_combine$grade)
str(df_combine)

#smooth trajectory curve
smoothed_curve = smooth.spline(x=df_combine$curve_x, y=df_combine$curve_y, spar = 0.95)
df_smoothed_curve = as.data.frame(cbind(smoothed_curve$x, smoothed_curve$y))

#make the plot

#for METH activity
library(ggplot2)
cairo_pdf(filename = "PCA_Trajectory_combine_METH_activity_with_curve.pdf", width = 8, height = 8)
ggplot(df_combine, aes(x = PC1, y = PC2, color=avg_signature)) +
  geom_point(show.legend = F, size=6)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red"))+
  geom_line(data=df_smoothed_curve, aes(V1, V2),color="black", size=1.2)+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  coord_fixed(ratio=1)+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()

cairo_pdf(filename = "PCA_Trajectory_combine_METH_activity_without_curve.pdf", width = 8, height = 8)
ggplot(df_combine, aes(x = PC1, y = PC2, color=avg_signature)) +
  geom_point(show.legend = F, size=6)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  coord_fixed(ratio=1)+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()

cairo_pdf(filename = "PCA_Trajectory_combine_METH_activity_with_legend.pdf", width = 8, height = 8)
ggplot(df_combine, aes(x = PC1, y = PC2, color=avg_signature)) +
  geom_point(size=6)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  coord_fixed(ratio=1)+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()

#meningioma settings

cairo_pdf(filename = "PCA_Trajectory_combine_MNG_settings_with_curve.pdf", width = 8, height = 8)
ggplot(df_combine, aes(x = PC1, y = PC2, color=setting_detail)) +
  geom_point(show.legend = F, size=6)+
  scale_color_manual(values=c(Dura ="royalblue3",
                              Leptomeninges="mediumpurple4",
                              Primary_no_rec = "#1CFDB2",
                              Primary_rec = "darkgreen",
                              Rec_1="#FFCB1B",
                              Rec_2="#D95F02"))+
  geom_line(data=df_smoothed_curve, aes(V1, V2),color="black", size=1.2)+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  coord_fixed(ratio=1)+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()


#molecular group

cairo_pdf(filename = "PCA_Trajectory_combine_group_with_curve.pdf", width = 8, height = 8)
ggplot(df_combine, aes(x = PC1, y = PC2, color=MCconsensus)) +
  geom_point(show.legend = F, size=6)+
  scale_color_manual(values=c(control ="grey60",
                              Merlinintact="royalblue2",
                              hypermetabolic = "forestgreen",
                              proliferative = "darkorange2",
                              Immuneenriched="red3"))+
  geom_line(data=df_smoothed_curve, aes(V1, V2),color="black", size=1.2)+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  coord_fixed(ratio=1)+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()

#WHO grading

cairo_pdf(filename = "PCA_Trajectory_combine_grade_with_curve.pdf", width = 8, height = 8)
ggplot(df_combine, aes(x = PC1, y = PC2, color=grade)) +
  geom_point(show.legend = F, size=6)+
  scale_color_manual(values=c("grey60", "#9EBCDA","#8C6BB1", "#810F7C"))+
  geom_line(data=df_smoothed_curve, aes(V1, V2),color="black", size=1.2)+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  coord_fixed(ratio=1)+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()


####make pseudotime plots versus setting and METH activity
#get data for psdeudotime/setting
pseudotime_combine = as.data.frame(slingPseudotime(crv1))
all(rownames(pseudotime_combine)==rownames(df_combine))
pseudotime_combine$pseudotime = pseudotime_combine$Lineage1
pseudotime_combine$avg_signature = avg_top2000_signature_combine$sample_avg_top2000
pseudotime_combine$METH_cluster = targets_combine$Cluster
pseudotime_combine$METH_cluster = factor(pseudotime_combine$METH_cluster, levels = c("1","2","3"))
pseudotime_combine$setting_detail = targets_combine$setting_detail


cairo_pdf(filename = "METH_activity_along_pseudotime.pdf", width = 7, height = 8)
ggplot(pseudotime_combine, aes(x = pseudotime, y = avg_signature, color=avg_signature)) +
  geom_point(show.legend = F, size=8)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()

cor.test(pseudotime_combine$pseudotime, pseudotime_combine$avg_signature)
#t = 14.444, df = 73, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.7875820 0.9099097
#sample estimates:
# cor 
#0.8606939 

cairo_pdf(filename = "Setting_detail_along_pseudotime.pdf", width = 7, height = 8)
ggplot(pseudotime_combine, aes(x = pseudotime, y = setting_detail, color=setting_detail)) +
  geom_point(show.legend = F, size=8)+
  scale_color_manual(values=c(Dura ="royalblue3",
                              Leptomeninges="mediumpurple4",
                              Primary_no_rec = "#1CFDB2",
                              Primary_rec = "darkgreen",
                              Rec_1="#FFCB1B",
                              Rec_2="#D95F02"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=1)
  )+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()


#hide label meninges before plotting
pseudotime_combine_subset = pseudotime_combine[-c(1:20),]
pseudotime_combine_subset = pseudotime_combine_subset[order(pseudotime_combine_subset$pseudotime, decreasing = T),]

cairo_pdf(filename = "METH_cluster_along_pseudotime.pdf", width = 7, height = 8)
ggplot(pseudotime_combine_subset, aes(x = pseudotime, y = METH_cluster, color=METH_cluster)) +
  geom_point(show.legend = F, size=8)+
  scale_color_manual(values=c("1" = "cyan4",
                              "2" = "violetred4"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=1)
  )+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )+
  scale_x_continuous(limits = c(0, 410))
dev.off()















