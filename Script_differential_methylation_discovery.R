library(sesame)
library(ggplot2)
library(BiocParallel)
library(parallel)
library(rockchalk)
getwd()

#path to idats
idat_dir = "Directory"

#get meta data
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_differential_methylation/meta")
sampledata = read.csv(file = "metadata.csv", header = T, row.names = 1)
sampledata = sampledata[,-1]

#employ sesame pipeline to get beta values
beta.discovery = openSesame(idat_dir)

#re-order beta dataframe columns (Basename) to match Basename order in sampledata
beta.discovery = beta.discovery[,match(sampledata$Basename,colnames(beta.discovery))]
all(sampledata$Basename %in% colnames(beta.discovery))
all(sampledata$Basename == colnames(beta.discovery))
beta.discovery = as.data.frame(beta.discovery)
colnames(beta.discovery) = rownames(sampledata)

#remove sex chromosome probes
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

keep <- !(rownames(beta.discovery) %in% annoEPICv2$Name[annoEPICv2$chr %in% 
                                                          c("chrX","chrY")])
table(keep)
beta.discovery = beta.discovery[keep,]

#remove all probes which are masked in at least one case
beta.discovery.flt = beta.discovery[complete.cases(beta.discovery),]


# differential methylation in discovery cohort
#first, make summarizedExperiment
library(SummarizedExperiment)
se.discovery <- SummarizedExperiment(as.matrix(beta.discovery.flt), colData = sampledata)
str(se.discovery)
dim(se.discovery)
colData(se.discovery)
assay(se.discovery)

#subset summarizedexperiment for individual contrasts
se.discovery_for_grade_2 = se.discovery[,!se.discovery$grading2016 == "3"]
str(se.discovery_for_grade_2)
colData(se.discovery_for_grade_2)$grading2016 <- relevel(factor(colData(se.discovery_for_grade_2)$grading2016), "1")
str(se.discovery_for_grade_2)

se.discovery_for_grade_3 = se.discovery[,!se.discovery$grading2016 == "2"]
str(se.discovery_for_grade_3)
colData(se.discovery_for_grade_3)$grading2016 <- relevel(factor(colData(se.discovery_for_grade_3)$grading2016), "1")
str(se.discovery_for_grade_3)

se.discovery_for_hypermetabolic = se.discovery[,!se.discovery$MCconsensus == "proliferative"]
colData(se.discovery_for_hypermetabolic)$MCconsensus <- relevel(factor(colData(se.discovery_for_hypermetabolic)$MCconsensus), "Merlinintact")
colData(se.discovery_for_hypermetabolic)$MCconsensus = combineLevels(se.discovery_for_hypermetabolic$MCconsensus,
                                                                     levs = c("Merlinintact", "Immuneenriched"), newLabel = "other")
colData(se.discovery_for_hypermetabolic)$MCconsensus <- relevel(factor(colData(se.discovery_for_hypermetabolic)$MCconsensus), "other")

se.discovery_for_proliferative = se.discovery[,!se.discovery$MCconsensus == "hypermetabolic"]
colData(se.discovery_for_proliferative)$MCconsensus <- relevel(factor(colData(se.discovery_for_proliferative)$MCconsensus), "Merlinintact")
colData(se.discovery_for_proliferative)$MCconsensus = combineLevels(se.discovery_for_proliferative$MCconsensus,
                                                                    levs = c("Merlinintact", "Immuneenriched"), newLabel = "other")
colData(se.discovery_for_proliferative)$MCconsensus <- relevel(factor(colData(se.discovery_for_proliferative)$MCconsensus), "other")


####determine differentially methylated probes for several predictors of interest (cluster, grading,status, invasion)
#####prepare predictor variables
#turn discrete contrast variables to factor and relevel
str(se.discovery$status)
colData(se.discovery)$status <- relevel(factor(colData(se.discovery)$status), "Primary")
str(se.discovery$status)
str(se.discovery$invasionhisto)
colData(se.discovery)$invasionhisto <- relevel(factor(colData(se.discovery)$invasionhisto), "no")
str(se.discovery$invasionhisto)
str(se.discovery$cluster_2000)
colData(se.discovery)$cluster_2000 <- relevel(factor(colData(se.discovery)$cluster_2000), "1")
str(se.discovery$cluster_2000)
colData(se.discovery)$MCgroup <- relevel(factor(colData(se.discovery)$MCgroup), "Merlin-intact")
str(se.discovery$MCgroup)
colData(se.discovery)$MCgroup2 = combineLevels(se.discovery$MCgroup, levs = c("Merlin-intact", "Immune-enriched"), newLabel = "other")
colData(se.discovery)$MCgroup2 <- relevel(factor(colData(se.discovery)$MCgroup2), "other")
str(se.discovery$MCgroup2)
str(se.discovery$MCsubtype)
colData(se.discovery)$MCsubtype <- relevel(factor(colData(se.discovery)$MCsubtype), "Benign")
str(se.discovery$MCsubtype)
colData(se.discovery)$MCsubtype2 = combineLevels(se.discovery$MCsubtype, levs = c("Intermediate", "Malignant"), newLabel = "Int/Mal")
str(se.discovery$MCsubtype2)
str(se.discovery$grading2016)
colData(se.discovery)$grading2016 <- relevel(factor(colData(se.discovery)$grading2016), "1")
str(se.discovery$grading2016)
colData(se.discovery)$grading = combineLevels(se.discovery$grading2016, levs = c("2", "3"), newLabel = "highgrade")
str(se.discovery$grading)
str(se.discovery)



###model methylation differences

####for cluster_2000
smry.discovery.cluster = DML(se.discovery, ~cluster_2000)
lm_results_discovery_cluster = summaryExtractTest(smry.discovery.cluster)
# adjust p values for status
padj.c2vsc1 = p.adjust(lm_results_discovery_cluster$Pval_cluster_20002, method = "BH", n = length(lm_results_discovery_cluster$Pval_cluster_20002))
lm_results_discovery_cluster = cbind(lm_results_discovery_cluster, padj.c2vsc1)

#visualize differences in methylation in volcano for cluster
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_differential_methylation/results")
cairo_pdf(filename = "Volcano_c2_vs_c1_probe_level.pdf", width = 8, height = 6)
ggplot(lm_results_discovery_cluster) + geom_point(aes(x=Est_cluster_20002, y=-log10(padj.c2vsc1), color=-log10(padj.c2vsc1), size=abs(Est_cluster_20002)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.05, 0)) + 
  scale_radius(range = c(0.05,5))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.5,0.25), limits=c(-0.5,0.5))+
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -0.2, yend = Inf), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = 0.2, yend = Inf), colour = "grey35", linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in c2 vs c1
#204 hypomethylated, 6,812 hypermethylated in c2
hypo_c2_vs_c1 = lm_results_discovery_cluster %>% dplyr::filter(Est_cluster_20002 < -0.2 & padj.c2vsc1 < 0.05)
hyper_c2_vs_c1 = lm_results_discovery_cluster %>% dplyr::filter(Est_cluster_20002 > 0.2 & padj.c2vsc1 < 0.05)





####model for invasion
smry.discovery.invasion = DML(se.discovery, ~invasionhisto)
lm_results_discovery_invasion = summaryExtractTest(smry.discovery.invasion)
str(lm_results_discovery_invasion)

# adjust p values for invasion
padj.invasion = p.adjust(lm_results_discovery_invasion$Pval_invasionhistoyes, method = "BH", n = length(lm_results_discovery_invasion$Pval_invasionhistoyes))
lm_results_discovery_invasion = cbind(lm_results_discovery_invasion, padj.invasion)

#visualize differences in methylation in volcano for invasion
cairo_pdf(filename = "Volcano_Invasion_probe_level.pdf", width = 8, height = 6)
ggplot(lm_results_discovery_invasion) + geom_point(aes(x=Est_invasionhistoyes, y=-log10(padj.invasion), color=-log10(padj.invasion), size=abs(Est_invasionhistoyes)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.05, 0)) + 
  scale_size(range = c(0.3,3))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.5,0.25), limits=c(-0.5,0.5))+
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -0.2, yend = Inf), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = 0.2, yend = Inf), colour = "grey35", linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()



#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in cases with invasion
# 17 hypomethylated,  52 hypermethylated probes in cases with invasion
hypo_invasion = lm_results_discovery_invasion %>% dplyr::filter(Est_invasionhistoyes < -0.2 & padj.invasion < 0.05)
hyper_invasion = lm_results_discovery_invasion %>% dplyr::filter(Est_invasionhistoyes > 0.2 & padj.invasion < 0.05)





####model for MCgroup2
smry.discovery.group = DML(se.discovery, ~MCgroup2)
lm_results_discovery_group = summaryExtractTest(smry.discovery.group)
str(lm_results_discovery_group)

# adjust p values for hypermitotic vs Merlin-intact/Immune-enriched
padj.hypermitotic = p.adjust(lm_results_discovery_group$Pval_MCgroup2Hypermitotic, method = "BH", n = length(lm_results_discovery_group$Pval_MCgroup2Hypermitotic))
lm_results_discovery_group = cbind(lm_results_discovery_group, padj.hypermitotic)

#visualize hypermitotic
cairo_pdf(filename = "Volcano_hypermitoticvsreference_probe_level.pdf", width = 8, height = 6)
ggplot(lm_results_discovery_group) + geom_point(aes(x=Est_MCgroup2Hypermitotic, y=-log10(padj.hypermitotic), color=-log10(padj.hypermitotic), size=abs(Est_MCgroup2Hypermitotic)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.05, 0)) + 
  scale_radius(range = c(0.5,5))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.5,0.25), limits=c(-0.5,0.5))+
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -0.2, yend = Inf), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = 0.2, yend = Inf), colour = "grey35", linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in hypermitotic over reference
# 262 hypomethylated,  3470 hypermethylated probes in hypermitotic vs reference
hypo_hypermitotic = lm_results_discovery_group %>% dplyr::filter(Est_MCgroup2Hypermitotic < -0.2 & padj.hypermitotic < 0.05)
hyper_hypermitotic = lm_results_discovery_group %>% dplyr::filter(Est_MCgroup2Hypermitotic > 0.2 & padj.hypermitotic < 0.05)



####model for grading highgrade (2+3) vs 1
smry.discovery.grade = DML(se.discovery, ~grading)
lm_results_discovery_grade = summaryExtractTest(smry.discovery.grade)
str(lm_results_discovery_grade)

# adjust p values for high grade vs grade1
padj.gradehigh = p.adjust(lm_results_discovery_grade$Pval_gradinghighgrade, method = "BH", n = length(lm_results_discovery_grade$Pval_gradinghighgrade))
lm_results_discovery_grade = cbind(lm_results_discovery_grade, padj.gradehigh)

#visualize grade3 associated probes
cairo_pdf(filename = "Volcano_gradehighvsreference_probe_level.pdf", width = 8, height = 6)
ggplot(lm_results_discovery_grade) + geom_point(aes(x=Est_gradinghighgrade, y=-log10(padj.gradehigh), color=-log10(padj.gradehigh), size=abs(Est_gradinghighgrade)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.05, 0)) + 
  scale_size(range = c(0.3,3))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.5,0.25), limits=c(-0.5,0.5))+
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -0.2, yend = Inf), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = 0.2, yend = Inf), colour = "grey35", linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in high grade vs grade1
# 60 hypomethylated,  705 hypermethylated probes in in high grade vs grade1
hypo_gradehigh = lm_results_discovery_grade %>% dplyr::filter(Est_gradinghighgrade < -0.2 & padj.gradehigh < 0.05)
hyper_gradehigh = lm_results_discovery_grade %>% dplyr::filter(Est_gradinghighgrade > 0.2 & padj.gradehigh < 0.05)





####model for grading grade 3 vs 1+2
smry.discovery.grade2 = DML(se.discovery, ~grading)
lm_results_discovery_grade2 = summaryExtractTest(smry.discovery.grade2)
str(lm_results_discovery_grade2)

# adjust p values for grade3 vs reference
padj.grade3 = p.adjust(lm_results_discovery_grade2$Pval_grading23, method = "BH", n = length(lm_results_discovery_grade2$Pval_grading23))
lm_results_discovery_grade2 = cbind(lm_results_discovery_grade2, padj.grade3)

#visualize grade3 associated probes
cairo_pdf(filename = "Volcano_grade3vs12_probe_level.pdf", width = 8, height = 6)
ggplot(lm_results_discovery_grade2) + geom_point(aes(x=Est_grading23, y=-log10(padj.grade3), color=-log10(padj.grade3), size=abs(Est_grading23)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.05, 0)) + 
  scale_size(range = c(0.3,3))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.5,0.25), limits=c(-0.5,0.5))+
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -0.2, yend = Inf), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = 0.2, yend = Inf), colour = "grey35", linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade3 over reference
# 6939 hypomethylated,  18350 hypermethylated probes in grade3 vs reference
hypo_grade3 = lm_results_discovery_grade2 %>% dplyr::filter(Est_grading23 < -0.2 & padj.grade3 < 0.05)
hyper_grade3 = lm_results_discovery_grade2 %>% dplyr::filter(Est_grading23 > 0.2 & padj.grade3 < 0.05)


#model differences grade2 vs grade1
smry.discovery.grade2vs1 = DML(se.discovery_for_grade_2, ~grading2016)
lm_results_discovery_grade2vs1 = summaryExtractTest(smry.discovery.grade2vs1)
str(lm_results_discovery_grade2vs1)

# adjust p values for grade2 vs grade1
padj.grade2 = p.adjust(lm_results_discovery_grade2vs1$Pval_grading20162, method = "BH", n = length(lm_results_discovery_grade2vs1$Pval_grading20162))
lm_results_discovery_grade2vs1 = cbind(lm_results_discovery_grade2vs1, padj.grade2)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade2 over grade1
# 1 hypomethylated,  45 hypermethylated probes in grade2 vs grade1
hypo_grade2 = lm_results_discovery_grade2vs1 %>% dplyr::filter(Est_grading20162 < -0.2 & padj.grade2 < 0.05)
hyper_grade2 = lm_results_discovery_grade2vs1 %>% dplyr::filter(Est_grading20162 > 0.2 & padj.grade2 < 0.05)



#model differences grade3 vs grade1
smry.discovery.grade3vs1 = DML(se.discovery_for_grade_3, ~grading2016)
lm_results_discovery_grade3vs1 = summaryExtractTest(smry.discovery.grade3vs1)
str(lm_results_discovery_grade3vs1)

# adjust p values for grade3 vs grade1
padj.grade3 = p.adjust(lm_results_discovery_grade3vs1$Pval_grading20163, method = "BH", n = length(lm_results_discovery_grade3vs1$Pval_grading20163))
lm_results_discovery_grade3vs1 = cbind(lm_results_discovery_grade3vs1, padj.grade3)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade3 over grade1
# 5517 hypomethylated,  7801 hypermethylated probes in grade3 vs grade1
hypo_grade3 = lm_results_discovery_grade3vs1 %>% dplyr::filter(Est_grading20163 < -0.2 & padj.grade3 < 0.05)
hyper_grade3 = lm_results_discovery_grade3vs1 %>% dplyr::filter(Est_grading20163 > 0.2 & padj.grade3 < 0.05)


#model differences hypermetabolic vs NF2-WT/Immunogenic
smry.discovery.hypermetabolic = DML(se.discovery_for_hypermetabolic, ~MCconsensus)
lm_results_discovery_hypermetabolic = summaryExtractTest(smry.discovery.hypermetabolic)
str(lm_results_discovery_hypermetabolic)

# adjust p values for hypermetabolic vs reference
padj.hypermetabolic = p.adjust(lm_results_discovery_hypermetabolic$Pval_MCconsensushypermetabolic, method = "BH", n = length(lm_results_discovery_hypermetabolic$Pval_MCconsensushypermetabolic))
lm_results_discovery_hypermetabolic = cbind(lm_results_discovery_hypermetabolic, padj.hypermetabolic)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in hypermetabolic vs reference
#  hypomethylated,   hypermethylated probes in grade3 vs reference
hypo_hypermetabolic = lm_results_discovery_hypermetabolic %>% dplyr::filter(Est_MCconsensushypermetabolic < -0.2 & padj.hypermetabolic < 0.05)
hyper_hypermetabolic = lm_results_discovery_hypermetabolic %>% dplyr::filter(Est_MCconsensushypermetabolic > 0.2 & padj.hypermetabolic < 0.05)







#####gene enrichment in differentially methylated probes
#cluster c2 vs c1
genes.hyper.c2 = testEnrichment(hyper_c2_vs_c1$Probe_ID, KYCG_buildGeneDBs(rownames(beta.discovery.flt)))
genes.hyper.c2.sig = genes.hyper.c2 %>% dplyr::filter(FDR < 0.01)
genes.hyper.c2.sig = genes.hyper.c2.sig$gene_name
genes.hyper.c2.sig.df = as.data.frame(genes.hyper.c2.sig)
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_differential_methylation/results")
write.csv(genes.hyper.c2, file = "Cluster2_hypermethylated_genes.csv")
write.csv(genes.hyper.c2.sig, file = "Cluster2_hypermethylated_genes_sig.csv")

#test for manhatten plot
library(qqman)
library(dplyr)
library(plyr)

#make Manhatten plot for hypermethylated genes in cluster c2
data.Man = as.data.frame(annoEPICv2$chr)
data.Man$geneID = annoEPICv2$GencodeV41_Name
data.Man$probe = annoEPICv2$Name
data.Man$gene = annoEPICv2$UCSC_RefGene_Name
data.Man$position = annoEPICv2$pos
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_differential_methylation/results")
write.csv(data.Man, file = "EPIC_anno.csv")

#get data
setwd("/Users/lab/Desktop/Meningioma/T26_discovery_differential_methylation/data")
c2_hyper_Man_data = read.csv(file="c2_hypermethylated_Manhatten_mod.csv", header = T)
colnames(c2_hyper_Man_data)[which(names(c2_hyper_Man_data)=="Bpcum")] = "BPcum"

#make the plot
library(ggplot2)

axisdf = read.csv(file="axisdf.csv")

cairo_pdf(filename = "Manhatten_genes_enriched_hypermethylated_c2.pdf", width = 9, height = 5)
ggplot(c2_hyper_Man_data, aes(x=BPcum, y=P)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.3, size=3) +
  scale_color_manual(values = rep(c("grey23", "grey60"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 1))+
  ylim(0,63)+
  
  # Add highlighted points
  geom_point(data=subset(c2_hyper_Man_data, PCDHG=="yes"), color="red", size=6,alpha=1) +
  geom_point(data=subset(c2_hyper_Man_data, PCDHB=="yes"), color="darkorange", size=6,alpha=1) +
  geom_point(data=subset(c2_hyper_Man_data, PCDHA=="yes"), color="purple3", size=6,alpha=1) +
  geom_point(data=subset(c2_hyper_Man_data, HOXA=="yes"), color="deepskyblue1", size=6) +
  geom_point(data=subset(c2_hyper_Man_data, HOXB=="yes"), color="deepskyblue2", size=6) +
  geom_point(data=subset(c2_hyper_Man_data, HOXC=="yes"), color="deepskyblue3", size=6) +
  geom_point(data=subset(c2_hyper_Man_data, HOXD=="yes"), color="deepskyblue4", size=6) +
  
  # Custom the theme:
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
dev.off()


###############validation of recall of c2 hyper probes in UCSF dataset for grade, group, and c2 cluster metric
#path to idats
idat_dir = "/Users/lab/Desktop/Meningioma/Raleigh_cohort/data"

#get sample data
setwd("/Users/lab/Desktop/Meningioma/Raleigh_cohort/data")
sampledata_UCSF = read.csv(file = "targets_UCSF_predicted.csv", header = T, row.names = 1)
sampledata_UCSF = sampledata_UCSF[,-1]

#employ sesame pipeline to get beta values
beta.UCSF = openSesame(idat_dir)

#re-order beta dataframe columns (Basename) to match Basename order in sampledata
beta.UCSF = beta.UCSF[,match(sampledata_UCSF$Basename,colnames(beta.UCSF))]
all(sampledata_UCSF$Basename %in% colnames(beta.UCSF))
all(sampledata_UCSF$Basename == colnames(beta.UCSF))
beta.UCSF = as.data.frame(beta.UCSF)
colnames(beta.UCSF) = sampledata_UCSF$ID


#remove all probes which are masked in at least one case
beta.UCSF.flt = beta.UCSF[complete.cases(beta.UCSF),]
save.image()

# differential methylation in discovery cohort
#first, make summarizedExperiment
library(SummarizedExperiment)
se.UCSF <- SummarizedExperiment(as.matrix(beta.UCSF.flt), colData = sampledata_UCSF)
str(se.UCSF)

#make contrasts
str(se.UCSF$cluster)
colData(se.UCSF)$cluster <- relevel(factor(colData(se.UCSF)$cluster), "1")
str(se.UCSF)


#subset summarizedexperiment for individual contrasts
se.UCSF_for_grade_2 = se.UCSF[,!se.UCSF$grade == "3"]
str(se.UCSF_for_grade_2)
colData(se.UCSF_for_grade_2)$grade <- relevel(factor(colData(se.UCSF_for_grade_2)$grade), "1")
str(se.UCSF_for_grade_2)

se.UCSF_for_grade_3 = se.UCSF[,!se.UCSF$grade == "2"]
str(se.UCSF_for_grade_3)
colData(se.UCSF_for_grade_3)$grade <- relevel(factor(colData(se.UCSF_for_grade_3)$grade), "1")
str(se.UCSF_for_grade_3)

se.UCSF_for_hypermetabolic = se.UCSF[,!se.UCSF$MCconsensus == "proliferative"]
str(se.UCSF_for_hypermetabolic)
colData(se.UCSF_for_hypermetabolic)$MCconsensus <- relevel(factor(colData(se.UCSF_for_hypermetabolic)$MCconsensus), "Merlinintact")
colData(se.UCSF_for_hypermetabolic)$MCconsensus = combineLevels(se.UCSF_for_hypermetabolic$MCconsensus,
                                                                levs = c("Merlinintact", "Immuneenriched"), newLabel = "other")
colData(se.UCSF_for_hypermetabolic)$MCconsensus <- relevel(factor(colData(se.UCSF_for_hypermetabolic)$MCconsensus), "other")
str(se.UCSF_for_hypermetabolic)

se.UCSF_for_proliferative = se.UCSF[,!se.UCSF$MCconsensus == "hypermetabolic"]
colData(se.UCSF_for_proliferative)$MCconsensus <- relevel(factor(colData(se.UCSF_for_proliferative)$MCconsensus), "Merlinintact")
colData(se.UCSF_for_proliferative)$MCconsensus = combineLevels(se.UCSF_for_proliferative$MCconsensus,
                                                               levs = c("Merlinintact", "Immuneenriched"), newLabel = "other")
colData(se.UCSF_for_proliferative)$MCconsensus <- relevel(factor(colData(se.UCSF_for_proliferative)$MCconsensus), "other")
str(se.UCSF_for_proliferative)




########model differences grade2 vs grade1
smry.UCSF.grade2vs1 = DML(se.UCSF_for_grade_2, ~grade)
lm_results_UCSF_grade2vs1 = summaryExtractTest(smry.UCSF.grade2vs1)
str(lm_results_UCSF_grade2vs1)

# adjust p values for grade2 vs grade1
padj.grade2_UCSF = p.adjust(lm_results_UCSF_grade2vs1$Pval_grade2, method = "BH", n = length(lm_results_UCSF_grade2vs1$Pval_grade2))
lm_results_UCSF_grade2vs1 = cbind(lm_results_UCSF_grade2vs1, padj.grade2_UCSF)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade2 over grade1
# 1 hypomethylated,  18 hypermethylated probes in grade2 vs grade1
hypo_grade2_UCSF = lm_results_UCSF_grade2vs1 %>% dplyr::filter(Est_grade2 < -0.2 & padj.grade2_UCSF < 0.05)
hyper_grade2_UCSF = lm_results_UCSF_grade2vs1 %>% dplyr::filter(Est_grade2 > 0.2 & padj.grade2_UCSF < 0.05)
write.csv(hyper_grade2_UCSF, file = "Grade2_hypermethylated_probes_sig_UCSF.csv")
save.image()





########model differences grade3 vs grade1
smry.UCSF.grade3vs1 = DML(se.UCSF_for_grade_3, ~grade)
lm_results_UCSF_grade3vs1 = summaryExtractTest(smry.UCSF.grade3vs1)
str(lm_results_UCSF_grade3vs1)

# adjust p values for grade2 vs grade1
padj.grade3_UCSF = p.adjust(lm_results_UCSF_grade3vs1$Pval_grade3, method = "BH", n = length(lm_results_UCSF_grade3vs1$Pval_grade3))
lm_results_UCSF_grade3vs1 = cbind(lm_results_UCSF_grade3vs1, padj.grade3_UCSF)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade2 over grade1
# 1 hypomethylated,  18 hypermethylated probes in grade2 vs grade1
hypo_grade3_UCSF = lm_results_UCSF_grade3vs1 %>% dplyr::filter(Est_grade3 < -0.2 & padj.grade3_UCSF < 0.05)
hyper_grade3_UCSF = lm_results_UCSF_grade3vs1 %>% dplyr::filter(Est_grade3 > 0.2 & padj.grade3_UCSF < 0.05)
write.csv(hyper_grade3_UCSF, file = "Grade3_hypermethylated_probes_sig_UCSF.csv")




#model differences hypermetabolic vs NF2-WT/Immunogenic in UCSF
smry.UCSF.hypermetabolic = DML(se.UCSF_for_hypermetabolic, ~MCconsensus)
lm_results_UCSF_hypermetabolic = summaryExtractTest(smry.UCSF.hypermetabolic)
str(lm_results_UCSF_hypermetabolic)

# adjust p values for hypermetabolic vs reference
padj.hypermetabolic_UCSF = p.adjust(lm_results_UCSF_hypermetabolic$Pval_MCconsensushypermetabolic, method = "BH", n = length(lm_results_UCSF_hypermetabolic$Pval_MCconsensushypermetabolic))
lm_results_UCSF_hypermetabolic = cbind(lm_results_UCSF_hypermetabolic, padj.hypermetabolic_UCSF)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in hypermetabolic vs reference
#  12 hypomethylated,   523 hypermethylated probes in grade3 vs reference
hypo_hypermetabolic_UCSF = lm_results_UCSF_hypermetabolic %>% dplyr::filter(Est_MCconsensushypermetabolic < -0.2 & padj.hypermetabolic_UCSF < 0.05)
hyper_hypermetabolic_UCSF = lm_results_UCSF_hypermetabolic %>% dplyr::filter(Est_MCconsensushypermetabolic > 0.2 & padj.hypermetabolic_UCSF < 0.05)
write.csv(hyper_hypermetabolic_UCSF, file = "Hypermetabolic_hypermethylated_probes_sig_UCSF.csv")
save.image()



#model differences hypermetabolic vs NF2-WT/Immunogenic in UCSF
smry.UCSF.proliferative = DML(se.UCSF_for_proliferative, ~MCconsensus)
lm_results_UCSF_proliferative = summaryExtractTest(smry.UCSF.proliferative)
str(lm_results_UCSF_proliferative)

# adjust p values for hypermetabolic vs reference
padj.proliferative_UCSF = p.adjust(lm_results_UCSF_proliferative$Pval_MCconsensusproliferative, method = "BH", n = length(lm_results_UCSF_proliferative$Pval_MCconsensusproliferative))
lm_results_UCSF_proliferative = cbind(lm_results_UCSF_proliferative, padj.proliferative_UCSF)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in hypermetabolic vs reference
#  1877 hypomethylated, 7339 hypermethylated probes in grade3 vs reference
hypo_proliferative_UCSF = lm_results_UCSF_proliferative %>% dplyr::filter(Est_MCconsensusproliferative < -0.2 & padj.proliferative_UCSF < 0.05)
hyper_proliferative_UCSF = lm_results_UCSF_proliferative %>% dplyr::filter(Est_MCconsensusproliferative > 0.2 & padj.proliferative_UCSF < 0.05)
write.csv(hyper_proliferative_UCSF, file = "Proliferative_hypermethylated_probes_sig_UCSF.csv")
save.image()


#model differences c2 vs c1 in UCSF
smry.UCSF.cluster = DML(se.UCSF, ~cluster)
lm_results_UCSF_cluster = summaryExtractTest(smry.UCSF.cluster)
str(lm_results_UCSF_cluster)

# adjust p values for c2 vs c1
padj.cluster_UCSF = p.adjust(lm_results_UCSF_cluster$Pval_cluster2, method = "BH", n = length(lm_results_UCSF_cluster$Pval_cluster2))
lm_results_UCSF_cluster = cbind(lm_results_UCSF_cluster, padj.cluster_UCSF)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in c2 vs c1
#  681 hypomethylated, 8757 hypermethylated probes in grade3 vs reference
hypo_cluster_UCSF = lm_results_UCSF_cluster %>% dplyr::filter(Est_cluster2 < -0.2 & padj.cluster_UCSF < 0.05)
hyper_cluster_UCSF = lm_results_UCSF_cluster %>% dplyr::filter(Est_cluster2 > 0.2 & padj.cluster_UCSF < 0.05)
write.csv(hyper_cluster_UCSF, file = "Cluster2_hypermethylated_probes_sig_UCSF.csv")
save.image()


####visualize recall of c2 hypermethylated probes in contrasts based on grade and group
#get data
Recall_discovery = read.csv(file="Recall_c2_probes_discovery.csv", header = T)
str(Recall_discovery)
Recall_discovery$condition = factor(Recall_discovery$condition, levels = c("Grade2", "Grade3",
                                                                           "Hypermetabolic",
                                                                           "Proliferative"))
str(Recall_discovery)
cairo_pdf(filename = "Recall_discovery_c2_methylated_3.pdf", width = 5.2, height = 7.5)
ggplot(Recall_discovery, 
       aes(condition, recall, fill=condition))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#8C6BB1","#810F7C","forestgreen","darkorange2")) +
  scale_y_continuous(trans="sqrt")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

Recall_UCSF = read.csv(file="Recall_c2_probes_UCSF.csv", header = T)
str(Recall_UCSF)
Recall_UCSF$condition = factor(Recall_UCSF$condition, levels = c("Grade2", "Grade3",
                                                                 "Hypermetabolic",
                                                                 "Proliferative","c2"))
cairo_pdf(filename = "Recall_UCSF_c2_methylated.pdf", width = 5.1, height = 7.5)
ggplot(Recall_UCSF, 
       aes(condition, recall, fill=condition))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#8C6BB1","#810F7C","forestgreen","darkorange2","violetred4")) +
  scale_y_continuous(trans="sqrt")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()



#make track view for genomic regions or genes
#cluster
cairo_pdf(filename = "Track_view_PCDHG@.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',141330685,141512979, beta.discovery.flt_ordered,
  show.probeNames = FALSE)
dev.off()

cairo_pdf(filename = "Track_view_PCDHB@.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',141051394,141248234, beta.discovery.flt_ordered,
  show.probeNames = FALSE)
dev.off()

cairo_pdf(filename = "Track_view_PCDHA@.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',140786136,141012344, beta.discovery.flt_ordered,
  show.probeNames = FALSE)
dev.off()

cairo_pdf(filename = "Track_view_PCDH_clusters.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',140786136,141512979, beta.discovery.flt_ordered,
  show.probeNames = FALSE, show.sampleNames = FALSE, heat.height=0.3)
dev.off()

save.image()

cairo_pdf(filename = "Track_view_HOXA.pdf", width = 8, height = 6)
visualizeRegion(
  'chr7',27041605,27172274, beta.discovery.flt_ordered,
  show.probeNames = FALSE)
dev.off()


#gene-specific track view

visualizeGene('PITX1', beta.discovery.flt_ordered,show.probeNames = FALSE)

cairo_pdf(filename = "Track_view_PITX1.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',135023734,135039000, beta.discovery.flt_ordered,
  show.probeNames = FALSE)
dev.off()


visualizeGene('PCDHGC3', beta.discovery.flt_ordered,show.probeNames = FALSE)

cairo_pdf(filename = "Track_view_PCDHGC3-5.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',141465000,141515100, beta.discovery.flt_ordered,
  show.probeNames = FALSE)
dev.off()


visualizeGene('SIM2', beta.discovery.flt_ordered,show.probeNames = FALSE)





