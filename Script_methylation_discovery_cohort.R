library(minfi)
library(RColorBrewer)
library(limma)
library(sesame)
library(gdata)
library(ggplot2)
library(devtools)
library(ggstatsplot)
library(gginnards)
library(forcats)
library(ggpubr)
library(dplyr)
library(ggExtra)

#read in data (n=231 MNG)
idat_dir = "Directory_discovery"
targets <- read.metharray.sheet(idat_dir, pattern="targets_discovery.csv")
rgSet <- read.metharray.exp(base = idat_dir, targets=targets)
rgSet
sampleNames(rgSet) <- targets$ID
rgSet
annotation(rgSet) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
rgSet

###Normalization
mSetSq <- preprocessNoob(rgSet) 
#no normalization control
mSetRaw <- preprocessRaw(rgSet)
gmSetRaw = mapToGenome(mSetRaw)

#filtering
pvalspOOBAH = openSesame(idat_dir, func = pOOBAH, return.pval = TRUE)
detP = as.data.frame(pvalspOOBAH)
detP <- detP[,match(targets$Basename, colnames(detP))]
colnames(detP) = targets$ID
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
keep <- rowSums(detP < 0.05) == ncol(mSetSq) 
table(keep)
#based on detP, 375,314 probes are excluded, 561,676 kept

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

#make MethylSet a GenomicRatioSet
mSetSqFlt = mapToGenome(mSetSqFlt)
mSetSqFlt = ratioConvert(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
#probes after SNP drop (n=548,339)

# tag sex chromosome probes for removal
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
data("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

keep <- !(featureNames(mSetSqFlt) %in% annoEPICv2$Name[annoEPICv2$chr %in% 
                                                         c("chrX","chrY")])
table(keep)
#10,235 probes in sex chromosomes

#visualize effect of gender probes on MDS
pal6 = brewer.pal(11, "PiYG")[c(2,11)]

setwd("/Users/lab/Desktop/Meningioma/T26_Discovery_QC_cluster/results")
tiff("MDS_effect_gender_probes.tiff", res=300,width=5120,height=2048)
par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=10000, pch = 19,cex = 3,gene.selection="common", 
        col=pal6[factor(targets$gender)], main="With Sex CHR Probes")
legend("topright",  inset = c(-0.08, 0),legend=levels(factor(targets$gender)), text.col=pal6,
       cex=1, box.lty=0,xpd = TRUE)
plotMDS(getM(mSetSqFlt[keep,]), top=10000, pch = 19,cex = 3,gene.selection="common", 
        col=pal6[factor(targets$gender)],main="Without Sex CHR Probes")
legend("topright",  inset = c(-0.08, 0),legend=levels(factor(targets$gender)), text.col=pal6,
       cex=1, box.lty=0,xpd = TRUE)
dev.off()

#sex plot suggest major variation
#clear separation of samples based on gender-specific probes, thus exclude gender probes
#final probe number n=538,104
mSetSqFlt
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt

#get B values
bVals <- getBeta(mSetSqFlt)
bVals = as.data.frame(bVals)

#determine minimal bVals subset to capture main components
#make subset of most variable beta values
bVals.sub = bVals
bVals.sub$var = apply(bVals.sub,1,var)
bVals.sub <- bVals.sub[order(bVals.sub$var, decreasing = TRUE),]
bVals.sub = bVals.sub[,-232]
bVals.sub.10k = bVals.sub[1:10000,]
bVals.sub.50k = bVals.sub[1:50000,]
bVals.sub.100k = bVals.sub[1:100000,]
bVals.sub.200k = bVals.sub[1:200000,]

#make for all betas or top variable subsets
bVals.pca = prcomp(t(bVals), center = T, scale. = F)
bVals.sub.10k.pca = prcomp(t(bVals.sub.10k), center = T, scale. = F)
bVals.sub.50k.pca = prcomp(t(bVals.sub.50k), center = T, scale. = F)
bVals.sub.100k.pca = prcomp(t(bVals.sub.100k), center = T, scale. = F)
bVals.sub.200k.pca = prcomp(t(bVals.sub.200k), center = T, scale. = F)

setwd("/Users/lab/Desktop/Meningioma/T26_Discovery_QC_cluster/results")
tiff("Scree_PCs_top10k.tiff", res=300,width=1560,height=2048)
fviz_screeplot(bVals.sub.10k.pca)
dev.off()

tiff("Scree_PCs_top50k.tiff", res=300,width=1560,height=2048)
fviz_screeplot(bVals.sub.50k.pca)
dev.off()

tiff("Scree_PCs_top100k.tiff", res=300,width=1560,height=2048)
fviz_screeplot(bVals.sub.100k.pca)
dev.off()

tiff("Scree_PCs_top200k.tiff", res=300,width=1560,height=2048)
fviz_screeplot(bVals.sub.200k.pca)
dev.off()

#make pca plots for various conditions for 50k probes
df <- cbind(targets, bVals.sub.50k.pca$x[,1:3])
df$grading2016 = as.character(df$grading2016)
df$grading2021 = as.character(df$grading2021)
df$cluster_2000 = factor(df$cluster_2000, levels = c("1","2"))
df[df=="SMARCE1-altered"]<-"SMARCE1altered"
df[df=="Immune-enriched"]<-"Immuneenriched"
df[df=="Merlin-intact"]<-"Merlinintact"
df$grading2021 = factor(df$grading2021, levels = c("1", "2", "3"))


tiff("PCA_50k_margin_MCconsensus_PC1-2.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC2",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("forestgreen","red3","royalblue2","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("test.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC2",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("forestgreen","red3","royalblue2","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_MCconsensus_PC1-3.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC3",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("forestgreen","red3","royalblue2","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_grade_PC1-2.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC2",
  color = "grading2021", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_grade_PC1-3.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC3",
  color = "grading2021", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_risk_PC1-2.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC2",
  color = "risk_score", size = 5, alpha = 0.8,
  palette = c("firebrick1","purple1","dodgerblue2"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_risk_PC1-3.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC3",
  color = "risk_score", size = 5, alpha = 0.8,
  palette = c("firebrick1","purple1","dodgerblue2"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_status_PC1-2.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC2",
  color = "status", size = 5, alpha = 0.8,
  palette = c("#1B9E77","#D95F02"),
  margin.params = list(fill = "status", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()
tiff("PCA_50k_margin_status_PC1-3.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC3",
  color = "status", size = 5, alpha = 0.8,
  palette = c("#1B9E77","#D95F02"),
  margin.params = list(fill = "status", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_MCsubtype_PC1-2.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC2",
  color = "MCsubtype", size = 5, alpha = 0.8,
  palette = c("#386CB0","#F0027F","#BF5B17","#666666"),
  margin.params = list(fill = "MCsubtype", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()
tiff("PCA_50k_margin_MCsubtype_PC1-3.tiff", res=300,width=2460,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC3",
  color = "MCsubtype", size = 5, alpha = 0.8,
  palette = c("#386CB0","#F0027F","#BF5B17","#666666"),
  margin.params = list(fill = "MCsubtype", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_cluster_PC1-2.tiff", res=300,width=3560,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC2",
  color = "cluster_2000", size = 5, alpha = 0.8,
  palette = c("cyan4","violetred4"),
  margin.params = list(fill = "cluster_2000", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_cluster_PC1-3.tiff", res=300,width=3560,height=2048)
ggscatterhist(
  df, x = "PC1", y = "PC3",
  color = "cluster_2000", size = 5, alpha = 0.8,
  palette = c("cyan4","violetred4"),
  margin.params = list(fill = "cluster_2000", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

###########get probes that contribute most to PC1, as PC1 best captures difference between high grade and low grade
var = get_pca_var(bVals.pca)
head(var$contrib,10)

var.contrib = as.data.frame(var$contrib)
var.contrib.sortPC1 = var.contrib[order(var.contrib$Dim.1, decreasing = TRUE),]

top_PC1 = data.frame(PC1=rownames(var.contrib.sortPC1))

top1000.features = c(top_PC1$PC1[1:1000])
top2000.features = c(top_PC1$PC1[1:2000])
top3000.features = c(top_PC1$PC1[1:3000])

write.csv(top2000.features, file="top2000_features.csv")

#subset beta values by top features for PC1
beta.PC1.top1000 = bVals %>% dplyr::filter(rownames(bVals) %in% top1000.features)
beta.PC1.top2000 = bVals %>% dplyr::filter(rownames(bVals) %in% top2000.features)
beta.PC1.top3000 = bVals %>% dplyr::filter(rownames(bVals) %in% top3000.features)


#unsupervised hierarchical clustering

T26_samples = as.data.frame(targets$ID)
colnames(T26_samples) = "ID"
rownames(T26_samples) = T26_samples$ID
all(rownames(T26_samples) %in% targets$ID)
all(rownames(T26_samples) == targets$ID)
T26_samples$risk = targets$risk_score
T26_samples$MCconsensus = targets$MCconsensus
T26_samples$MCsubtype = targets$MCsubtype
T26_samples$grade = targets$grading2021
T26_samples$grade = as.character(T26_samples$grade)
T26_samples$status = targets$status
T26_samples$gender = targets$gender
T26_samples[T26_samples=="SMARCE1-altered"]<-"SMARCE1altered"
T26_samples[T26_samples=="Immune-enriched"]<-"Immuneenriched"
T26_samples[T26_samples=="Merlin-intact"]<-"Merlinintact"
T26_samples = T26_samples[,-1]


#determine optimal number fo clusters
library(factoextra)
set.seed(123)
tiff("Elbow_clusters_top1000.tiff", res=300,width=2560,height=1548)
fviz_nbclust(beta.PC1.top1000, kmeans,  method = "wss")+
  geom_vline(xintercept = 2, linetype=2)
dev.off()

set.seed(123)
tiff("Elbow_clusters_top2000.tiff", res=300,width=2560,height=1548)
fviz_nbclust(beta.PC1.top2000, kmeans,  method = "wss")+
  geom_vline(xintercept = 2, linetype=2)
dev.off()

set.seed(123)
tiff("Elbow_clusters_top3000.tiff", res=300,width=2560,height=1548)
fviz_nbclust(beta.PC1.top3000, kmeans,  method = "wss")+
  geom_vline(xintercept = 2, linetype=2)
dev.off()

#heatmaps

ann_colors3 = list(gender = c(M="#276419", F="#C51B7D"),
                   grade=c("1" = "#9EBCDA", "2" = "#8C6BB1", "3" = "#810F7C" ),
                   MCsubtype =c(Benign = "#386CB0", Intermediate = "#F0027F", 
                                Malignant = "#BF5B17",SMARCE1altered = "#666666"),
                   status = c(Primary="#1B9E77", Recurrence="#D95F02"),
                   MCgroup = c(Hypermitotic = "red",
                               Immuneenriched = "darkmagenta",
                               Merlinintact = "blue3"),
                   MCconsensus = c(Immuneenriched = "red3",
                                   Merlinintact = "royalblue2",
                                   hypermetabolic = "forestgreen",
                                   proliferative = "darkorange2"))

ann_colors4 = list(gender = c(M="#276419", F="#C51B7D"),
                   grade=c("1" = "#9EBCDA", "2" = "#8C6BB1", "3" = "#810F7C" ),
                   MCsubtype =c(Benign = "#386CB0", Intermediate = "#F0027F", 
                                Malignant = "#BF5B17",SMARCE1altered = "#666666"),
                   risk =c(low="dodgerblue2", intermediate = "purple1",
                           high = "firebrick1"),
                   status = c(Primary="#1B9E77", Recurrence="#D95F02"),
                   MCconsensus = c(Immuneenriched = "red3",
                                   Merlinintact = "royalblue2",
                                   hypermetabolic = "forestgreen",
                                   proliferative = "darkorange2"))



#final plots
library(pheatmap)

tiff("Heat_cluster_PC1_top1000.tiff", res=300,width=2560,height=2048)
pheatmap(beta.PC1.top1000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors4,
         cutree_cols = 2, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top2000.tiff", res=300,width=2560,height=2048)
pheatmap(beta.PC1.top2000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors4,
         cutree_cols = 2, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top3000.tiff", res=300,width=2560,height=2048)
pheatmap(beta.PC1.top3000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors4,
         cutree_cols = 2, fontsize = 8)
dev.off()

#chisquare for primary/recurrence
library(ggstatsplot)
#get clustering results

cluster_1000 = pheatmap(beta.PC1.top1000)
x = cutree(clust_1000$tree_col, k=2)
targets = cbind(targets, cluster_1000 = x)

cluster_2000 = pheatmap(beta.PC1.top2000)
x = cutree(cluster_2000$tree_col, k=2)
targets = cbind(targets, cluster_2000 = x)

cluster_3000 = pheatmap(beta.PC1.top3000)
x = cutree(cluster_3000$tree_col, k = 2)
targets = cbind(targets, cluster_3000 = x)

write.csv(targets, file="PC1_cluster_results_new_grade.csv")

#load file with re-named clusters
targets = read.csv(file = "PC1_cluster_results_new_grade.csv", header = T)
clusters = read.csv(file = "PC1_cluster_results_new_grade.csv", header = T)

#clust_1000
length(which(clusters$status == "Recurrence" & clusters$cluster_1000 == 1))
length(which(clusters$status == "Recurrence" & clusters$cluster_1000 == 2))
length(which(clusters$status == "Primary" & clusters$cluster_1000 == 1))
length(which(clusters$status == "Primary" & clusters$cluster_1000 == 2))

dat1000 <- data.frame(
  c(139, 19),
  c(29, 44),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat1000) <- c("Primary", "Recurrence")
dat1000

x <- c()
for (row in rownames(dat1000)) {
  for (col in colnames(dat1000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat1000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df1000 <- as.data.frame(x)
colnames(df1000) <- c("clustering", "status")
df1000

set.seed(123)
test1000 <- chisq.test(table(df1000))

tiff("Chisquare_clust1000_status_new.tiff", res=300,width=1000,height=2048)
ggbarstats(df1000, status, clustering,
           legend.position = "none",
           results.subtitle = FALSE,
           subtitle = paste0(test1000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#D95F02","#1B9E77"))+
  theme(legend.position = "none")
dev.off()


#clust_2000
length(which(clusters$status == "Recurrence" & clusters$cluster_2000 == 1))
length(which(clusters$status == "Recurrence" & clusters$cluster_2000 == 2))
length(which(clusters$status == "Primary" & clusters$cluster_2000 == 1))
length(which(clusters$status == "Primary" & clusters$cluster_2000 == 2))

dat2000 <- data.frame(
  c(141, 17),
  c(29, 44),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat2000) <- c("Primary", "Recurrence")
dat2000

x <- c()
for (row in rownames(dat2000)) {
  for (col in colnames(dat2000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat2000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df2000 <- as.data.frame(x)
colnames(df2000) <- c("clustering", "status")
df2000

set.seed(123)
test2000 <- chisq.test(table(df2000))

tiff("Chisquare_clust2000_status.tiff", res=300,width=1000,height=2048)
ggbarstats(df2000, status, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test2000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#D95F02","#1B9E77"))+
  theme(legend.position = "none")
dev.off()

#clust_3000
length(which(clusters$status == "Recurrence" & clusters$cluster_3000 == 1))
length(which(clusters$status == "Recurrence" & clusters$cluster_3000 == 2))
length(which(clusters$status == "Primary" & clusters$cluster_3000 == 1))
length(which(clusters$status == "Primary" & clusters$cluster_3000 == 2))

dat3000 <- data.frame(
  c(139, 19),
  c(29, 44),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat3000) <- c("Primary", "Recurrence")
dat3000

x <- c()
for (row in rownames(dat3000)) {
  for (col in colnames(dat3000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat3000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df3000 <- as.data.frame(x)
colnames(df3000) <- c("clustering", "status")
df3000

set.seed(123)
test3000 <- chisq.test(table(df3000))

tiff("Chisquare_clust3000_status.tiff", res=300,width=1000,height=2048)
ggbarstats(df3000, status, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test3000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#D95F02","#1B9E77"))+
  theme(legend.position = "none")
dev.off()

#add chisquare for grade
#clust1000
length(which(clusters$grading2021 == "1" & clusters$cluster_1000 == 2))
length(which(clusters$grading2021 == "2" & clusters$cluster_1000 == 2))
length(which(clusters$grading2021 == "3" & clusters$cluster_1000 == 2))
length(which(clusters$grading2021 == "1" & clusters$cluster_1000 == 1))
length(which(clusters$grading2021 == "2" & clusters$cluster_1000 == 1))
length(which(clusters$grading2021 == "3" & clusters$cluster_1000 == 1))

dat1000 <- data.frame(
  c(76, 7),
  c(89, 40),
  c(3,16),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat1000) <- c("1", "2", "3")
dat1000

x <- c()
for (row in rownames(dat1000)) {
  for (col in colnames(dat1000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat1000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df1000 <- as.data.frame(x)
colnames(df1000) <- c("clustering", "grade")
df1000

test1000 <- fisher.test(table(df1000))

tiff("Chisquare_clust1000_grade.tiff", res=300,width=1000,height=2048)
ggbarstats(df1000, grade, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test1000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#810F7C", "#8C6BB1", "#9EBCDA"))+
  theme(legend.position = "none")
dev.off()

#clust2000
length(which(clusters$grading2021 == "1" & clusters$cluster_2000 == 2))
length(which(clusters$grading2021 == "2" & clusters$cluster_2000 == 2))
length(which(clusters$grading2021 == "3" & clusters$cluster_2000 == 2))
length(which(clusters$grading2021 == "1" & clusters$cluster_2000 == 1))
length(which(clusters$grading2021 == "2" & clusters$cluster_2000 == 1))
length(which(clusters$grading2021 == "3" & clusters$cluster_2000 == 1))

dat2000 <- data.frame(
  c(77, 6),
  c(90, 39),
  c(3,16),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat2000) <- c("1", "2", "3")
dat2000

x <- c()
for (row in rownames(dat2000)) {
  for (col in colnames(dat2000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat2000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df2000 <- as.data.frame(x)
colnames(df2000) <- c("clustering", "grade")
df2000

test2000 <- fisher.test(table(df2000))

tiff("Chisquare_clust2000_grade.tiff", res=300,width=1000,height=2048)
ggbarstats(df2000, grade, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test2000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#810F7C", "#8C6BB1", "#9EBCDA"))+
  theme(legend.position = "none")
dev.off()

#clust3000
length(which(clusters$grading2021 == "1" & clusters$cluster_3000 == 2))
length(which(clusters$grading2021 == "2" & clusters$cluster_3000 == 2))
length(which(clusters$grading2021 == "3" & clusters$cluster_3000 == 2))
length(which(clusters$grading2021 == "1" & clusters$cluster_3000 == 1))
length(which(clusters$grading2021 == "2" & clusters$cluster_3000 == 1))
length(which(clusters$grading2021 == "3" & clusters$cluster_3000 == 1))

dat3000 <- data.frame(
  c(76, 7),
  c(89, 40),
  c(3,16),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat3000) <- c("1", "2", "3")
dat3000

x <- c()
for (row in rownames(dat3000)) {
  for (col in colnames(dat3000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat3000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df3000 <- as.data.frame(x)
colnames(df3000) <- c("clustering", "grade")
df3000

test3000 <- fisher.test(table(df3000))

tiff("Chisquare_clust3000_grade.tiff", res=300,width=1000,height=2048)
ggbarstats(df3000, grade, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test3000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#810F7C", "#8C6BB1", "#9EBCDA"))+
  theme(legend.position = "none")
dev.off()

#add chisquare for Heidelberg groups 
#clust1000
length(which(clusters$MCsubtype == "Benign" & clusters$cluster_1000 == 2))
length(which(clusters$MCsubtype == "Intermediate" & clusters$cluster_1000 == 2))
length(which(clusters$MCsubtype == "Malignant" & clusters$cluster_1000 == 2))
length(which(clusters$MCsubtype == "Benign" & clusters$cluster_1000 == 1))
length(which(clusters$MCsubtype == "Intermediate" & clusters$cluster_1000 == 1))
length(which(clusters$MCsubtype == "Malignant" & clusters$cluster_1000 == 1))

dat1000 <- data.frame(
  c(137, 0),
  c(31, 55),
  c(0,7),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat1000) <- c("Benign", "Intermediate", "Malignant")
dat1000

x <- c()
for (row in rownames(dat1000)) {
  for (col in colnames(dat1000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat1000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df1000 <- as.data.frame(x)
colnames(df1000) <- c("clustering", "MCsubtype")
df1000

test1000 <- fisher.test(table(df1000))

tiff("Chisquare_clust1000_MCsubtype.tiff", res=300,width=1000,height=2048)
ggbarstats(df1000, MCsubtype, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test1000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#BF5B17", "#F0027F", "#386CB0"))+
  theme(legend.position = "none")
dev.off()

#clust2000
length(which(clusters$MCsubtype == "Benign" & clusters$cluster_2000 == 2))
length(which(clusters$MCsubtype == "Intermediate" & clusters$cluster_2000 == 2))
length(which(clusters$MCsubtype == "Malignant" & clusters$cluster_2000 == 2))
length(which(clusters$MCsubtype == "Benign" & clusters$cluster_2000 == 1))
length(which(clusters$MCsubtype == "Intermediate" & clusters$cluster_2000 == 1))
length(which(clusters$MCsubtype == "Malignant" & clusters$cluster_2000 == 1))

dat2000 <- data.frame(
  c(137, 0),
  c(33, 53),
  c(0,7),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat2000) <- c("Benign", "Intermediate", "Malignant")
dat2000

x <- c()
for (row in rownames(dat2000)) {
  for (col in colnames(dat2000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat2000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df2000 <- as.data.frame(x)
colnames(df2000) <- c("clustering", "MCsubtype")
df2000

test2000 <- fisher.test(table(df2000))

tiff("Chisquare_clust2000_MCsubtype.tiff", res=300,width=1000,height=2048)
ggbarstats(df2000, MCsubtype, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test2000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#BF5B17", "#F0027F", "#386CB0"))+
  theme(legend.position = "none")
dev.off()

#clust3000
length(which(clusters$MCsubtype == "Benign" & clusters$cluster_3000 == 2))
length(which(clusters$MCsubtype == "Intermediate" & clusters$cluster_3000 == 2))
length(which(clusters$MCsubtype == "Malignant" & clusters$cluster_3000 == 2))
length(which(clusters$MCsubtype == "Benign" & clusters$cluster_3000 == 1))
length(which(clusters$MCsubtype == "Intermediate" & clusters$cluster_3000 == 1))
length(which(clusters$MCsubtype == "Malignant" & clusters$cluster_3000 == 1))

dat3000 <- data.frame(
  c(137, 0),
  c(31, 55),
  c(0,7),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat3000) <- c("Benign", "Intermediate", "Malignant")
dat3000

x <- c()
for (row in rownames(dat3000)) {
  for (col in colnames(dat3000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat3000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df3000 <- as.data.frame(x)
colnames(df3000) <- c("clustering", "MCsubtype")
df3000

test3000 <- fisher.test(table(df3000))

tiff("Chisquare_clust3000_MCsubtype.tiff", res=300,width=1000,height=2048)
ggbarstats(df3000, MCsubtype, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test3000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#BF5B17", "#F0027F", "#386CB0"))+
  theme(legend.position = "none")
dev.off()

#add chisquare for Raleigh groups 
#clust1000

length(which(clusters$MCconsensus == "Merlinintact" & clusters$cluster_1000 == 2))
length(which(clusters$MCconsensus == "Immuneenriched" & clusters$cluster_1000 == 2))
length(which(clusters$MCconsensus == "hypermetabolic" & clusters$cluster_1000 == 2))
length(which(clusters$MCconsensus == "proliferative" & clusters$cluster_1000 == 2))
length(which(clusters$MCconsensus == "Merlinintact" & clusters$cluster_1000 == 1))
length(which(clusters$MCconsensus == "Immuneenriched" & clusters$cluster_1000 == 1))
length(which(clusters$MCconsensus == "hypermetabolic" & clusters$cluster_1000 == 1))
length(which(clusters$MCconsensus == "proliferative" & clusters$cluster_1000 == 1))

dat1000 <- data.frame(
  c(85, 1),
  c(51, 2),
  c(29,17),
  c(3, 43),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat1000) <- c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative")
dat1000

x <- c()
for (row in rownames(dat1000)) {
  for (col in colnames(dat1000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat1000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df1000 <- as.data.frame(x)
colnames(df1000) <- c("clustering", "MCconsensus")
df1000
df1000$MCconsensus = factor(df1000$MCconsensus, levels =  c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

test1000 <- chisq.test(table(df1000))

tiff("Chisquare_clust1000_MCconsensus.tiff", res=300,width=1000,height=2048)
ggbarstats(df1000, MCconsensus, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test1000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("darkorange2","forestgreen","red3","royalblue2"))+
  theme(legend.position = "none")
dev.off()

#clust2000
length(which(clusters$MCconsensus == "Merlinintact" & clusters$cluster_2000 == 2))
length(which(clusters$MCconsensus == "Immuneenriched" & clusters$cluster_2000 == 2))
length(which(clusters$MCconsensus == "hypermetabolic" & clusters$cluster_2000 == 2))
length(which(clusters$MCconsensus == "proliferative" & clusters$cluster_2000 == 2))
length(which(clusters$MCconsensus == "Merlinintact" & clusters$cluster_2000 == 1))
length(which(clusters$MCconsensus == "Immuneenriched" & clusters$cluster_2000 == 1))
length(which(clusters$MCconsensus == "hypermetabolic" & clusters$cluster_2000 == 1))
length(which(clusters$MCconsensus == "proliferative" & clusters$cluster_2000 == 1))

dat2000 <- data.frame(
  c(85, 1),
  c(52, 1),
  c(30,16),
  c(3, 43),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat2000) <- c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative")
dat2000

x <- c()
for (row in rownames(dat2000)) {
  for (col in colnames(dat2000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat2000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df2000 <- as.data.frame(x)
colnames(df2000) <- c("clustering", "MCconsensus")
df2000
df2000$MCconsensus = factor(df2000$MCconsensus, levels =  c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

test2000 <- chisq.test(table(df2000))

tiff("Chisquare_clust2000_MCconsensus.tiff", res=300,width=1000,height=2048)
ggbarstats(df2000, MCconsensus, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test2000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("darkorange2","forestgreen","red3","royalblue2"))+
  theme(legend.position = "none")
dev.off()

#clust3000
length(which(clusters$MCconsensus == "Merlinintact" & clusters$cluster_3000 == 2))
length(which(clusters$MCconsensus == "Immuneenriched" & clusters$cluster_3000 == 2))
length(which(clusters$MCconsensus == "hypermetabolic" & clusters$cluster_3000 == 2))
length(which(clusters$MCconsensus == "proliferative" & clusters$cluster_3000 == 2))
length(which(clusters$MCconsensus == "Merlinintact" & clusters$cluster_3000 == 1))
length(which(clusters$MCconsensus == "Immuneenriched" & clusters$cluster_3000 == 1))
length(which(clusters$MCconsensus == "hypermetabolic" & clusters$cluster_3000 == 1))
length(which(clusters$MCconsensus == "proliferative" & clusters$cluster_3000 == 1))

dat3000 <- data.frame(
  c(84, 2),
  c(52, 1),
  c(29,17),
  c(3, 43),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat3000) <- c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative")
dat3000

x <- c()
for (row in rownames(dat3000)) {
  for (col in colnames(dat3000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat3000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df3000 <- as.data.frame(x)
colnames(df3000) <- c("clustering", "MCconsensus")
df3000
df3000$MCconsensus = factor(df3000$MCconsensus, levels =  c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

test3000 <- chisq.test(table(df3000))

tiff("Chisquare_clust3000_MCconsensus.tiff", res=300,width=1000,height=2048)
ggbarstats(df3000, MCconsensus, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test3000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("darkorange2","forestgreen","red3","royalblue2"))+
  theme(legend.position = "none")
dev.off()


#add chisquare for risk_score
#clust1000
length(which(clusters$risk_score == "low" & clusters$cluster_1000 == 2))
length(which(clusters$risk_score == "intermediate" & clusters$cluster_1000 == 2))
length(which(clusters$risk_score == "high" & clusters$cluster_1000 == 2))
length(which(clusters$risk_score == "low" & clusters$cluster_1000 == 1))
length(which(clusters$risk_score == "intermediate" & clusters$cluster_1000 == 1))
length(which(clusters$risk_score == "high" & clusters$cluster_1000 == 1))

dat1000 <- data.frame(
  c(118, 2),
  c(49, 38),
  c(1,23),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat1000) <- c("low", "intermediate", "high")
dat1000

x <- c()
for (row in rownames(dat1000)) {
  for (col in colnames(dat1000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat1000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df1000 <- as.data.frame(x)
colnames(df1000) <- c("clustering", "risk")
df1000
df1000$risk = factor(df1000$risk, levels =  c("low", "intermediate", "high"))

test1000 <- fisher.test(table(df1000))

tiff("Chisquare_clust1000_risk.tiff", res=300,width=1000,height=2048)
ggbarstats(df1000, risk, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test1000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("firebrick1","purple1", "dodgerblue2"))+
  theme(legend.position = "none")
dev.off()


#clust2000
length(which(clusters$risk_score == "low" & clusters$cluster_2000 == 2))
length(which(clusters$risk_score == "intermediate" & clusters$cluster_2000 == 2))
length(which(clusters$risk_score == "high" & clusters$cluster_2000 == 2))
length(which(clusters$risk_score == "low" & clusters$cluster_2000 == 1))
length(which(clusters$risk_score == "intermediate" & clusters$cluster_2000 == 1))
length(which(clusters$risk_score == "high" & clusters$cluster_2000 == 1))

dat2000 <- data.frame(
  c(118, 2),
  c(51, 36),
  c(1,23),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat2000) <- c("low", "intermediate", "high")
dat2000

x <- c()
for (row in rownames(dat2000)) {
  for (col in colnames(dat2000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat2000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df2000 <- as.data.frame(x)
colnames(df2000) <- c("clustering", "risk")
df2000
df2000$risk = factor(df2000$risk, levels =  c("low", "intermediate", "high"))

test2000 <- fisher.test(table(df2000))

tiff("Chisquare_clust2000_risk.tiff", res=300,width=1000,height=2048)
ggbarstats(df2000, risk, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test2000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("firebrick1","purple1", "dodgerblue2"))+
  theme(legend.position = "none")
dev.off()

#clust3000
length(which(clusters$risk_score == "low" & clusters$cluster_3000 == 2))
length(which(clusters$risk_score == "intermediate" & clusters$cluster_3000 == 2))
length(which(clusters$risk_score == "high" & clusters$cluster_3000 == 2))
length(which(clusters$risk_score == "low" & clusters$cluster_3000 == 1))
length(which(clusters$risk_score == "intermediate" & clusters$cluster_3000 == 1))
length(which(clusters$risk_score == "high" & clusters$cluster_3000 == 1))

dat3000 <- data.frame(
  c(118, 2),
  c(49, 38),
  c(1,23),
  row.names = c("cluster1", "cluster2"),
  stringsAsFactors = FALSE
)
colnames(dat3000) <- c("low", "intermediate", "high")
dat3000

x <- c()
for (row in rownames(dat3000)) {
  for (col in colnames(dat3000)) {
    x <- rbind(x, matrix(rep(c(row, col), dat3000[row, col]), ncol = 2, byrow = TRUE))
  }
}
df3000 <- as.data.frame(x)
colnames(df3000) <- c("clustering", "risk")
df3000
df3000$risk = factor(df3000$risk, levels =  c("low", "intermediate", "high"))

test3000 <- fisher.test(table(df3000))

tiff("Chisquare_clust3000_risk.tiff", res=300,width=1000,height=2048)
ggbarstats(df3000, risk, clustering,
           results.subtitle = FALSE,
           subtitle = paste0(test3000$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("firebrick1","purple1", "dodgerblue2"))+
  theme(legend.position = "none")
dev.off()













