library(DESeq2)
library(ggplot2)
library(pheatmap)
library(factoextra)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(forcats)
library(tibble)
library(DEGreport)
library(stringr)
library(ggvenn)
library(devtools)
library(canceR)

####get raw read counts data and meta data for all 26 samples
data_26 = read.table("data/gene_counts_26_samples.txt", header=T)
data_26 = data_26[,-1]
data_26 <- aggregate(data_26[-1], by = list(data_26$gene_name), FUN=sum)
colnames(data_26)[which(names(data_26) == "Group.1")] <- "gene_name"
str(data_26)

#round all expression values and name rownames
data_26[,-1] <-round(data_26[,-1],0)
row.names(data_26) <- data_26$gene_name
data_26 = data_26[,-1]

#read in meta data
meta_26 = read.table("meta/meta_26_samples.txt", header=T, row.names=1)
meta_26[] = lapply(meta_26, factor)
str(meta_26)

#check and change sample names
all(colnames(data_26) %in% meta_26$QBIC_ID)
all(colnames(data_26) == meta_26$QBIC_ID)
colnames(data_26) = rownames(meta_26)
all(colnames(data_26) == rownames(meta_26))


##start analyzing the entire meta_26##start analyzing the entire dataset
## Create DESeq2Dataset object
dds_26 <- DESeqDataSetFromMatrix(countData = data_26, colData = meta_26, 
                                  design = ~cluster)
View(counts(dds_26))
dim(dds_26)

###Normalization
# inserts size factors (gene length and seq depth) into dds
dds_26 <- estimateSizeFactors(dds_26)
sizeFactors(dds_26)

#filter out low/non expressed genes (genes out where less than 3 samples have norm counts greater/equal to 10, before 55.773 genes, after  21.003 genes)
dim(dds_26)
idx <- rowSums(counts(dds_26, normalized=TRUE) >= 10) >= 3
dds_26 <- dds_26[idx,]
dim(dds_26)


### Normalized data counts
normalized_counts_26 <- counts(dds_26, normalized=TRUE)


###Transformation prior to QC
rld_26 <- rlog(dds_26, blind = TRUE)


### Hierarchical clustering
## Extract the rlog matrix from the object
# assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
rld_mat_26 <- assay(rld_26)
rld_mat_df_26 = as.data.frame(rld_mat_26)
## Compute pairwise corrrelation values
rld_cor_26 <- cor(rld_mat_df_26)
### Plot heatmap
ann_colors = list(grade2021 = c("1"="#9EBCDA", "2"="#8C6BB1", "3"="#810F7C"),
                  cluster = c("1"="cyan4", "2"="violetred4"),
                  MCconsensus = c("hypermetabolic"="forestgreen", "Immuneenriched"="red3",
                                  "Merlinintact"="royalblue2", "proliferative"="darkorange2"),
                  risk_score= c("high"="firebrick1", "intermediate"="purple1", "low"="dodgerblue2"))
meta_sub_26 = meta_26[, c("grade2021", "cluster", "MCconsensus", "risk_score")]

cairo_pdf(filename = "Correlation_plot_26.pdf",width = 12, height = 10)
pheatmap(rld_cor_26, annotation = meta_sub_26, annotation_colors = ann_colors)
dev.off()


####PCA analysis full dataset
#make PCA plot for any PC
pca_26 <- prcomp(t(rld_mat_26), scale. = TRUE)
df_pca_26 <- cbind(meta_sub_26, pca_26$x)
str(df_pca_26)

cairo_pdf(filename = "Screeplot_26.pdf",width = 16, height = 5)
fviz_screeplot(pca_26, addlabels=TRUE)
dev.off()

cairo_pdf(filename = "Screeplot_all_clean.pdf",width = 16, height = 5)
fviz_screeplot(pca_all)
dev.off()


#make PCA plots
cairo_pdf(filename = "PC1-2_cluster_26.pdf",width = 8, height = 8)
ggplot(df_pca_26, aes(x=PC1, y=PC2, color = cluster))+
  geom_point(size=7) +
  scale_color_manual(values=c("1"="cyan4", "2"="violetred4"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

cairo_pdf(filename = "PC1-3_cluster_26.pdf",width = 8, height = 8)
ggplot(df_pca_26, aes(x=PC1, y=PC3, color = cluster))+
  geom_point(size=7) +
  scale_color_manual(values=c("1"="cyan4", "2"="violetred4"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

cairo_pdf(filename = "PC1-4_cluster_26.pdf",width = 8, height = 8)
ggplot(df_pca_26, aes(x=PC1, y=PC4, color = cluster))+
  geom_point(size=7) +
  scale_color_manual(values=c("1"="cyan4", "2"="violetred4"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

cairo_pdf(filename = "PC1-5_cluster_26.pdf",width = 8, height = 8)
ggplot(df_pca_26, aes(x=PC1, y=PC5, color = cluster))+
  geom_point(size=7) +
  scale_color_manual(values=c("1"="cyan4", "2"="violetred4"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


### Running DESeq2 for differential gene expression
##Fill further add slots to dds object
dds_26 <- DESeq(dds_26)
## Coefficients available for testing
resultsNames(dds_26)

#get results
res_table_cluster2 <- results(dds_26, name ="cluster_2_vs_1", cooksCutoff = FALSE)
summary(res_table_cluster2,alpha=0.05)

res_cluster2_df <- data.frame(res_table_cluster2)
write.table(res_table_Val_2h, file="DEG_Val_2h_across_lines.txt")


## Set thresholds for more stringency
padj.cutoff <- 0.05
lfc.cutoff <- 1

threshold_Val_2h <- res_table_Val_2h$padj < padj.cutoff & abs(res_table_Val_2h$log2FoldChange) > lfc.cutoff
length(which(threshold_Val_2h))
res_table_Val_2h$threshold <- threshold_Val_2h
threshold_Val_72h <- res_table_Val_72h$padj < padj.cutoff & abs(res_table_Val_72h$log2FoldChange) > lfc.cutoff
length(which(threshold_Val_72h))
res_table_Val_72h$threshold <- threshold_Val_72h

#show changes in Volcano plot
LFCs_Val_2h <- data.frame(res_table_Val_2h)
LFCs_Val_2h$gene = row.names(LFCs_Val_2h)
length(which(LFCs_Val_2h$threshold == "TRUE" & LFCs_Val_2h$log2FoldChange > 1))
Val_2h_up = LFCs_Val_2h %>% filter(LFCs_Val_2h$threshold == "TRUE" & LFCs_Val_2h$log2FoldChange > 1)
Val_2h_up = Val_2h_up$gene
length(which(LFCs_Val_2h$threshold == "TRUE" & LFCs_Val_2h$log2FoldChange < 1))
Val_2h_down = LFCs_Val_2h %>% filter(LFCs_Val_2h$threshold == "TRUE" & LFCs_Val_2h$log2FoldChange < 1)
Val_2h_down = Val_2h_down$gene


#subset to 14 good quality samples

####get raw read counts data and meta data for all 14 samples
data_14 = read.csv(file="data_14.csv")
data_14 = data_26[,meta_14$QBIC_ID]
data_14$gene = rownames(data_14)
data_subset = data_14[,c("T26_0094", "T26_0097", "T26_0061", "T26_0016", "T26_0103", "T26_0076", "T26_0175", "T26_0072", "T26_0058", "T26_0123")]
rownames(data_subset) = data_14$X

#get meta data
meta_14 = meta_26[c(2,4,5,6,9,10,12,14,17,18,20,22,25,26),]
meta_subset = meta_14[colnames(data_subset),]
write.csv(meta_subset, file="meta_subset.csv")

#get raw count data
data_subset = read.csv(file = "raw_counts.csv")

#get meta data
meta_subset = read.csv(file="meta.csv")


#short code to get raw read counts for GEO
all(rownames(meta_subset) %in% colnames(data_subset))
all(rownames(meta_subset) == colnames(data_subset))
colnames(data_subset) = meta_subset$Lab_ID
write.csv(data_subset, file="raw_reads.csv")

## Create DESeq2Dataset object
dds_sub <- DESeqDataSetFromMatrix(countData = data_subset, colData = meta_subset, 
                                 design = ~cluster)
dim(dds_sub)

###Normalization
# inserts size factors (gene length and seq depth) into dds
dds_sub <- estimateSizeFactors(dds_sub)
sizeFactors(dds_sub)

#filter out low/non expressed genes (genes out where less than 3 samples have norm counts greater/equal to 10, before 55.773 genes, after  14.783 genes)
dim(dds_sub)
idx <- rowSums(counts(dds_sub, normalized=TRUE) >= 10) >= 3
dds_sub <- dds_sub[idx,]
dim(dds_sub)


### Normalized data counts
normalized_counts_sub <- counts(dds_sub, normalized=TRUE)


###Transformation prior to QC
rld_sub <- rlog(dds_sub, blind = TRUE)

# assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
rld_mat_sub <- assay(rld_sub)
rld_mat_df_sub = as.data.frame(rld_mat_sub)


####PCA analysis full dataset
#make PCA plot for any PC

pca_sub <- prcomp(t(rld_mat_sub), scale. = TRUE)
df_pca_sub <- cbind(meta_subset, pca_sub$x)
str(df_pca_sub)

cairo_pdf(filename = "Screeplot_14.pdf",width = 16, height = 5)
fviz_screeplot(pca_sub, addlabels=TRUE)
dev.off()


#make PCA plots
cairo_pdf(filename = "PC1-2_cluster_26.pdf",width = 8, height = 8)
ggplot(df_pca_sub, aes(x=PC1, y=PC2, color = cluster))+
  geom_point(size=7) +
  scale_color_manual(values=c("1"="cyan4", "2"="violetred4"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


### Running DESeq2 for differential gene expression
##Fill further add slots to dds object
dds_sub <- DESeq(dds_sub)
## Coefficients available for testing
resultsNames(dds_sub)

#get results
res_table_cluster2_sub <- results(dds_sub, name ="cluster_2_vs_1")
summary(res_table_cluster2_sub,alpha=0.05)

res_cluster2_df_sub <- data.frame(res_table_cluster2_sub)
res_cluster2_df_sub = read.csv(file = "res_cluster2_df_sub.csv")
write.csv(res_cluster2_df_sub, file="LFC_cluster2_vs_cluster1.csv")


##make volcano indicating protocadherin cluster
res_cluster2_df_sub$gene = res_cluster2_df_sub$X
PCDHA = rownames(res_cluster2_df_sub)[grep("^PCDHA", rownames(res_cluster2_df_sub))]
PCDHB = rownames(res_cluster2_df_sub)[grep("^PCDHB", rownames(res_cluster2_df_sub))]
PCDHG = rownames(res_cluster2_df_sub)[grep("^PCDHG", rownames(res_cluster2_df_sub))]

cairo_pdf(filename = "Volcano_cluster2vs1_hypermethylated_genes_4.pdf",width = 8, height = 6)
res_cluster2_df_sub %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point(color="grey30")+
  geom_point(data = res_cluster2_df_sub %>% filter(hypermethylated == "yes"), color = "deeppink2", size=5) +
  xlim(c(-10,10)) +
  ylim(c(0,10))+
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "grey", size=1.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()

cairo_pdf(filename = "Volcano_cluster2vs1_PCDH_cluster_4.pdf",width = 8, height = 6)
res_cluster2_df_sub %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point(color="grey30")+
  geom_point(data = res_cluster2_df_sub %>% filter(gene %in% PCDHA), color = "purple3", size=5) +
  geom_point(data = res_cluster2_df_sub %>% filter(gene %in% PCDHG), color = "red", size=5) +
  geom_point(data = res_cluster2_df_sub %>% filter(gene %in% PCDHB), color = "darkorange", size=5) +
  xlim(c(-10,10)) +
  ylim(c(0,10))+
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "grey", size=1.5)+
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+
  theme_classic()
dev.off()









################get data on n = 185 cases from Raleigh paper GSE183653

####get raw read counts data and meta data for all 26 samples
data_Raleigh = read.table("data/Raleigh_RNAseq.txt", header=T)

#convert Gene_IDs to symbol
library(org.Hs.eg.db)
library(AnnotationDbi)

mapped <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = as.character(data_Raleigh$Gene_ID),
                                columns = c("SYMBOL"),
                                keytype = "ENTREZID")
write.csv(mapped, file="ID_to_symbol.csv")

#re-load Raleigh data with gene symbols
data_Raleigh = read.table("data/Raleigh_RNAseq.txt", header=T)
data_Raleigh <- aggregate(data_Raleigh[-1], by = list(data_Raleigh$gene_symbol), FUN=sum)
colnames(data_Raleigh)[which(names(data_Raleigh) == "Group.1")] <- "gene_name"
rownames(data_Raleigh) = data_Raleigh$gene_name
data_Raleigh = data_Raleigh[,-1]

#get meta for Raleigh
meta_Raleigh = read.table("meta/Raleigh_meta.txt", header=T, row.names=1)
meta_Raleigh[] = lapply(meta_Raleigh, factor)
str(meta_Raleigh)

## Create DESeq2Dataset object
dds_Raleigh <- DESeqDataSetFromMatrix(countData = data_Raleigh, colData = meta_Raleigh, 
                                 design = ~Cluster)
dim(dds_Raleigh)

###Normalization
# inserts size factors (gene length and seq depth) into dds
dds_Raleigh <- estimateSizeFactors(dds_Raleigh)
sizeFactors(dds_Raleigh)

#filter out low/non expressed genes (genes out where less than 10 samples have norm counts greater/equal to 10, before 37,691 genes, after  22.831 genes)
dim(dds_Raleigh)
idx <- rowSums(counts(dds_Raleigh, normalized=TRUE) >= 10) >= 10
dds_Raleigh <- dds_Raleigh[idx,]
dim(dds_Raleigh)


### Normalized data counts
normalized_counts_Raleigh <- counts(dds_Raleigh, normalized=TRUE)


###Transformation prior to QC
vst_Raleigh <- vst(dds_Raleigh, blind = TRUE)

# assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
vst_mat_Raleigh <- assay(vst_Raleigh)

####PCA analysis full dataset
#make PCA plot for any PC

pca_Raleigh <- prcomp(t(vst_mat_Raleigh), scale. = TRUE)
df_pca_Raleigh <- cbind(meta_Raleigh, pca_Raleigh$x)
str(df_pca_Raleigh)

cairo_pdf(filename = "Screeplot_14.pdf",width = 16, height = 5)
fviz_screeplot(pca_Raleigh, addlabels=TRUE)
dev.off()

cairo_pdf(filename = "Screeplot_all_clean.pdf",width = 16, height = 5)
fviz_screeplot(pca_all)
dev.off()


#make PCA plots
cairo_pdf(filename = "PC1-2_cluster_26.pdf",width = 8, height = 8)
ggplot(df_pca_Raleigh, aes(x=PC1, y=PC2, color = grade))+
  geom_point(size=7) +
  theme_classic()
dev.off()



scale_color_manual(values=c("1"="cyan4", "2"="violetred4"))+

  theme(legend.position = "none")


### Running DESeq2 for differential gene expression
##Fill further add slots to dds object
dds_Raleigh <- DESeq(dds_Raleigh)
## Coefficients available for testing
resultsNames(dds_Raleigh)

#get results
res_table_cluster2_Raleigh <- results(dds_Raleigh, name ="Cluster_2_vs_1")
summary(res_table_cluster2_Raleigh,alpha=0.05)

res_cluster2_Raleigh_df <- data.frame(res_table_cluster2_Raleigh)
write.table(res_table_Val_2h, file="DEG_Val_2h_across_lines.txt")


