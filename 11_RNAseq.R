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
library(ggrepel)


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
write.table(normalized_counts_Raleigh, file="Raleigh_norm_expression.txt", sep = "\t")


###Transformation prior to QC
vst_Raleigh <- vst(dds_Raleigh, blind = TRUE)

# assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
vst_mat_Raleigh <- assay(vst_Raleigh)

####PCA analysis full dataset
#make PCA plot for any PC

pca_Raleigh <- prcomp(t(vst_mat_Raleigh), scale. = TRUE)
df_pca_Raleigh <- cbind(meta_Raleigh, pca_Raleigh$x)
str(df_pca_Raleigh)

cairo_pdf(filename = "Screeplot_Raleigh.pdf",width = 16, height = 5)
fviz_screeplot(pca_Raleigh, addlabels=TRUE)
dev.off()


#make PCA plots
cairo_pdf(filename = "PC1-2_cluster_grade.pdf",width = 8, height = 8)
ggplot(df_pca_Raleigh, aes(x=PC1, y=PC2, color = grade))+
  geom_point(size=7) +
  theme_classic()
dev.off()

cairo_pdf(filename = "PC1-2_cluster_grade.pdf",width = 8, height = 8)
ggplot(df_pca_Raleigh, aes(x=PC1, y=PC2, color = Cluster))+
  geom_point(size=7) +
  scale_color_manual(values=c("1"="cyan4", "2"="violetred4"))+
  theme_classic()
dev.off()


### Running DESeq2 for differential gene expression
##Fill further add slots to dds object
dds_Raleigh <- DESeq(dds_Raleigh)
## Coefficients available for testing
resultsNames(dds_Raleigh)

#get results
res_table_cluster2_Raleigh <- results(dds_Raleigh, name ="Cluster_2_vs_1")
summary(res_table_cluster2_Raleigh,alpha=0.05)

res_cluster2_Raleigh_df <- data.frame(res_table_cluster2_Raleigh)
write.table(res_cluster2_Raleigh_df, file="Raleigh_cluster_2_vs_1.txt", sep = "\t")


#make the plots
Raleigh_FC = read.csv(file="Raleigh_cluster_2_vs_1.csv")
names(Raleigh_FC)[names(Raleigh_FC) == "X"] <- "gene"

#all hypermethylted genes
cairo_pdf(filename = "Volcano_Raleigh_cluster2vs1_hypermethylated_genes.pdf",width = 8, height = 6)
Raleigh_FC %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point(color="grey30")+
  geom_point(data = Raleigh_FC %>% filter(hypermethylated == "yes"), color = "deeppink2", size=5) +
  xlim(c(-6,8)) +
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

#help to find genes
Raleigh_FC %>% filter(gene == "PITX1") %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj),label=gene))+
  geom_point(color="grey30")+
  geom_text()+
  xlim(c(-6,8)) +
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

#protocadherin genes


cairo_pdf(filename = "Volcano_Raleigh_PCDH_clusters.pdf",width = 8, height = 6)
Raleigh_FC %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point(color="grey30")+
  geom_point(data = Raleigh_FC %>% filter(gene %in% PCDHA), color = "purple3", size=5) +
  geom_point(data = Raleigh_FC %>% filter(gene %in% PCDHG), color = "red", size=5) +
  geom_point(data = Raleigh_FC %>% filter(gene %in% PCDHB), color = "darkorange", size=5) +
  xlim(c(-6,8)) +
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

#help to find genes
Raleigh_FC %>% filter(gene %in% PCDHB) %>%
  ggplot(aes(x=log2FoldChange, y= -log10(padj),label=gene))+
  geom_point(color="grey30")+
  geom_text()+
  xlim(c(-6,8)) +
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


#make correlation plots individually for each PCDH cluster

#PCDHA
df_PCDHA = read.csv(file="Raleigh_PCDHA.csv")

cairo_pdf(filename = "Scatter_Raleigh_PCDHA.pdf",width = 8, height = 6)
ggplot(df_PCDHA, aes(x = AVG_clust1, y = AVG_clust2, label = gene)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_text_repel(size = 5, fontface = "italic") +
  labs(
    x = "Average expression in cluster 1",
    y = "Average expression in cluster 2",
    title = "Scatter plot of cluster averages with regression line"
  ) +
  theme_minimal()
dev.off()

#PCDHB
df_PCDHB = read.csv(file="Raleigh_PCDHB.csv")

cairo_pdf(filename = "Scatter_Raleigh_PCDHB.pdf",width = 8, height = 6)
ggplot(df_PCDHB, aes(x = AVG_clust1, y = AVG_clust2, label = gene)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_text_repel(size = 5, fontface = "italic") +
  labs(
    x = "Average expression in cluster 1",
    y = "Average expression in cluster 2",
    title = "Scatter plot of cluster averages with regression line"
  ) +
  theme_minimal()
dev.off()




#PCDHG
df_PCDHG = read.csv(file="Raleigh_PCDHG.csv")

cairo_pdf(filename = "Scatter_Raleigh_PCDHG.pdf",width = 8, height = 6)
ggplot(df_PCDHG, aes(x = AVG_clust1, y = AVG_clust2, label = gene)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  geom_text_repel(size = 5, fontface = "italic") +
  labs(
    x = "Average expression in cluster 1",
    y = "Average expression in cluster 2",
    title = "Scatter plot of cluster averages with regression line"
  ) +
  theme_minimal()
dev.off()







