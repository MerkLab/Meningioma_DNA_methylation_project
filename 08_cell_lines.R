library(sesame)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(pheatmap)



###get data on validation cohort (longitudinal + non-recurrent, both EPICv2)
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/T26_discovery_celllines"
targets = read.csv(file="targets_discovery_celllines.csv")
betas.celllines = openSesame(idat_dir)
betas.celllines <- betas.celllines[complete.cases(betas.celllines), ]

# tag sex chromosome probes for removal
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

keep <- !(rownames(betas.celllines) %in% annoEPICv2$Name[annoEPICv2$chr %in% 
                                                         c("chrX","chrY")])
table(keep)
betas.celllines = betas.celllines[keep,]

#######make PCA analysis
#get top 10k most variable probes as for discovery
#make subset of most variable beta values
bVals.sub = betas.celllines
bVals.sub = as.data.frame(bVals.sub)
bVals.sub$var = apply(bVals.sub,1,var)
bVals.sub <- bVals.sub[order(bVals.sub$var, decreasing = TRUE),]
bVals.sub = bVals.sub[,-242]
bVals.sub.10k = bVals.sub[1:10000,]


###projecting cell lines in pca space from discovery cohort
###use PCA function from factormineR that contains active and supplementary individuals
###discovery samples are used as active individuals
###those determine the prinicipal components
###cell lines are used as supplementary components
meth.pca = PCA(t(bVals.sub.10k), ind.sup = 1:10,graph=FALSE)
fviz_eig(meth.pca, addlabels = TRUE)

coord.meth.pca = meth.pca$ind$coord
coord.meth.pca = as.data.frame(coord.meth.pca)
supp.ind.coord = as.data.frame(meth.pca$ind.sup$coord)
coord.meth.pca = rbind(supp.ind.coord, coord.meth.pca)
rownames(targets) = targets$ID
all(rownames(targets)==rownames(coord.meth.pca))


#make pca plots for various conditions for 10k probes
df <- cbind(targets, coord.meth.pca)
write.csv(df, file="PCA_discovery_celllines.csv")
df = read.csv(file="PCA_discovery_celllines.csv")

cairo_pdf(filename = "PCA_10k_disc_celllines_PC1-2.pdf", width = 5, height = 5)
ggplot(df, aes(x=Dim.1, y=Dim.2, color = clustering, size = 6, shape = clustering))+
  geom_jitter(data = df %>% filter(clustering == "METHlow"), shape=16, size=5, alpha=0.7,show.legend = F)+
  geom_jitter(data = df %>% filter(clustering == "METHhigh"), shape=16, size=5, alpha=0.7,show.legend = F)+
  geom_jitter(data = df %>% filter(clustering == "BEN_MEN"), shape=25, size=7,fill="black", alpha=0.6,show.legend = F)+
  geom_jitter(data = df %>% filter(clustering == "HBL52"), shape=17, size=7, alpha=0.6,show.legend = F)+
  geom_jitter(data = df %>% filter(clustering == "IOMM_LEE"), shape=19, size=7, alpha=0.6,show.legend = F)+  
  geom_jitter(data = df %>% filter(clustering == "KT21"), shape=18, size=9, alpha=0.6,show.legend = F)+
  geom_jitter(data = df %>% filter(clustering == "NCH93"), shape=15, size=7, alpha=0.6,show.legend = F)+
  scale_color_manual(values=c(BEN_MEN = "black",
                              IOMM_LEE = "black",
                              HBL52 = "black",
                              NCH93 = "black",
                              KT21 = "black",
                              METHhigh = "violetred4",
                              METHlow = "cyan4"))+
  scale_shape_manual(values = c(BEN_MEN = 0,
                                IOMM_LEE = 1,
                                HBL52 = 2,
                                NCH93 = 5,
                                KT21 = 6,
                                METHhigh = 16,
                                METHlow = 16)) +
  theme_classic()
dev.off()

cairo_pdf(filename = "PCA_10k_disc_celllines_PC1-3.pdf", width = 5, height = 5)
ggplot(df, aes(x=Dim.1, y=Dim.3, color = clustering, size = 6, shape = clustering))+
  geom_jitter(data = df %>% filter(clustering == "METHlow"), shape=16, size=5, alpha=0.7)+
  geom_jitter(data = df %>% filter(clustering == "METHhigh"), shape=16, size=5, alpha=0.7)+
  geom_jitter(data = df %>% filter(clustering == "BEN_MEN"), shape=25, size=7,fill="black", alpha=0.6)+
  geom_jitter(data = df %>% filter(clustering == "HBL52"), shape=17, size=7, alpha=0.6)+
  geom_jitter(data = df %>% filter(clustering == "IOMM_LEE"), shape=19, size=7, alpha=0.6)+  
  geom_jitter(data = df %>% filter(clustering == "KT21"), shape=18, size=9, alpha=0.6)+
  geom_jitter(data = df %>% filter(clustering == "NCH93"), shape=15, size=7, alpha=0.6)+
  scale_color_manual(values=c(BEN_MEN = "black",
                              IOMM_LEE = "black",
                              HBL52 = "black",
                              NCH93 = "black",
                              KT21 = "black",
                              METHhigh = "violetred4",
                              METHlow = "cyan4"))+
  scale_shape_manual(values = c(BEN_MEN = 0,
                                IOMM_LEE = 1,
                                HBL52 = 2,
                                NCH93 = 5,
                                KT21 = 6,
                                METHhigh = 16,
                                METHlow = 16)) +
  theme_classic()
dev.off()



####make boxplot to show hypermethylated METHhigh probes in discovery and cell lines
#get hyper probes only and subset beta values
hyper_probes = read.csv(file="probes_hyper_cluster_disc.csv", header = T)
hyper_probes = hyper_probes$Probe_ID
bVals = as.data.frame(betas.celllines)
bVals.hyper.probes = subset(bVals, rownames(bVals) %in% hyper_probes)
all(colnames(bVals.hyper.probes) == targets$Basename)
bVals.hyper.probes = bVals.hyper.probes[,match(targets$Basename, colnames(bVals.hyper.probes))]
all(colnames(bVals.hyper.probes) == targets$Basename)
colnames(bVals.hyper.probes) = targets$clustering
write.csv(bVals.hyper.probes, file="hyper_probes_beta_discovery_cellines.csv")


#do the plot, with only single averages for the cell lines as horizontal lines
avg_beta_hyper_probes_disc = read.csv(file="avg_hyper_probes.csv", header = T)

cairo_pdf(filename = "Hyper_probes_box_disc_cl_avg.pdf", width = 3.2, height = 7)
ggplot(avg_beta_hyper_probes_disc, 
       aes(x=factor(cluster, levels = c("METHlow", "METHhigh")), y=betas))+
  geom_point(position = position_jitter(), alpha=0.95, aes(color=betas))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red3"), limits=c(0,1))+
  geom_hline(yintercept=0.4808, linetype="dotted", color = "#A56D90", size=2)+
  geom_hline(yintercept=0.2751, linetype="dotted", color = "#1257BA", size=2)+
  geom_hline(yintercept=0.7595, linetype="dotted", color = "#F25050", size=2)+
  geom_hline(yintercept=0.7464, linetype="dotted", color = "#F45353", size=2)+
  geom_hline(yintercept=0.6986, linetype="dotted", color = "#FD6666", size=2)+
  theme_classic()+
  theme(legend.position = "none")
dev.off()



#visualize genomic regions and get probe IDs and their beta values
#after getting the probes for regions, re-do heatmap using pheatmap to have the same color range
all(colnames(bVals) == targets$Basename)
bVals = bVals[,match(targets$Basename, colnames(bVals))]
colnames(bVals) = targets$ID
bVals_cl = bVals[,1:10]

# extract cell line names (remove _A / _B)
cell_lines <- sub("_[AB]$", "", colnames(bVals_cl))

# average replicates
bVals_merged <- sapply(unique(cell_lines), function(cl) {
  rowMeans(bVals_cl[, cell_lines == cl, drop = FALSE])
})

# keep row names
bVals_merged <- as.data.frame(bVals_merged)
rownames(bVals_merged) <- rownames(bVals_cl)

#order from benign to malignant
bVals_merged = bVals_merged[,c(1,2,3,5,4)]

#PCDHA@ cluster
PCDHA_probes = visualizeRegion("chr5",140786136,141012344,bVals_merged, draw = F)
write.csv(PCDHA_probes,file ="Cell_lines_PCDHA_probes_betas.csv")

cairo_pdf(filename = "PCDHA_betas_cell_lines.pdf", width = 12, height = 3)
pheatmap(t(PCDHA_probes), color = colorRampPalette(c("navy","dodgerblue3","indianred1","red3"))(100),
         cluster_rows=FALSE, cluster_cols=FALSE,show_colnames = F, border_color = NA)
dev.off()


#PCDHB@ cluster
PCDHB_probes = visualizeRegion("chr5",141051394,141248234,bVals_merged, draw = F)
write.csv(PCDHB_probes,file ="Cell_lines_PCDHB_probes_betas.csv")

cairo_pdf(filename = "PCDHB_betas_cell_lines.pdf", width = 12, height = 3)
pheatmap(t(PCDHB_probes), color = colorRampPalette(c("navy","dodgerblue3","indianred1","red3"))(100),
         cluster_rows=FALSE, cluster_cols=FALSE,show_colnames = F, border_color = NA)
dev.off()

#PCDHG@ cluster
PCDHG_probes = visualizeRegion("chr5",141330685,141512979,bVals_merged, draw = F)
write.csv(PCDHG_probes,file ="Cell_lines_PCDHG_probes_betas.csv")

cairo_pdf(filename = "PCDHG_betas_cell_lines.pdf", width = 12, height = 3)
pheatmap(t(PCDHG_probes), color = colorRampPalette(c("navy","dodgerblue3","indianred1","red3"))(100),
         cluster_rows=FALSE, cluster_cols=FALSE,show_colnames = F, border_color = NA)
dev.off()



###boxplot for PCDH clusters
#get data
data_boxplot_PCDH = read.csv(file="data_longrange_cluster.csv")
str(data_boxplot_PCDH)
data_boxplot_PCDH$cell_line = factor(data_boxplot_PCDH$cell_line, levels = c("H","B","N","K","I"))

# grouped boxplot PCDH clusters
cairo_pdf(filename = "Boxplot_Celllines_PCDH_clusters.pdf", width = 5, height = 5)
ggplot(data_boxplot_PCDH, aes(x=cluster, y=values, fill=cell_line)) + 
  geom_boxplot(outliers = F)+
  scale_fill_manual(values=c("#1257BA","#A56D90","#FD6666","#F45353","#F25050"))+
  theme_classic()
dev.off()



save.image(file = "my_work_space.RData")
load("my_work_space.RData")


