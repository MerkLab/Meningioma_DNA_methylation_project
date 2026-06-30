library(sesame)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(FactoMineR)
library(factoextra)
library(dplyr)
library(ggplot2)
library(here)

###clone figshare repository including the given folder structure to your local directory
###populate raw data subfolders with data from either GEO or figshare itself
###refer to the README.txt file at https://github.com/MerkLab/Meningioma_DNA_methylation_project for further information


###get data (discovery cohort and cell lines)
idat_dir = here("data", "GSE304093")
setwd(here("data", "general_datasets"))
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
bVals.sub.50k = bVals.sub[1:10000,]


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

setwd(here("results"))
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

setwd(here("results"))
cairo_pdf(filename = "PCA_10k_disc_celllines_PC1-3.pdf", width = 5, height = 5)
ggplot(df, aes(x=Dim.1, y=Dim.3, color = clustering, size = 6, shape = clustering))+
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


####make boxplot to show hypermethylated METHhigh probes in discovery and cell lines
#get hyper probes only and subset beta values
setwd(here("data", "08_cell_lines"))
hyper_probes = read.csv(file="probes_hyper_cluster_disc.csv", header = T)
hyper_probes = hyper_probes$Probe_ID
bVals = as.data.frame(betas.celllines)
bVals.hyper.probes = subset(bVals, rownames(bVals) %in% hyper_probes)
all(colnames(bVals.hyper.probes) == targets$Basename)
bVals.hyper.probes = bVals.hyper.probes[,match(targets$Basename, colnames(bVals.hyper.probes))]
all(colnames(bVals.hyper.probes) == targets$Basename)
colnames(bVals.hyper.probes) = targets$clustering
setwd(here("results"))
write.csv(bVals.hyper.probes, file="hyper_probes_beta_discovery_cellines.csv")


#do the plot, with only single averages for the cell lines as horizontal lines
setwd(here("data", "08_cell_lines"))
avg_beta_hyper_probes_disc = read.csv(file="avg_hyper_probes.csv", header = T)

setwd(here("results"))
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
setwd(here("results"))
write.csv(PCDHA_probes,file ="Cell_lines_PCDHA_probes_betas.csv")

setwd(here("results"))
cairo_pdf(filename = "PCDHA_betas_cell_lines.pdf", width = 12, height = 3)
pheatmap(t(PCDHA_probes), color = colorRampPalette(c("navy","dodgerblue3","indianred1","red3"))(100),
         cluster_rows=FALSE, cluster_cols=FALSE,show_colnames = F, border_color = NA)
dev.off()


#PCDHB@ cluster
PCDHB_probes = visualizeRegion("chr5",141051394,141248234,bVals_merged, draw = F)
setwd(here("results"))
write.csv(PCDHB_probes,file ="Cell_lines_PCDHB_probes_betas.csv")

setwd(here("results"))
cairo_pdf(filename = "PCDHB_betas_cell_lines.pdf", width = 12, height = 3)
pheatmap(t(PCDHB_probes), color = colorRampPalette(c("navy","dodgerblue3","indianred1","red3"))(100),
         cluster_rows=FALSE, cluster_cols=FALSE,show_colnames = F, border_color = NA)
dev.off()

#PCDHG@ cluster
PCDHG_probes = visualizeRegion("chr5",141330685,141512979,bVals_merged, draw = F)
setwd(here("results"))
write.csv(PCDHG_probes,file ="Cell_lines_PCDHG_probes_betas.csv")

setwd(here("results"))
cairo_pdf(filename = "PCDHG_betas_cell_lines.pdf", width = 12, height = 3)
pheatmap(t(PCDHG_probes), color = colorRampPalette(c("navy","dodgerblue3","indianred1","red3"))(100),
         cluster_rows=FALSE, cluster_cols=FALSE,show_colnames = F, border_color = NA)
dev.off()



###boxplot for PCDH clusters
#get data
setwd(here("data", "08_cell_lines"))
data_boxplot_PCDH = read.csv(file="data_longrange_cluster.csv")
str(data_boxplot_PCDH)
data_boxplot_PCDH$cell_line = factor(data_boxplot_PCDH$cell_line, levels = c("H","B","N","K","I"))

# grouped boxplot PCDH clusters
setwd(here("results"))
cairo_pdf(filename = "Boxplot_Celllines_PCDH_clusters.pdf", width = 5, height = 5)
ggplot(data_boxplot_PCDH, aes(x=cluster, y=values, fill=cell_line)) + 
  geom_boxplot(outliers = F)+
  scale_fill_manual(values=c("#1257BA","#A56D90","#FD6666","#F45353","#F25050"))+
  theme_classic()
dev.off()




###work on FACS data
###parental meningioma cells MFI and percentage b-catenin cyto + nuclear

#get MFI data for parental cells
setwd(here("data", "08_cell_lines"))
parental_MFI = read.csv(file="Parental_MFI.csv", header = T)
parental_MFI$line = factor(parental_MFI$line, levels = c("HBL", "BEN", "NCH", "Lee", "KT"))
parental_MFI$group = factor(parental_MFI$group, levels = c("cytplasmic", "nuclear"))
str(parental_MFI)

#plot with datapoints
# Define the manual order of your samples
desired_order <- c("HBL", "BEN", "NCH", "Lee", "KT")  # replace with your order

# Set factor levels manually
parental_MFI$line <- factor(parental_MFI$line, levels = desired_order)

# Create numeric x-axis for spacing
parental_MFI$x_pos <- as.numeric(parental_MFI$line) * 1.5  # increase multiplier for more spacing

# Define consistent dodge for groups
dodge <- position_dodge(width = 0.6)

# Plot
setwd(here("results"))
cairo_pdf(filename = "Parental_MFI_points.pdf", width = 3, height = 7)
ggplot(parental_MFI, aes(x = x_pos, y = MFI_FC, color = group)) +
  # Raw points with slight jitter
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = dodge$width
    ),
    size = 4,
    alpha = 0.8
  ) +
  # Error bars on top of mean
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = dodge,
    width = 0.5,
    linewidth = 1.4
  ) +
  # Mean points
  geom_point(
    stat = "summary",
    fun = mean,
    position = dodge,
    size = 3
  ) +
  scale_color_manual(values = c("magenta2")) +
  
  #add horizontal line
  geom_hline(yintercept = 1, linetype = "dashed", color = "deepskyblue1", size=2)+
  
  # Map numeric x back to labels with manual order
  scale_x_continuous(
    breaks = 1:length(desired_order) * 1.5,
    labels = desired_order
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()


#get percentage data for parental cells
setwd(here("data", "08_cell_lines"))
parental_percent = read.csv(file="Parental_percentage.csv", header = T)
parental_percent$line = factor(parental_percent$line, levels = c("HBL", "BEN", "NCH", "Lee", "KT"))
parental_percent$group = factor(parental_percent$group, levels = c("cytplasmic", "nuclear"))
str(parental_percent)

#plot with datapoints shown
desired_order <- c("HBL", "BEN", "NCH", "Lee", "KT")

# Set factor levels manually
parental_percent$line <- factor(parental_percent$line, levels = desired_order)

# Create numeric x-axis for spacing
parental_percent$x_pos <- as.numeric(parental_percent$line) * 1.5

# Define dodge for groups
dodge <- position_dodge(width = 0.6)

setwd(here("results"))
cairo_pdf(filename = "Parental_percent_points.pdf", width = 4, height = 7)
ggplot(parental_percent, aes(x = x_pos, y = MFI_FC, color = group)) +
  # Raw points with slight jitter
  geom_point(
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = dodge$width),
    size = 4,
    alpha = 0.8
  ) +
  # Error bars
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = dodge,
    width = 0.5,
    linewidth = 1.4
  ) +
  # Mean points
  geom_point(
    stat = "summary",
    fun = mean,
    position = dodge,
    size = 3
  ) +
  scale_color_manual(values = c("deepskyblue1", "magenta2")) +
  # Map numeric x back to labels in the correct order
  scale_x_continuous(
    breaks = 1:length(desired_order) * 1.5,   # <- fixed
    labels = desired_order                     # <- fixed
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()


###investigate changes of PCDHGC3 overexpressiion on b-catenin localisation

#####BEN-MEN-1 analysis with GFP or PCDHGC3 overexpression
#get MFI data
#get MFI data for GFP cells
setwd(here("data", "08_cell_lines"))
GFP_MFI = read.csv(file="GFP_MFI_BEN.csv", header = T)
GFP_MFI$line = factor(GFP_MFI$line, levels = c("BEN_parental","BEN_GFP", "BEN_C3"))
GFP_MFI$group = factor(GFP_MFI$group, levels = c("cytoplasmic", "nuclear"))
str(GFP_MFI)

#make plot with datapoints
# Define manual order of samples
desired_order <- c("BEN_parental", "BEN_GFP", "BEN_C3") 

# Set factor levels manually
GFP_MFI$line <- factor(GFP_MFI$line, levels = desired_order)

# Create numeric x-axis for spacing
GFP_MFI$x_pos <- as.numeric(GFP_MFI$line) * 1.5  # adjust multiplier for more spacing

# Define consistent dodge for groups
dodge <- position_dodge(width = 0.6)

# Plot
setwd(here("results"))
cairo_pdf(filename = "GFP_MFI_points_BEN.pdf", width = 4, height = 5)
ggplot(GFP_MFI, aes(x = x_pos, y = MFI_FC, color = group)) +
  # Raw data points with slight jitter
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = dodge$width
    ),
    size = 6,
    alpha = 0.8
  ) +
  # Error bars for mean ± SE
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = dodge,
    width = 0.5,
    linewidth = 1.4
  ) +
  # Mean points
  geom_point(
    stat = "summary",
    fun = mean,
    position = dodge,
    size = 3
  ) +
  scale_color_manual(values = c("magenta2")) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 2.5, by = 0.5)) +
  
  #add horizontal line
  geom_hline(yintercept = 1, linetype = "dashed", color = "deepskyblue1", size=2)+
  
  # Map numeric x back to labels in the manual order
  scale_x_continuous(
    breaks = 1:length(desired_order) * 1.5,
    labels = desired_order
  ) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()


#####IOMM-Lee analysis with GFP or PCDHGC3 overexpression
#get MFI data
#get MFI data for GFP cells
setwd(here("data", "08_cell_lines"))
GFP_MFI = read.csv(file="GFP_MFI_Lee.csv", header = T)
GFP_MFI$line = factor(GFP_MFI$line, levels = c("Lee_parental","Lee_GFP", "Lee_C3"))
GFP_MFI$group = factor(GFP_MFI$group, levels = c("cytplasmic", "nuclear"))
str(GFP_MFI)

#make plot with datapoints
# Define manual order of samples
desired_order <- c("Lee_parental", "Lee_GFP", "Lee_C3") 

# Set factor levels manually
GFP_MFI$line <- factor(GFP_MFI$line, levels = desired_order)

# Create numeric x-axis for spacing
GFP_MFI$x_pos <- as.numeric(GFP_MFI$line) * 1.5  # adjust multiplier for more spacing

# Define consistent dodge for groups
dodge <- position_dodge(width = 0.6)

# Plot
setwd(here("results"))
cairo_pdf(filename = "GFP_MFI_points.pdf", width = 5, height = 5)
ggplot(GFP_MFI, aes(x = x_pos, y = MFI_FC, color = group)) +
  # Raw data points with slight jitter
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.12,
      dodge.width = dodge$width
    ),
    size = 6,
    alpha = 0.8
  ) +
  # Error bars for mean ± SE
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = dodge,
    width = 0.5,
    linewidth = 1.4
  ) +
  # Mean points
  geom_point(
    stat = "summary",
    fun = mean,
    position = dodge,
    size = 3
  ) +
  scale_color_manual(values = c("magenta2")) +
  scale_y_continuous(limits = c(0, NA),
                     breaks = seq(0, 2.5, by = 0.5)) +
  
  #add horizontal line
  geom_hline(yintercept = 1, linetype = "dashed",color = "deepskyblue1", size=2)+
  
  # Map numeric x back to labels in the manual order
  scale_x_continuous(
    breaks = 1:length(desired_order) * 1.5,
    labels = desired_order
  ) +
  theme_classic() +
  theme(legend.position = "none")
dev.off()


###analyses b-catenin cytoplasmic and nuclear from ICC stains

##############make boxplot of total b-catenin intensity
###IOMM-LEE
setwd(here("data", "08_cell_lines"))
Lee_total = read.csv(file="Bcatenin_ICC_quartiles_total_Lee.csv")
Lee_total$MainGroup <- factor(Lee_total$MainGroup, levels = c("Lee-GFP", "Lee-C3"))

setwd(here("results"))
cairo_pdf(filename = "Boxplots_ICC_Bcatenin_total_3.pdf", width = 5, height = 7)
ggplot(Lee_total, aes(x = MainGroup, y = Value, fill=MainGroup)) +
  geom_jitter(aes(color = MainGroup),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.8, alpha = 0.6) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("Lee-GFP" = "violetred4", "Lee-C3" = "steelblue4")) +
  scale_color_manual(values = c("Lee-GFP" = "violetred4", "Lee-C3" = "steelblue4")) +
  labs(title = "Boxplot title", x = "Main Group", y = "Value") +
  theme_classic()
dev.off()


##############make boxplot of b-catenin ratio nucleus to cytoplasm
setwd(here("data", "08_cell_lines"))
Lee_ratio = read.csv(file="Bcatenin_ICC_quartiles_ratios_Lee.csv")
Lee_ratio$MainGroup <- factor(Lee_ratio$MainGroup, levels = c("Lee-GFP", "Lee-C3"))

setwd(here("results"))
cairo_pdf(filename = "Boxplots_ICC_Bcatenin_ratios_new.pdf", width = 7, height = 7)
ggplot(Lee_ratio, aes(x = MainGroup, y = Value, fill=SubGroup)) +
  geom_jitter(aes(color = SubGroup),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.8, alpha = 0.6) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("Q1" = "cyan2", "Q2" = "turquoise4", 
                               "Q3" = "violetred2", "Q4" = "darkmagenta")) +
  scale_color_manual(values = c("Q1" = "cyan2", "Q2" = "turquoise4", 
                                "Q3" = "violetred2", "Q4" = "darkmagenta")) +
  geom_hline(yintercept = 0.88, linetype = "dashed", color = "cyan2",linewidth = 1.25)+
  geom_hline(yintercept = 1.13, linetype = "dashed", color = "turquoise4",linewidth = 1.25)+
  geom_hline(yintercept = 1.46, linetype = "dashed", color = "violetred2",linewidth = 1.25)+
  geom_hline(yintercept = 2.15, linetype = "dashed", color = "darkmagenta",linewidth = 1.25)+
  labs(title = "Boxplot title", x = "Main Group", y = "Value") +
  theme_classic()
dev.off()

wilcox.test(Value ~ MainGroup, data = Lee_total)
wilcox.test(Value ~ MainGroup, data = Lee_ratio)



#do tha same for BEN-MEN-1
setwd(here("data", "08_cell_lines"))
BEN_total = read.csv(file="Bcatenin_ICC_quartiles_total_BEN.csv")
BEN_total$MainGroup <- factor(BEN_total$MainGroup, levels = c("BEN_GFP", "BEN_C3"))

setwd(here("results"))
cairo_pdf(filename = "Boxplots_ICC_Bcatenin_total_BEN.pdf", width = 5, height = 7)
ggplot(BEN_total, aes(x = MainGroup, y = Value, fill=MainGroup)) +
  geom_jitter(aes(color = MainGroup),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.8, alpha = 0.6) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("BEN_GFP" = "violetred4", "BEN_C3" = "steelblue4")) +
  scale_color_manual(values = c("BEN_GFP" = "violetred4", "BEN_C3" = "steelblue4")) +
  labs(title = "Boxplot title", x = "Main Group", y = "Value") +
  theme_classic()
dev.off()

##############make boxplot of b-catenin ratio nucleus to cytoplasm
setwd(here("data", "08_cell_lines"))
BEN_ratio = read.csv(file="Bcatenin_ICC_quartiles_ratios_BEN.csv")
BEN_ratio$MainGroup <- factor(BEN_ratio$MainGroup, levels = c("BEN_GFP", "BEN_C3"))

setwd(here("results"))
cairo_pdf(filename = "Boxplots_ICC_Bcatenin_ratios_BEN.pdf", width = 7, height = 7)
ggplot(BEN_ratio, aes(x = MainGroup, y = Value, fill=SubGroup)) +
  geom_jitter(aes(color = SubGroup),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              size = 1.8, alpha = 0.6) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("Q1" = "cyan2", "Q2" = "turquoise4", 
                               "Q3" = "violetred2", "Q4" = "darkmagenta")) +
  scale_color_manual(values = c("Q1" = "cyan2", "Q2" = "turquoise4", 
                                "Q3" = "violetred2", "Q4" = "darkmagenta")) +
  geom_hline(yintercept = 0.12538226, linetype = "dashed", color = "cyan2",linewidth = 1.25)+
  geom_hline(yintercept = 0.37512742, linetype = "dashed", color = "turquoise4",linewidth = 1.25)+
  geom_hline(yintercept = 0.62487258, linetype = "dashed", color = "violetred2",linewidth = 1.25)+
  geom_hline(yintercept = 0.87512742, linetype = "dashed", color = "darkmagenta",linewidth = 1.25)+
  labs(title = "Boxplot title", x = "Main Group", y = "Value") +
  theme_classic()
dev.off()

wilcox.test(Value ~ MainGroup, data = BEN_total)
wilcox.test(Value ~ MainGroup, data = BEN_ratio)

