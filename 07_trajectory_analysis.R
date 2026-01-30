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
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(grDevices)
library(lme4)
library(lmerTest) 
library(ordinal)


######make Sankey plot for longitudinal cohort
Disc_Sankey = read.csv(file="Disc_For_Sankey.csv")

node_order <- c(
  # Setting
  "Primary_no_rec", "Primary_rec", "Rec_1", "Rec_2",
  # WHO grades
  "1", "2", "3",
  # Toronto (desired order)
  "Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"
)

nodes <- data.frame(
  name = node_order,
  stringsAsFactors = FALSE
)

nodes$label <- ""

nodes$layer <- NA_integer_

nodes$layer[nodes$name %in% c("Primary_no_rec", "Primary_rec", "Rec_1", "Rec_2")] <- 0
nodes$layer[nodes$name %in% c("1", "2", "3")] <- 1
nodes$layer[nodes$name %in% c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative")] <- 2

nodes <- nodes %>%
  arrange(layer, factor(name, levels = node_order))



links_1 <- Disc_Sankey %>%
  group_by(Setting, WHO_grade) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(Setting, nodes$name) - 1,
    target = match(WHO_grade, nodes$name) - 1
  ) %>%
  select(source, target, value)

links_2 <- Disc_Sankey %>%
  group_by(WHO_grade, Toronto) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(WHO_grade, nodes$name) - 1,
    target = match(Toronto, nodes$name) - 1
  ) %>%
  select(source, target, value)

links <- bind_rows(links_1, links_2)

# Setting
nodes$group <- NA_character_
nodes$group[nodes$name == "Primary_no_rec"]     <- "Pri_no"
nodes$group[nodes$name == "Primary_rec"]   <- "Pri_yes"
nodes$group[nodes$name == "Rec_1"]   <- "Rec1"
nodes$group[nodes$name == "Rec_2"]    <- "Rec2"

# WHO grades
nodes$group[nodes$name == "1"] <- "Grade1"
nodes$group[nodes$name == "2"] <- "Grade2"
nodes$group[nodes$name == "3"] <- "Grade3"

# Toronto
nodes$group[nodes$name == "Merlinintact"]     <- "Merlin"
nodes$group[nodes$name == "Immuneenriched"]   <- "Immune"
nodes$group[nodes$name == "hypermetabolic"]   <- "Hyper"
nodes$group[nodes$name == "proliferative"]    <- "Prolif"

stopifnot(!any(is.na(nodes$group)))



my_colour_scale <- '
d3.scaleOrdinal()
  .domain([
    "Grade1","Grade2","Grade3",
    "Merlin","Immune","Hyper","Prolif",
    "Pri_no", "Pri_yes", "Rec1", "Rec2"
  ])
  .range([
    "#9EBCDA", "#8C6BB1", "#810F7C",
    "#436eee", "#cd0000", "#228b22", "#ee7600",
    "#1CFDB2", "#006400", "#FFCB1B", "#D95F02"
  ])
'

links$linkGroup <- nodes$group[links$source + 1]


stopifnot("linkGroup" %in% colnames(links))


sankeyNetwork(
  Links = as.data.frame(links),
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value  = "value",
  NodeID = "name",
  
  NodeGroup = "group",
  LinkGroup = "linkGroup",
  
  colourScale = my_colour_scale,
  fontSize = 14,
  nodeWidth = 50,
  sinksRight = FALSE
)


sankey = sankeyNetwork(
  Links = as.data.frame(links),
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value  = "value",
  NodeID = "label",
  
  NodeGroup = "group",
  LinkGroup = "linkGroup",
  
  colourScale = my_colour_scale,
  fontSize = 14,
  nodeWidth = 50,
  sinksRight = FALSE
)

saveWidget(sankey, "sankey.html", selfcontained = TRUE)
webshot(
  "sankey.html",
  file = "sankey_discovery.pdf",
  vwidth = 1600,
  vheight = 1000
)





dat_mol_group <- data.frame(
  c(11,2,1,0),
  c(1,0,0,0),
  c(1,8,9,3),
  c(0,8,8,3),
  row.names = c("Pri_no", "Pri_yes", "Rec_1", "Rec_2"),
  stringsAsFactors = FALSE
)
colnames(dat_mol_group) <- c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative")
dat_mol_group

x <- c()
for (row in rownames(dat_mol_group)) {
  for (col in colnames(dat_mol_group)) {
    x <- rbind(x, matrix(rep(c(row, col), dat_mol_group[row, col]), ncol = 2, byrow = TRUE))
  }
}
df_mol_group <- as.data.frame(x)
colnames(df_mol_group) <- c("setting", "MCconsensus")
df_mol_group
df_mol_group$MCconsensus = factor(df_mol_group$MCconsensus, levels =  c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

test_mol_group <- chisq.test(table(df_mol_group))

tiff("Chisquare_validation_mol_group.tiff", res=300,width=1000,height=2048)
ggbarstats(df_mol_group, MCconsensus, setting,
           results.subtitle = FALSE,
           subtitle = paste0(test_mol_group$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("darkorange2","forestgreen","red3","royalblue2"))+
  theme(legend.position = "none")
dev.off()

#make the test with repeated measures
write.csv(df_mol_group, file="data_longi_group.csv")
test_data_group = read.csv(file="data_longi_group.csv", header = T)


test_data_group$MCconsensus_ord <- factor(
  test_data_group$MCconsensus,
  levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"),
  ordered = TRUE
)
test_data_group$setting <- factor(
  test_data_group$setting,
  levels = c("Pri_no", "Pri_yes", "Rec_1", "Rec_2"),
  ordered = TRUE 
)


model <- clmm(
  MCconsensus_ord ~ setting + (1 | patient_ID),
  data = test_data_group
)

# LRT against null
model_null <- clmm(MCconsensus_ord ~ 1 + (1 | patient_ID), data = test_data_group)
anova(model_null, model)





#for WHO grade
library(ggstatsplot)

dat_mol_grade <- data.frame(
  c(10,5,0,0),
  c(3,13,12,5),
  c(0,0,6,1),
  row.names = c("Pri_no", "Pri_yes", "Rec_1", "Rec_2"),
  stringsAsFactors = FALSE
)
colnames(dat_mol_grade) <- c("grade1", "grade2", "grade3")
dat_mol_grade


x <- c()
for (row in rownames(dat_mol_grade)) {
  for (col in colnames(dat_mol_grade)) {
    x <- rbind(x, matrix(rep(c(row, col), dat_mol_grade[row, col]), ncol = 2, byrow = TRUE))
  }
}
df_mol_grade <- as.data.frame(x)
colnames(df_mol_grade) <- c("setting", "grade")
df_mol_grade
df_mol_grade$grade = factor(df_mol_grade$grade, levels =  c("grade1", "grade2", "grade3"))

test_mol_grade <- chisq.test(table(df_mol_grade))

tiff("Chisquare_validation_grade.tiff", res=300,width=1000,height=2048)
ggbarstats(df_mol_grade, grade, setting,
           results.subtitle = FALSE,
           subtitle = paste0(test_mol_grade$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#810F7C", "#8C6BB1","#9EBCDA"))+
  theme(legend.position = "none")
dev.off()

#make the test with repeated measures
write.csv(df_mol_grade, file="data_longi_grade.csv")
test_data_grade = read.csv(file="data_longi_grade.csv", header = T)
test_data_grade$grade <- factor(
  test_data_grade$grade,
  levels = c("grade1", "grade2", "grade3"),
  ordered = TRUE
)

test_data_grade$setting <- factor(
  test_data_grade$setting,
  levels = c("Pri_no", "Pri_yes", "Rec_1", "Rec_2"),
  ordered = TRUE
)


model <- clmm(
  grade ~ setting + (1 | patient_ID),
  data = test_data_grade
)

summary(model)
model_null <- clmm(
  grade ~ 1 + (1 | patient_ID),
  data = test_data_grade
)

anova(model_null, model)













####get data for normal meninges, longitudinal cohort, and non-recurrent controls
####combine only probes present on EPIC

###get data on validation cohort (longitudinal + non-recurrent, both EPICv2)
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/T26_validation"
meta_val = read.csv(file="targets_validation_with_predictions.csv")

beta.val = openSesame(idat_dir)
colnames(beta.val)
beta.val = beta.val[,match(meta_val$Basename, colnames(beta.val))]
all(colnames(beta.val) %in% meta_val$Basename)
colnames(beta.val) = meta_val$ID

#remove sex chromosome probes masked probes
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
keep <- !(rownames(beta.val) %in% annoEPICv2$Name[annoEPICv2$chr %in% 
                                                          c("chrX","chrY")])
table(keep)
dim(beta.val)
beta.val = beta.val[keep,]
dim(beta.val) 
beta.val = beta.val[complete.cases(beta.val),]
#552.488 probes left

#convert to EPIC probes
beta.val = mLiftOver(beta.val, "EPIC")
dim(beta.val)
rownames(beta.val)
beta.val = beta.val[complete.cases(beta.val),]
#434.789 probes left EPIC compatible



#get normal meninges data
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/IDAT_normal_meninges"
meta_meninges <- read.metharray.sheet(idat_dir, pattern="targets_normal_meninges.csv")

beta.men = openSesame(idat_dir)
colnames(beta.men)
all(colnames(beta.men) %in% meta_meninges$Basename)
colnames(beta.men) = meta_meninges$ID
dim(beta.men)
beta.men = beta.men[complete.cases(beta.men),]
dim(beta.men)

#reduce validation and meninges beta values to same rownames, bring to same order, and combine
common_rows = intersect(rownames(beta.val), rownames(beta.men))
#428.961 probes common to both
beta.val.common = beta.val[common_rows,,drop=FALSE]
beta.men.common = beta.men[common_rows,,drop=FALSE]
all(rownames(beta.val.common) %in% rownames(beta.men.common))
all(rownames(beta.val.common) == rownames(beta.men.common))

beta.combine = cbind(beta.men.common,beta.val.common)
dim(beta.combine)



###check average METH cluster signature activity in meninges and MNG, avg per sample
###methylation cluster signature defined as top 1000 probes from PC1 used fo clustering
clust_sig = read.csv(file="top1000_PC1_EPICv2.csv", header = T)
clust_sig = clust_sig$x
clust_sig = mLiftOver(clust_sig, "EPIC")

#subset combined bVal data by top2000 probes
beta.combine = as.data.frame(beta.combine)
beta.combine.top = beta.combine %>% dplyr::filter(rownames(beta.combine) %in% clust_sig)
#calculate avg per sample
sample_avg_top1000 = apply(beta.combine.top,2,mean)
#make dataframe for plot
df_sample_signature_cluster = as.data.frame(sample_avg_top1000)
meta_temp = meta_val[,-3]
names(meta_temp)[names(meta_temp) == "Mcconsensus"] <- "MCconsensus"
all(colnames(meta_meninges) == colnames(meta_temp))
meta_combine = rbind(meta_meninges, meta_temp)
all(names(sample_avg_top1000) == meta_combine$ID)
df_sample_signature_cluster$setting_detail = meta_combine$setting_detail

cairo_pdf(filename = "Boxplot_cluster_probes_avg_combined_data.pdf", width = 8, height = 6)
ggplot(df_sample_signature_cluster, 
       aes(setting_detail, sample_avg_top1000))+
  geom_point(position = position_jitter(width=0.1), alpha=0.95, aes(color=sample_avg_top1000), size=6)+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  scale_y_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.85))+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red3"), limits=c(0,1))+
  theme_classic()+
  theme(axis.line=element_line(size=0.5))+
  theme(legend.position = "none")
dev.off()

#linear mixed-effect model
all(meta_combine$ID == rownames(df_sample_signature_cluster))
df_sample_signature_cluster$patient_ID = meta_combine$patient

# Fit mixed-effects model
model <- lmer(
  sample_avg_top1000 ~ setting_detail + (1 | patient_ID),
  data = df_sample_signature_cluster
)

# Overall test for setting_detail
anova(model)



##########check development of cluster signature in several samples from recurrence patients
###add sample avg of top1000 to meta_combine
all(names(sample_avg_top1000) == meta_combine$ID)
meta_combine$cluster_avg = sample_avg_top1000
#subset by only recurrent cases (primaries and recurrences)
targets_rec = meta_combine[meta_combine$Cohort=="longitudinal",]

cairo_pdf(filename = "Lineplot_cluster_signature_development_patient.pdf", width = 5, height = 6)
ggplot(targets_rec, aes(x=setting_detail, y = cluster_avg, group = patient))+
  geom_line(alpha=0.8)+
  geom_point(aes(color=cluster_avg),size=6)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red3"), limits=c(0,1))+
  theme_classic()+
  theme(axis.line=element_line(size=0.5))+
  theme(legend.position = "none")
dev.off()

save.image()



###########check differential methylation, only include probes that are setting (meninges, Pri_no, Pri_yes, Recurrence) or MNG specific, then do PCA
###check meta and beta data
all(meta_combine$ID == colnames(beta.combine))

#make summarized experiment
se.combine <- SummarizedExperiment(as.matrix(beta.combine), colData = meta_combine)
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

#differential methylation for setting contrast
smry_combine = DML(se.combine, ~setting)
smry_combine
test_result_combine = summaryExtractTest(smry_combine)
colnames(test_result_combine)


#differential methylation for only contrast primary_no_rec, primary_rec, and recurrences (MNG probes)
#generate probes specific for only MNG, without meninges
#make subsets of bates and meta first
beta.MNG = beta.combine[,-c(1:20)]
meta_MNG = meta_combine[-c(1:20),]
#make summarized experiment
se.MNG <- SummarizedExperiment(as.matrix(beta.MNG), colData = meta_MNG)
str(se.MNG)
colData(se.MNG)$setting = relevel(factor(colData(se.MNG)$setting),"Primary_no_rec")
se_MNG_ok = (checkLevels(assay(se.MNG), colData(se.MNG)$setting))
sum(se_MNG_ok)
#no NAs

#differential methylation for MNG contrast
smry_MNG = DML(se.MNG, ~setting)
test_result_MNG = summaryExtractTest(smry_MNG)
colnames(test_result_MNG)

save.image(file = "my_work_space.Rdata")

#check setting specific probes (113,047)
test_result_setting = test_result_combine %>%
  mutate(setting_specific =
           ifelse(FPval_setting < 0.01 & Eff_setting > 0.1, TRUE, FALSE)) 
setting_probes_combine = test_result_setting %>% dplyr::filter(setting_specific=="TRUE")
dim(setting_probes_combine)


#check MNG specific probes (19,142)
test_result_MNG = test_result_MNG %>%
  mutate(setting_specific_MNG =
           ifelse(FPval_setting < 0.01 & Eff_setting > 0.1, TRUE, FALSE))
setting_MNG_probes = test_result_MNG %>% dplyr::filter(setting_specific_MNG=="TRUE")
dim(setting_MNG_probes)


#####subset
#make beta values subset of setting and setting_detail specific probes, and do PCA
bVals_combine_setting_probes = beta.combine %>% filter(rownames(beta.combine) %in% (setting_probes_combine$Probe_ID))
dim(bVals_combine_setting_probes)

bVals_combine_setting_probes_MNG = beta.combine %>% filter(rownames(beta.combine) %in% (setting_MNG_probes$Probe_ID))
dim(bVals_combine_setting_probes_MNG)

#make PCA
bVals_combine_setting_PCA = prcomp(t(bVals_combine_setting_probes), scale. = T)
bVals_combine_setting_MNG_PCA = prcomp(t(bVals_combine_setting_probes_MNG), scale. = T)

df_PCA_setting_all <- cbind(meta_combine, bVals_combine_setting_PCA$x[,1:3])
df_PCA_setting_MNG <- cbind(meta_combine, bVals_combine_setting_MNG_PCA$x[,1:3])


#scree plots

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
#optimal number of clusters k for all setting probes
sub <- bVals_combine_setting_probes[sample(nrow(bVals_combine_setting_probes), 2000), ]

tiff("Elbow_clusters_combine_probes.tiff", res=300,width=1548,height=1548)
fviz_nbclust(sub, kmeans,  method = "wss")+
  geom_vline(xintercept = 2, linetype=2)
dev.off()
#optimal k=2


#perform kmeans clustering for all setting probes
k=2
km_res = kmeans(t(bVals_combine_setting_probes), centers = k, nstart = 25)
df_PCA_setting_all$kmeans = as.factor(km_res$cluster)
str(df_PCA_setting_all)

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
  scale_linetype_manual(values = c("solid", "dashed"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.8))+
  theme(legend.position = "none")
dev.off()



#optimal number of clusters k for MNG setting probes
sub <- bVals_combine_setting_probes_MNG[sample(nrow(bVals_combine_setting_probes_MNG), 2000), ]

tiff("Elbow_clusters_combine_MNG_probes.tiff", res=300,width=1548,height=1548)
fviz_nbclust(sub, kmeans,  method = "wss")+
  geom_vline(xintercept = 2, linetype=2)
dev.off()

#perform kmeans clustering for MNG probes
k=2
km_res = kmeans(t(bVals_combine_setting_probes_MNG), centers = k, nstart = 25)
df_PCA_setting_MNG$kmeans = as.factor(km_res$cluster)
str(df_PCA_setting_MNG)

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
#differential methylation for MNG settings vs normal meninges including pvalue adjustment

se_kmeans <- SummarizedExperiment(as.matrix(beta.combine), colData = meta_combine)

smry_meninges_contrast = DML(se_kmeans, ~setting)
test_result_meninges_contrast = summaryExtractTest(smry_meninges_contrast)
colnames(test_result_meninges_contrast)

padj_Pri_no_rec = p.adjust(test_result_meninges_contrast$Pval_settingPrimary_no_rec, method = "BH", n = length(test_result_meninges_contrast$Pval_settingPrimary_no_rec))
padj_Pri_rec = p.adjust(test_result_meninges_contrast$Pval_settingPrimary_rec, method = "BH", n = length(test_result_meninges_contrast$Pval_settingPrimary_rec))
padj_Rec = p.adjust(test_result_meninges_contrast$Pval_settingRecurrence, method = "BH", n = length(test_result_meninges_contrast$Pval_settingRecurrence))
test_result_meninges_contrast = cbind(test_result_meninges_contrast, padj_Pri_no_rec, padj_Pri_rec, padj_Rec)

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
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% clust_sig), 
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
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% clust_sig), 
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
  geom_point(data= test_result_meninges_contrast %>% filter(Probe_ID %in% clust_sig), 
             aes(x=Est_settingRecurrence, y=-log10(padj_Rec)), color="red3", alpha=0.7, shape=16)+
  scale_radius(range = c(0.05,5))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.6,0.25), limits=c(-0.5,0.7))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth = 0.8))+
  theme(legend.position = "none")
dev.off()



###perfom slingshot trajectory
##subset all data by MNG probes and make matrix
bVals_combine_MNG_probes = beta.combine %>% filter(rownames(beta.combine) %in% (setting_MNG_probes$Probe_ID))
bVals_combine_sling = as.matrix(bVals_combine_MNG_probes)

#make singlecellexperiment
sce_combine <- SingleCellExperiment(assays = List(counts = bVals_combine_sling))
assays(sce_combine)$norm <- bVals_combine_sling

#perform PCA and add values to singlecellexperiment
pca_combine_sling <- prcomp(t(assays(sce_combine)$norm), scale. = TRUE)
rd1_combine <- pca_combine_sling$x[,1:2]
rd1_combine=as.data.frame(rd1_combine)
rd1_combine$PC2 = rd1_combine$PC2 * -1
rd1_combine=as.matrix(rd1_combine)

plot(rd1_combine, col = rgb(0,0,0,.5), pch=16, asp = 1)
reducedDims(sce_combine) <- SimpleList(PCA = rd1_combine)

#add information on clusters
meta_combine$cluster_development = c(rep(1,20), rep(2,13), rep(3,42))
all(meta_combine$ID == sce_combine@colData@rownames)

cl1_combine = meta_combine$cluster_development
names(cl1_combine) = meta_combine$ID
cl1_combine
all(colData(sce_combine)$rownames == names(cl1_combine))
colData(sce_combine)$GMM <- cl1_combine

#check clusters on PCA
library(RColorBrewer)
plot(rd1_combine, col = brewer.pal(9,"Set1")[cl1_combine], pch=16, asp = 1)

#perform sling
sce_combine <- slingshot(sce_combine, clusterLabels = 'GMM', reducedDim = 'PCA')
summary(sce_combine$slingPseudotime_1)

#show trajectory in plot
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

###use ggplot to visualize the trajectory inference
#make dataframe
df_combine = data.frame(rd1_combine, "cl1_combine" = as.character(cl1_combine))
all(rownames(df_combine) == meta_combine$ID)
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
all(rownames(df_combine)==meta_combine$ID)
df_combine= df_combine[match(meta_combine$ID,rownames(df_combine)),]
all(rownames(df_combine)==meta_combine$ID)

df_combine$avg_signature = meta_combine$cluster_avg
df_combine$MCconsensus = meta_combine$MCconsensus
df_combine$setting_detail = meta_combine$setting_detail
df_combine$grade = meta_combine$grade
df_combine$grade = factor(df_combine$grade, levels = c("0","1","2","3"))
rownames(df_combine) = df_combine$X
str(df_combine)

#smooth trajectory curve
smoothed_curve = smooth.spline(x=df_combine$curve_x, y=df_combine$curve_y, spar = 0.95)
df_smoothed_curve = as.data.frame(cbind(smoothed_curve$x, smoothed_curve$y))

###########make trajectory plots
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
pseudotime_combine$avg_signature = df_combine$avg_signature
pseudotime_combine$METH_cluster = df_combine$cluster
pseudotime_combine$METH_cluster = factor(pseudotime_combine$METH_cluster, levels = c("0","1","2"))
pseudotime_combine$setting_detail = df_combine$setting_detail

#METH acitvity
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

#Pearson's product-moment correlation

#data:  pseudotime_combine$pseudotime and pseudotime_combine$avg_signature
#t = 16.728, df = 73, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
 #0.8317020 0.9296221
#sample estimates:
 #    cor 
#0.890558 


#MNG settings
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


#hide label meninges before plotting METH cluster
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












#####Patient-aware sensitivity analysis
#####Collapse samples within patient (pseudo-bulk / centroid approach)
#collapse on PCA positions
pcs_df <- as.data.frame(rd1_combine)
all(rownames(pcs_df) == meta_combine$ID)
pcs_df$patient_id <- meta_combine$patient

#collapse
pcs_centroid <- aggregate(
  pcs_df[, 1:3],
  by = list(patient_id = pcs_df$patient_id),
  FUN = mean
)
rownames(pcs_centroid) = pcs_centroid$patient_id

#make PCAs a matrix for getLineages
pcs_centroid_sling = as.matrix(pcs_centroid[,c(1:2)])

#make cluster named numeric
cl1_combine_centroid = pcs_centroid$centroid_ID
names(cl1_combine_centroid) = pcs_centroid$patient_id.1

#define linages by start and end point
lin1 <- getLineages(pcs_centroid_sling, cl1_combine_centroid, start.clus="1", end.clus="3")
lin1

#get ciurve
crv1 <- getCurves(lin1,extend="pc1")
crv1

slingPseudotime(crv1)


#ggplot for plot
#make dataframe
df_combine = data.frame(pcs_centroid_sling, "cl1_combine_centroid" = as.character(cl1_combine_centroid))

#make curve data a dataframe
crv1_plot = SlingshotDataSet(crv1)
crv1_data = as.data.frame(crv1_plot@curves$Lineage1$s)

df_combine$curve_x = crv1_data$PC1
df_combine$curve_y = crv1_data$PC2

#add methylation signature and setting information
all(rownames(df_combine)==rownames(pcs_centroid))

df_combine$avg_signature = pcs_centroid$avg_cluster
df_combine$centroid_ID = pcs_centroid$centroid_ID

#smooth trajectory curve
smoothed_curve = smooth.spline(x=df_combine$curve_x, y=df_combine$curve_y, spar = 0.95)
df_smoothed_curve = as.data.frame(cbind(smoothed_curve$x, smoothed_curve$y))

###########make trajectory plots
#for METH activity
cairo_pdf(filename = "PCA_Trajectory_centroid_combine_METH_activity_with_curve.pdf", width = 8.5, height = 6)
ggplot(df_combine, aes(x = PC1, y = PC2, color=avg_signature)) +
  geom_point(show.legend = F, size=6)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red"))+
  geom_line(data=df_smoothed_curve, aes(V1, V2),color="black", size=1.2)+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()

ggplot(df_combine, aes(x = PC1, y = PC2, color=avg_signature)) +
  geom_point(show.legend = F, size=6)


#for patient-level centroid conditions
cairo_pdf(filename = "PCA_Trajectory_centroid_combine_conditions_with_curve.pdf", width = 8.5, height = 6)
ggplot(df_combine, aes(x = PC1, y = PC2, color=centroid_ID_chr)) +
  geom_point(show.legend = F, size=6)+
  scale_color_manual(values = c(meninges="dodgerblue3",
                                    benign_cases="#1CFDB2",
                                    recurrent_cases="#FF9933"))+
  geom_line(data=df_smoothed_curve, aes(V1, V2),color="black", size=1.2)+
  theme_classic()+
  theme(axis.line = element_line(linewidth=0.9))+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()


####make pseudotime plots versus setting and METH activity
#get data for psdeudotime/setting
pseudotime_combine = as.data.frame(slingPseudotime(crv1))
pseudotime_combine$pseudotime = pseudotime_combine$Lineage1
pseudotime_combine$avg_signature = df_combine$avg_signature
pseudotime_combine$centroid_ID = df_combine$centroid_ID_chr
pseudotime_combine$centroid_ID <- factor(
  pseudotime_combine$centroid_ID,
  levels = c("meninges", "benign_cases", "recurrent_cases")
)


#METH acitvity
cairo_pdf(filename = "METH_activity_along_pseudotime_centroids.pdf", width = 7, height = 8)
ggplot(pseudotime_combine, aes(x = pseudotime, y = avg_signature, color=avg_signature)) +
  geom_point(show.legend = F, size=8)+
  scale_color_gradientn(colours = c("navy","dodgerblue3","indianred1","red"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()



cor.test(pseudotime_combine$pseudotime, pseudotime_combine$avg_signature)


#centroid IDs
cairo_pdf(filename = "Centroid_ID_along_pseudotime.pdf", width = 7, height = 8)
ggplot(pseudotime_combine, aes(x = pseudotime, y = centroid_ID, color=centroid_ID)) +
  geom_point(show.legend = F, size=8)+
  scale_color_manual(values = c(meninges="dodgerblue3",
                                benign_cases="#1CFDB2",
                                recurrent_cases="#FF9933"))+
  theme_classic()+
  theme(axis.line = element_line(linewidth=1)
  )+
  theme(axis.text.x=element_blank(),axis.title.x = element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()
  )
dev.off()


save.image(file="my_work_space.RData")
load("my_work_space.RData")


