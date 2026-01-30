library(ggplot2)
library(factoextra)
library(ggpubr)
library(VennDiagram)
library(grid)
library(NbClust)
library(pheatmap)
library(ConsensusClusterPlus)
library(tidyverse)
library(cluster)
library(mclust)
library(parallel)
library(dendextend)
library(mclust)
library(dplyr)
library(networkD3)
library(htmlwidgets)
library(webshot2)
library(sesame)
library(ggstatsplot)


###make sankey for disovery cohort

Disc_Sankey = read.csv(file="Disc_For_Sankey.csv")

nodes <- data.frame(
  name = unique(c(
    Disc_Sankey$WHO_grade,
    Disc_Sankey$setting,
    Disc_Sankey$Toronto,
    Disc_Sankey$DKFZ,
    Disc_Sankey$risk
  )),
  stringsAsFactors = FALSE
)
nodes$label <- ""

links_1 <- Disc_Sankey %>%
  group_by(WHO_grade, setting) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(WHO_grade, nodes$name) - 1,
    target = match(setting, nodes$name) - 1
  ) %>%
  select(source, target, value)

links_2 <- Disc_Sankey %>%
  group_by(setting, Toronto) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(setting, nodes$name) - 1,
    target = match(Toronto, nodes$name) - 1
  ) %>%
  select(source, target, value)

links_3 <- Disc_Sankey %>%
  group_by(Toronto, DKFZ) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(Toronto, nodes$name) - 1,
    target = match(DKFZ, nodes$name) - 1
  ) %>%
  select(source, target, value)

links_4 <- Disc_Sankey %>%
  group_by(DKFZ, risk) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(DKFZ, nodes$name) - 1,
    target = match(risk, nodes$name) - 1
  ) %>%
  select(source, target, value)

links_4 = links_4[-7,]

links <- bind_rows(links_1, links_2, links_3, links_4)

nodes$group <- NA_character_

# WHO grades
nodes$group[nodes$name == "1"] <- "Grade1"
nodes$group[nodes$name == "2"] <- "Grade2"
nodes$group[nodes$name == "3"] <- "Grade3"

# setting
nodes$group[nodes$name == "Primary"]     <- "Primary"
nodes$group[nodes$name == "Recurrence"]   <- "Recurrence"

# Toronto
nodes$group[nodes$name == "Merlinintact"]     <- "Merlin"
nodes$group[nodes$name == "Immuneenriched"]   <- "Immune"
nodes$group[nodes$name == "hypermetabolic"]   <- "Hyper"
nodes$group[nodes$name == "proliferative"]    <- "Prolif"

# DKFZ
nodes$group[nodes$name == "Benign"]  <- "Benign"
nodes$group[nodes$name == "Intermediate"] <- "Intermediate"
nodes$group[nodes$name == "SMARCE1altered"] <- "SMARCE1"
nodes$group[nodes$name == "Malignant"] <- "Malignant"

# risk
nodes$group[nodes$name == "low"]  <- "low"
nodes$group[nodes$name == "intermediate"] <- "intermediate"
nodes$group[nodes$name == "high"] <- "high"

stopifnot(!any(is.na(nodes$group)))



my_colour_scale <- '
d3.scaleOrdinal()
  .domain([
    "Grade1","Grade2","Grade3",
    "Primary", "Recurrence",
    "Merlin","Immune","Hyper","Prolif",
    "Benign", "Intermediate", "SMARCE1", "Malignant",
    "low","intermediate", "high"
  ])
  .range([
    "#9EBCDA", "#8C6BB1", "#810F7C",
    "#1B9E77", "#D95F02",
    "#436eee", "#cd0000", "#228b22", "#ee7600",
    "#386CB0", "#F0027F", "#666666", "#BF5B17",
    "#1c86ee", "#9b30ff", "#ff3030"             
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





#get beta values and meta
betas <- readRDS("Betas_discovery_preprocessed_corrected.rds")
meta = read.csv(file="meta_discovery_current.csv")
rownames(meta)=meta$ID


#subset to top variance
betas_sub = as.data.frame(betas)
str(betas_sub)
betas_sub$var = apply(betas_sub,1,var)
betas_sub <- betas_sub[order(betas_sub$var, decreasing = TRUE),]
betas_sub = betas_sub[,-232]
betas_10k = betas_sub[1:10000,]
betas_50k = betas_sub[1:50000,]
betas_200k = betas_sub[1:200000,]


#perform PCA
betas_10k_pca = prcomp(t(betas_10k), scale. = T)
betas_50k_pca = prcomp(t(betas_50k), scale. = T)
betas_200k_pca = prcomp(t(betas_200k), scale. = T)


#make scree plot

fviz_screeplot(betas_10k_pca,addlabels = T)
tiff("Scree_PCs_top10k.tiff", res=300,width=2100,height=1400)
fviz_screeplot(betas_10k_pca)
dev.off()

fviz_screeplot(betas_50k_pca,addlabels = T)
tiff("Scree_PCs_top50k.tiff", res=300,width=2100,height=1400)
fviz_screeplot(betas_50k_pca)
dev.off()

fviz_screeplot(betas_200k_pca,addlabels = T)
tiff("Scree_PCs_top200k.tiff", res=300,width=2100,height=1400)
fviz_screeplot(betas_200k_pca)
dev.off()


#make PCA plot with margins
df_10 <- cbind(meta, betas_10k_pca$x[,1:4])
df_50 <- cbind(meta, betas_50k_pca$x[,1:4])
df_200 <- cbind(meta, betas_200k_pca$x[,1:4])

df_10$grading2021_new = factor(df_10$grading2021_new)
df_50$grading2021_new = factor(df_50$grading2021_new)
df_200$grading2021_new = factor(df_200$grading2021_new)

df_10$mol_group = factor(df_10$mol_group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
df_50$mol_group = factor(df_50$mol_group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
df_200$mol_group = factor(df_200$mol_group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

df_10$risk_score_new = factor(df_10$risk_score_new, levels = c("low", "intermediate", "high"))
df_50$risk_score_new = factor(df_50$risk_score_new, levels = c("low", "intermediate", "high"))
df_200$risk_score_new = factor(df_200$risk_score_new, levels = c("low", "intermediate", "high"))


#for CNS WHO grade, several subsets, PC1-4

tiff("PCA_10k_margin_grade_PC1-2_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC2",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_10k_margin_grade_PC1-3_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC3",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_10k_margin_grade_PC1-4_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC4",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


tiff("PCA_50k_margin_grade_PC1-2_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC2",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_grade_PC1-3_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC3",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_grade_PC1-4_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC4",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


tiff("PCA_200k_margin_grade_PC1-2_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC2",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_200k_margin_grade_PC1-3_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC3",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_200k_margin_grade_PC1-4_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC4",
  color = "grading2021_new", size = 5, alpha = 0.8,
  palette = c("#9EBCDA","#8C6BB1","#810F7C"),
  margin.params = list(fill = "grading2021", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

####calculate weigthed pairwise mean overlaps of 3 grades along all four PCs for three probe selections
###this will serve as a single metric for the capability of both probe selection and PCs to separate the clinically relevant molecular groups
###PCs, 10k data (run iteratively for the first 4 PCs)
rownames(meta) = meta$ID
all(rownames(betas_10k_pca$x) == rownames(meta))
pc_scores <- betas_10k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$grading2021_new
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["1"]]
g2 = pc1_groups[["2"]]
g3 = pc1_groups[["3"]]


#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl



###PCs, 50k data (run iteratively for the first 4 PCs)
pc_scores <- betas_50k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$grading2021_new
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["1"]]
g2 = pc1_groups[["2"]]
g3 = pc1_groups[["3"]]

#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl


###PCs, 200k data (run iteratively for the first 4 PCs)
pc_scores <- betas_200k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$grading2021_new
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["1"]]
g2 = pc1_groups[["2"]]
g3 = pc1_groups[["3"]]

#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl


#for molecular group, several subsets, PC1-4

tiff("PCA_10k_margin_group_PC1-2.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC2",
  color = "mol_group", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "mol_group", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_10k_margin_group_PC1-3.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC3",
  color = "mol_group", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "mol_group", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_10k_margin_group_PC1-4.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC4",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


tiff("PCA_50k_margin_group_PC1-2.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC2",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_group_PC1-3.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC3",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_group_PC1-4.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC4",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


tiff("PCA_200k_margin_group_PC1-2.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC2",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_200k_margin_group_PC1-3.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC3",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_200k_margin_group_PC1-4.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC4",
  color = "MCconsensus", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MCconsensus", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


####calculate weigthed pairwise mean overlaps of the 4 molecular groups along all four PCs for three probe selections
###this will serve as a single metric for the capability of both probe selection and PCs to separate the clinically relevant molecular groups

###PCs, 10k data (run iteratively for the first 4 PCs)
all(rownames(betas_10k_pca$x) == rownames(meta))
pc_scores <- betas_10k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$MCconsensus
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["Merlinintact"]]
g2 = pc1_groups[["Immuneenriched"]]
g3 = pc1_groups[["hypermetabolic"]]
g4 = pc1_groups[["proliferative"]]

#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3, G4 = g4)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl



###PCs, 50k data (run iteratively for the first 4 PCs)
pc_scores <- betas_50k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$MCconsensus
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["Merlinintact"]]
g2 = pc1_groups[["Immuneenriched"]]
g3 = pc1_groups[["hypermetabolic"]]
g4 = pc1_groups[["proliferative"]]

#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3, G4 = g4)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl


###PCs, 200k data (run iteratively for the first 4 PCs)
pc_scores <- betas_200k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$MCconsensus
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["Merlinintact"]]
g2 = pc1_groups[["Immuneenriched"]]
g3 = pc1_groups[["hypermetabolic"]]
g4 = pc1_groups[["proliferative"]]

#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3, G4 = g4)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl




#for integrated risk, several subsets, PC1-4

tiff("PCA_10k_margin_risk_PC1-2_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC2",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_10k_margin_risk_PC1-3_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC3",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_10k_margin_risk_PC1-4_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_10, x = "PC1", y = "PC4",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


tiff("PCA_50k_margin_risk_PC1-2_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC2",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_risk_PC1-3_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC3",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_50k_margin_risk_PC1-4_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_50, x = "PC1", y = "PC4",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


tiff("PCA_200k_margin_risk_PC1-2_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC2",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_200k_margin_risk_PC1-3_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC3",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_200k_margin_risk_PC1-4_new.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_200, x = "PC1", y = "PC4",
  color = "risk_score_new", size = 5, alpha = 0.8,
  palette = c("dodgerblue2","purple1","firebrick1"),
  margin.params = list(fill = "risk_score", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


###PCs, 10k data (run iteratively for the first 4 PCs)
pc_scores <- betas_10k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$risk_score_new
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["low"]]
g2 = pc1_groups[["intermediate"]]
g3 = pc1_groups[["high"]]


#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl


###PCs, 50k data (run iteratively for the first 4 PCs)
pc_scores <- betas_50k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$risk_score_new
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["low"]]
g2 = pc1_groups[["intermediate"]]
g3 = pc1_groups[["high"]]


#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl


###PCs, 200k data (run iteratively for the first 4 PCs)
pc_scores <- betas_200k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta$risk_score_new
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["low"]]
g2 = pc1_groups[["intermediate"]]
g3 = pc1_groups[["high"]]


#get weigthed pairwise mean overlap
groups <- list(G1 = g1, G2 = g2, G3 = g3)
ns <- sapply(groups, length)
overlap_coef <- function(x, y, n = 1024) {
  r <- range(c(x, y))
  d1 <- density(x, from = r[1], to = r[2], n = n)
  d2 <- density(y, from = r[1], to = r[2], n = n)
  
  dx <- d1$x[2] - d1$x[1]
  sum(pmin(d1$y, d2$y)) * dx
}
pairs <- combn(names(groups), 2, simplify = FALSE)

pair_ovl <- sapply(pairs, function(p) {
  overlap_coef(groups[[p[1]]], groups[[p[2]]])
})
pair_weights <- sapply(pairs, function(p) {
  ns[p[1]] * ns[p[2]]
})
weighted_mean_ovl <- sum(pair_ovl * pair_weights) / sum(pair_weights)
weighted_mean_ovl


save.image()

########using top 10.000 most variable probes provides best separation of clinically relevant categories (molecular groups as well as integrated risk score) along PC1
####for probe selection and clustering approaches, 10.000 probe-set will be used
###########get probes that contribute most to PC1, as PC1 best captures difference between high grade and low grade

###first, get top2.000 probes contributing to PC1 for each probe selection (10k, 50k, and 200k), and check overlap of those using a venn diagramm

var_10k = get_pca_var(betas_10k_pca)
head(var_10k$contrib,10)
var_contrib_10k = as.data.frame(var_10k$contrib)
var_contrib_sortPC1 = var_contrib_10k[order(var_contrib_10k$Dim.1, decreasing = TRUE),]
top_PC1_10k = data.frame(PC1=rownames(var_contrib_sortPC1))
top2000_probes_PC1_10k = c(top_PC1_10k$PC1[1:2000])
top10000_probes_PC1_10k = c(top_PC1_10k$PC1[1:10000])

var_50k = get_pca_var(betas_50k_pca)
head(var_50k$contrib,10)
var_contrib_50k = as.data.frame(var_50k$contrib)
var_contrib_sortPC1 = var_contrib_50k[order(var_contrib_50k$Dim.1, decreasing = TRUE),]
top_PC1_50k = data.frame(PC1=rownames(var_contrib_sortPC1))
top2000_probes_PC1_50k = c(top_PC1_50k$PC1[1:2000])
top10000_probes_PC1_50k = c(top_PC1_50k$PC1[1:10000])

var_200k = get_pca_var(betas_200k_pca)
head(var_200k$contrib,10)
var_contrib_200k = as.data.frame(var_200k$contrib)
var_contrib_sortPC1 = var_contrib_200k[order(var_contrib_200k$Dim.1, decreasing = TRUE),]
top_PC1_200k = data.frame(PC1=rownames(var_contrib_sortPC1))
top2000_probes_PC1_200k = c(top_PC1_200k$PC1[1:2000])
top10000_probes_PC1_200k = c(top_PC1_200k$PC1[1:10000])

#make Venn diagramm for overlap of top probes contributing to PC1
venn_list <- list(
  A = top10000_probes_PC1_10k,
  B = top10000_probes_PC1_50k,
  C = top10000_probes_PC1_200k
)




venn.plot <- venn.diagram(
  x = list(
    A = top10000_probes_PC1_10k,
    B = top10000_probes_PC1_50k,
    C = top10000_probes_PC1_200k
  ),
  filename = NULL,
  fill = c("red", "magenta3", "dodgerblue3"),   # <-- group colors
  alpha = 0.5,                        # <-- transparency → overlaps blend
  cex = 1,                            # <-- no numbers
  cat.cex = 2                         # <-- labels visible
)

grid.draw(venn.plot)


###from here on, work with 10k probe data selection
#get top contributing probes 
var_10k = get_pca_var(betas_10k_pca)
head(var_10k$contrib,10)
var_contrib_10k = as.data.frame(var_10k$contrib)
var_contrib_sortPC1 = var_contrib_10k[order(var_contrib_10k$Dim.1, decreasing = TRUE),]
top_PC1_10k = data.frame(PC1=rownames(var_contrib_sortPC1))

#####get selected subsets of probes associated with PC1 in the 10k dataset

top1000_probes_PC1_10k = c(top_PC1_10k$PC1[1:1000])
top1500_probes_PC1_10k = c(top_PC1_10k$PC1[1:1500])
top2500_probes_PC1_10k = c(top_PC1_10k$PC1[1:2500])
top4000_probes_PC1_10k = c(top_PC1_10k$PC1[1:4000])
top6000_probes_PC1_10k = c(top_PC1_10k$PC1[1:6000])


#subset beta values by top features for PC1
df_beta = as.data.frame(betas)

beta_top1000 = df_beta %>% dplyr::filter(rownames(df_beta) %in% top1000_probes_PC1_10k)
beta_top1500 = df_beta %>% dplyr::filter(rownames(df_beta) %in% top1500_probes_PC1_10k)
beta_top2500 = df_beta %>% dplyr::filter(rownames(df_beta) %in% top2500_probes_PC1_10k)
beta_top4000 = df_beta %>% dplyr::filter(rownames(df_beta) %in% top4000_probes_PC1_10k)
beta_top6000 = df_beta %>% dplyr::filter(rownames(df_beta) %in% top6000_probes_PC1_10k)



#check number of clusters for subsets of PC1-associated probes
#use combination of Elbow method and Silhouette Analysis
#scale data

beta_top1000_scaled = scale(beta_top1000)
beta_top1500_scaled = scale(beta_top1500)
beta_top2500_scaled = scale(beta_top2500)
beta_top4000_scaled = scale(beta_top4000)
beta_top6000_scaled = scale(beta_top6000)


# Elbow method
fviz_nbclust(beta_top500_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2, color = "steelblue")+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(beta_top500_scaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")



# Elbow method
fviz_nbclust(beta_top1000_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2, color = "steelblue")+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(beta_top1000_scaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")


# Elbow method
fviz_nbclust(beta_top2000_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2, color = "steelblue")+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(beta_top2000_scaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")


# Elbow method
fviz_nbclust(beta_top3000_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2, color = "steelblue")+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(beta_top3000_scaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")


# Elbow method
fviz_nbclust(beta_top4000_scaled, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2, color = "steelblue")+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(beta_top4000_scaled, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")


#######show several silhouette plots in one plot

probe_list <- list(
  Top1000 = beta_top1000_scaled,
  Top1500 = beta_top1500_scaled,
  Top2500 = beta_top2500_scaled,
  Top4000 = beta_top4000_scaled,
  Top6000 = beta_top6000_scaled
)

compute_sil_width <- function(data, k.max = 10) {
  sapply(2:k.max, function(k) {
    km <- kmeans(data, centers = k, nstart = 25)
    ss <- silhouette(km$cluster, dist(data))
    mean(ss[, "sil_width"])
  })
}

sil_df <- lapply(names(probe_list), function(name) {
  widths <- compute_sil_width(probe_list[[name]], k.max = 10)
  data.frame(
    k = 2:10,
    sil_width = widths,
    probe_set = name
  ) %>%
    add_row(k = 1, sil_width = 0, probe_set = name) %>%  # Add k=1 with 0
    arrange(k)
}) %>% bind_rows()

#make the plot
ggplot(sil_df, aes(x = k, y = sil_width, color = probe_set)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:10) +
  scale_color_brewer(palette = "Dark2") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "black", size = 1) +
  labs(
    title = "Average Silhouette Width by Probe Set",
    x = "Number of Clusters (k)",
    y = "Average Silhouette Width",
    color = "Probe Set"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),           # remove all grid lines
    axis.line = element_line(color = "black"),  # show axes
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  )


#same for Elbow plot
compute_wss <- function(data, k.max = 10) {
  sapply(1:k.max, function(k) {
    km <- kmeans(data, centers = k, nstart = 10)
    km$tot.withinss
  })
}

wss_df <- lapply(names(probe_list), function(name) {
  wss <- compute_wss(probe_list[[name]], k.max = 10)
  data.frame(
    k = 1:10,
    wss = wss,
    probe_set = name
  )
}) %>% bind_rows()

ggplot(wss_df, aes(x = k, y = wss, color = probe_set)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_y_continuous(trans = "log1p")+
  scale_x_continuous(breaks = 1:10) +
  scale_color_brewer(palette = "Dark2") +
  geom_vline(xintercept = 2, linetype = "dotted", color = "black", size = 1) +
  labs(
    title = "Elbow Plot: WSS by Probe Set",
    x = "Number of Clusters (k)",
    y = "Total Within-Cluster Sum of Squares (WSS)",
    color = "Probe Set"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),           # remove grid lines
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black")
  )


###consensus clustering for subsets of probes to find cluster stability

set.seed(123)
out_dir <- "consensus_out"
dir.create(out_dir, showWarnings = FALSE)

run_consensus_for_n <- function(n_probes, meth_mat, probes_ranked, out_dir = "consensus_out") {
  
  # 1. Select top probes
  selected_probes <- probes_ranked[1:n_probes]
  data_subset <- meth_mat[selected_probes, , drop = FALSE]
  
  # 2. Create output folder (e.g., "consensus_out/cc_100_probes")
  title_name <- paste0("cc_", n_probes, "_probes")
  out_path <- file.path(out_dir, title_name)
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  
  # 3. Temporarily move into output directory
  old_wd <- getwd()
  setwd(out_path)
  
  # 4. Run ConsensusClusterPlus
  cc <- ConsensusClusterPlus(
    as.matrix(data_subset),
    maxK = 3,                   
    reps = 100,
    pItem = 0.8,
    pFeature = 0.8,
    clusterAlg = "hc",
    distance = "euclidean",
    innerLinkage = "ward.D2",
    finalLinkage = "ward.D2",
    seed = 123,
    title = title_name,        
    writeTable = TRUE,             
    plot = "pdf"                  
  )
  
  # 5. Return to original working directory
  setwd(old_wd)
  
  # 6. Extract K=2 results
  cc2 <- cc[[2]]
  
  # 7. Return tidy tibble
  tibble(
    sample    = names(cc2$consensusClass),
    cluster   = cc2$consensusClass,
    stability = cc2$ml,
    n_probes  = n_probes
  )
}

meth_mat <- betas_10k
pca <- prcomp(t(meth_mat), center = TRUE, scale. = TRUE)
pc1_loadings <- pca$rotation[, 1]
var_pc1 <- summary(pca)$importance[2, 1]

probes_ranked <- abs(pc1_loadings) %>% sort(decreasing = TRUE) %>% names()
length(probes_ranked)

probe_set_sizes <- c(1000, 1500, 2500, 4000, 6000)
probe_set_sizes <- probe_set_sizes[probe_set_sizes <= length(probes_ranked)]
probe_set_sizes

results_list <- list()
for (n in probe_set_sizes) {
  message("Running consensus for ", n, " probes ...")
  res <- run_consensus_for_n(n, meth_mat, probes_ranked, out_dir = out_dir)
  results_list[[as.character(n)]] <- res
}


results_all <- bind_rows(results_list)


stability_summary <- results_all |>
  group_by(n_probes) |>
  summarise(
    mean_stability = mean(stability),
    median_stability = median(stability),
    sd_stability = sd(stability)
  )

stability_summary

tiff("Cluster_stability_summary.tiff", res=300,width=2560,height=2048)
ggplot(
  data = stability_summary,
  mapping = aes(x = n_probes, y = median_stability, group = 1)) +
  geom_line(linewidth = 1, color = "steelblue") +
  geom_point(size = 4, color = "darkred") +
  scale_x_continuous(breaks = stability_summary$n_probes) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
                             axis.ticks = element_line(color = "black")
  ) +
  labs(
    title = "Cluster stability vs number of PC1-associated probes",
    x = "Number of probes",
    y = "Median CCP stability"
  )
dev.off()



####check similarities of cluster assignments for distinct probes subset
#extract cluster assignments

cluster_list <- list(
  clust_1000  = results_list$`1000`$cluster,
  clust_1500  = results_list$`1500`$cluster,
  clust_2500 = results_list$`2500`$cluster,
  clust_4000 = results_list$`4000`$cluster,
  clust_6000 = results_list$`6000`$cluster
)

sample_order <- names(cluster_list[[1]])
cluster_list <- lapply(cluster_list, function(x) x[sample_order])

#similarity of assignments across subsets
subset_names <- names(cluster_list)
n <- length(cluster_list)

similarity_matrix <- matrix(0, n, n,
                            dimnames=list(subset_names, subset_names))

for (i in 1:n) {
  for (j in 1:n) {
    similarity_matrix[i,j] <- adjustedRandIndex(
      cluster_list[[i]],
      cluster_list[[j]]
    )
  }
}


pheatmap(similarity_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         main = "Similarity of Cluster Assignments Across Probe Subset Sizes",
         fontsize_number = 10)


#unsupervised hierarchical clustering

T26_samples = as.data.frame(meta$ID)
colnames(T26_samples) = "ID"
rownames(T26_samples) = T26_samples$ID
all(rownames(T26_samples) %in% meta$ID)
all(rownames(T26_samples) == meta$ID)
T26_samples$risk = meta$risk_score
T26_samples$risk_new = meta$risk_score_new
T26_samples$mol_group = meta$mol_group
T26_samples$MCsubtype = meta$MCsubtype
T26_samples$grade = meta$grading2021
T26_samples$grade = as.character(T26_samples$grade)
T26_samples$grade_new = meta$grading2021_new
T26_samples$grade_new = as.character(T26_samples$grade_new)
T26_samples$status = meta$status
T26_samples$gender = meta$gender
T26_samples = T26_samples[,-1]
T26_samples = T26_samples[,c("risk_new","mol_group", "MCsubtype", "grade_new", "status", "gender")]




#heatmap

ann_colors = list(gender = c(M="#276419", F="#C51B7D"),
                  grade_new=c("1" = "#9EBCDA", "2" = "#8C6BB1", "3" = "#810F7C" ),
                  MCsubtype =c(Benign = "#386CB0", Intermediate = "#F0027F", 
                               Malignant = "#BF5B17",SMARCE1altered = "#666666"),
                  risk_new =c(low="dodgerblue2", intermediate = "purple1",
                          high = "firebrick1"),
                  status = c(Primary="#1B9E77", Recurrence="#D95F02"),
                  mol_group = c(Immuneenriched = "red3",
                                  Merlinintact = "royalblue2",
                                  hypermetabolic = "forestgreen",
                                  proliferative = "darkorange2"))


tiff("Heat_cluster_PC1_top500.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top500, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top1000_new.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top1000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()

#make top1000 heatmap cut by 3 clusters
tiff("Heat_cluster_PC1_top1000_new_3_clusters.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top1000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 3, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top2000.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top2000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top2500_new.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top2500, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top3000.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top3000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()


tiff("Heat_cluster_PC1_top4000.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top4000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top5000.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top5000, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()

tiff("Heat_cluster_PC1_top7500.tiff", res=300,width=2560,height=2048)
pheatmap(beta_top7500, color = colorRampPalette(c("navy", "white", "red"))(100),
         show_rownames = F,show_colnames = F,
         annotation = T26_samples, annotation_colors = ann_colors,clustering_method = "ward.D2",
         cutree_cols = 2, fontsize = 8)
dev.off()

#get clustering results and add to meta data

cluster_1000 = pheatmap(beta_top1000,clustering_method = "ward.D2")
x = cutree(cluster_1000$tree_col, k=2)
all(names(x) == rownames(meta))
meta = cbind(meta, cluster_1000_new = x)

cluster_2500 = pheatmap(beta_top2500,clustering_method = "ward.D2")
x = cutree(cluster_2500$tree_col, k=2)
all(names(x) == rownames(meta))
meta = cbind(meta, cluster_1000_3clusters = x)

#for 3 clusters
cluster_1000_3clusters = pheatmap(beta_top1000,clustering_method = "ward.D2")
x = cutree(cluster_1000_3clusters$tree_col, k=3)
all(names(x) == rownames(meta))
meta = cbind(meta, cluster_1000_new = x)

meta$cluster_1000_3clusters[meta$cluster_1000_3clusters == 1] <- "METHlow-high"
meta$cluster_1000_3clusters[meta$cluster_1000_3clusters== 2] <- "METHhigh"
meta$cluster_1000_3clusters[meta$cluster_1000_3clusters == 3] <- "METHlow-low"



meta$cluster_2000[meta$cluster_200 == 1] <- "METHlow"
meta$cluster_2000[meta$cluster_200 == 2] <- "METHhigh"

meta$cluster_1000_new[meta$cluster_1000_new == 1] <- "METHlow"
meta$cluster_1000_new[meta$cluster_1000_new == 2] <- "METHhigh"

meta$cluster_2500_new[meta$cluster_2500_new == 1] <- "METHlow"
meta$cluster_2500_new[meta$cluster_2500_new == 2] <- "METHhigh"

write.csv(meta, file="meta_with_clusters.csv")

###cluster assignments for top1000 and top25000 are identical
###different cluster assignments for 12 cases as compared to previous top2000 clustering, 11 of which were METHhigh and are now METHlow



#chisquare plots, only for top1000, top2500 are identical
#setting
#clust_1000
length(which(meta$status == "Primary" & meta$cluster_1000_new == "METHlow"))
length(which(meta$status == "Primary" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$status == "Recurrence" & meta$cluster_1000_new == "METHlow"))
length(which(meta$status == "Recurrence" & meta$cluster_1000_new == "METHhigh"))
dat1000 <- data.frame(
  c(144, 14),
  c(36, 37),
  row.names = c("METHlow", "METHhigh"),
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
df1000$clustering = factor(df1000$clustering, levels = c("METHlow", "METHhigh"))

set.seed(123)
test1000 <- chisq.test(table(df1000))

tiff("Chisquare_clust1000_status.tiff", res=300,width=1000,height=2048)
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


#clust_2500
length(which(meta$status == "Primary" & meta$cluster_2500_new == "METHlow"))
length(which(meta$status == "Primary" & meta$cluster_2500_new == "METHhigh"))
length(which(meta$status == "Recurrence" & meta$cluster_2500_new == "METHlow"))
length(which(meta$status == "Recurrence" & meta$cluster_2500_new == "METHhigh"))
dat2500 <- data.frame(
  c(144, 14),
  c(36, 37),
  row.names = c("METHlow", "METHhigh"),
  stringsAsFactors = FALSE
)
colnames(dat2500) <- c("Primary", "Recurrence")
dat2500


x <- c()
for (row in rownames(dat2500)) {
  for (col in colnames(dat2500)) {
    x <- rbind(x, matrix(rep(c(row, col), dat2500[row, col]), ncol = 2, byrow = TRUE))
  }
}
df2500 <- as.data.frame(x)
colnames(df2500) <- c("clustering", "status")
df2500
df2500$clustering = factor(df2500$clustering, levels = c("METHlow", "METHhigh"))

set.seed(123)
test2500 <- chisq.test(table(df2500))

tiff("Chisquare_clust2500_status.tiff", res=300,width=1000,height=2048)
ggbarstats(df2500, status, clustering,
           legend.position = "none",
           results.subtitle = FALSE,
           subtitle = paste0(test2500$p.value),
           label.args = list(
             alpha = 0,
             fill = NA,
             color = NA
           ))+
  scale_fill_manual(values = c("#D95F02","#1B9E77"))+
  theme(legend.position = "none")
dev.off()


##CNS WHO grade
#clust_1000
length(which(meta$grading2021_new == "1" & meta$cluster_1000_new == "METHlow"))
length(which(meta$grading2021_new == "2" & meta$cluster_1000_new == "METHlow"))
length(which(meta$grading2021_new == "3" & meta$cluster_1000_new == "METHlow"))
length(which(meta$grading2021_new == "1" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$grading2021_new == "2" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$grading2021_new == "3" & meta$cluster_1000_new == "METHhigh"))

dat1000 <- data.frame(
  c(79, 2),
  c(94, 32),
  c(7,17),
  row.names = c("METHlow", "METHhigh"),
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
df1000$clustering = factor(df1000$clustering, levels = c("METHlow", "METHhigh"))

test1000 <- fisher.test(table(df1000))

tiff("Chisquare_clust1000_grade_new.tiff", res=300,width=1000,height=2048)
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


#add chisquare for Heidelberg groups 
#clust1000
length(which(meta$MCsubtype == "Benign" & meta$cluster_1000_new == "METHlow"))
length(which(meta$MCsubtype == "Intermediate" & meta$cluster_1000_new == "METHlow"))
length(which(meta$MCsubtype == "Malignant" & meta$cluster_1000_new == "METHlow"))
length(which(meta$MCsubtype == "Benign" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$MCsubtype == "Intermediate" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$MCsubtype == "Malignant" & meta$cluster_1000_new == "METHhigh"))

dat1000 <- data.frame(
  c(137, 0),
  c(42, 44),
  c(0,7),
  row.names = c("METHlow", "METHhigh"),
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
df1000$clustering = factor(df1000$clustering, levels = c("METHlow", "METHhigh"))

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


#add chisquare for Raleigh groups 
#clust1000
length(which(meta$MCconsensus == "Merlinintact" & meta$cluster_1000_new == "METHlow"))
length(which(meta$MCconsensus == "Immuneenriched" & meta$cluster_1000_new == "METHlow"))
length(which(meta$MCconsensus == "hypermetabolic" & meta$cluster_1000_new == "METHlow"))
length(which(meta$MCconsensus == "proliferative" & meta$cluster_1000_new == "METHlow"))
length(which(meta$MCconsensus == "Merlinintact" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$MCconsensus == "Immuneenriched" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$MCconsensus == "hypermetabolic" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$MCconsensus == "proliferative" & meta$cluster_1000_new == "METHhigh"))

dat1000 <- data.frame(
  c(86, 0),
  c(53, 0),
  c(37,9),
  c(4, 42),
  row.names = c("METHlow", "METHhigh"),
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
df1000$clustering = factor(df1000$clustering, levels = c("METHlow", "METHhigh"))
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


#add chisquare for risk_score
#clust1000
length(which(meta$risk_score_new == "low" & meta$cluster_1000_new == "METHlow"))
length(which(meta$risk_score_new == "intermediate" & meta$cluster_1000_new == "METHlow"))
length(which(meta$risk_score_new == "high" & meta$cluster_1000_new == "METHlow"))
length(which(meta$risk_score_new == "low" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$risk_score_new == "intermediate" & meta$cluster_1000_new == "METHhigh"))
length(which(meta$risk_score_new == "high" & meta$cluster_1000_new == "METHhigh"))

dat1000 <- data.frame(
  c(119, 0),
  c(57, 29),
  c(4,22),
  row.names = c("METHlow", "METHhigh"),
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
df1000$clustering = factor(df1000$clustering, levels = c("METHlow", "METHhigh"))
df1000$risk = factor(df1000$risk, levels =  c("low", "intermediate", "high"))

test1000 <- fisher.test(table(df1000))

tiff("Chisquare_clust1000_risk_new.tiff", res=300,width=1000,height=2048)
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
save.image()


###make box whisker for probe-wise average of top 1000 PC1 probes across METHlow and METHhigh
all(rownames(meta) == colnames(beta_top1000))
meta2 <- meta[colnames(beta_top1000), ]

groups <- meta2$cluster_1000_new

row_means_wide <- sapply(
  split(colnames(beta_top1000), groups),
  function(samps) rowMeans(beta_top1000[, samps, drop = FALSE])
)

row_means <- do.call(
  rbind,
  lapply(
    colnames(row_means_wide),
    function(g) {
      data.frame(
        average = row_means_wide[, g],
        group   = g,
        row.names = NULL
      )
    }
  )
)

str(row_means)
row_means$group = factor(row_means$group, levels = c("METHlow", "METHhigh"))


cairo_pdf(filename = "Top1000_probe-wise_averages.pdf", width = 3, height = 8)
ggplot(row_means, 
       aes(group, average))+
  geom_point(position = position_jitter(), alpha=0.95, aes(color=average))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  scale_color_gradientn(colours = c("navy", "white", "red"), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#Discovery cohort
#visualize average methylation of top 1000 PC1 probes across grades and molecular groups, also split by cluster
#grades
grading <- meta$grading2021_new
row_means_wide <- sapply(
  split(colnames(beta_top1000), grading),
  function(sample_set) {
    rowMeans(beta_top1000[, sample_set, drop = FALSE])
  }
)
row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      grading2021_new = g,
      row.names = NULL
    )
  })
)

#make plot for averages per WHO grade
cairo_pdf(filename = "Top1000_avg_grade_new.pdf", width = 3.2, height = 7.5)
ggplot(row_means, 
       aes(grading2021_new, average, fill=grading2021_new))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#9EBCDA","#8C6BB1","#810F7C")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#split grade 2 by cluster
meta$cluster_grade <- paste(meta$cluster_1000_new, meta$grading2021_new, sep = "_")
cols_to_keep <- rownames(meta)[meta$grading2021_new == "2"]
beta_top1000_hyper <- beta_top1000[, cols_to_keep, drop = FALSE]
cluster_grades <- meta[cols_to_keep, "cluster_grade"]
row_means_wide <- sapply(
  split(colnames(beta_top1000_hyper), cluster_grades),
  function(sample_set) {
    rowMeans(beta_top1000_hyper[, sample_set, drop = FALSE])
  }
)

row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      grade_cluster = g,
      row.names = NULL
    )
  })
)
row_means$grade_cluster = factor(row_means$grade_cluster, levels = c("METHlow_2", "METHhigh_2"))

split_cols <- do.call(rbind, strsplit(as.character(row_means$grade_cluster), "_"))
split_cols <- data.frame(grade = split_cols[,1], 
                         cluster = split_cols[,2],
                         stringsAsFactors = FALSE)
row_means <- cbind(row_means, split_cols)
row_means$grade = factor(row_means$grade, levels = c("METHlow", "METHhigh"))
row_means$cluster = factor(row_means$cluster, levels = c("2"))


cairo_pdf(filename = "Top1000_average_grade-cluster_new.pdf", width = 2, height = 7.5)
ggplot(row_means, 
       aes(cluster, average, fill=interaction(cluster,grade), dodge=grade))+
  geom_point(position = position_jitterdodge(dodge.width = 0.85), alpha=0.95, color="grey60")+
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95))+
  scale_fill_manual(values=c("#8C6BB166","#8C6BB1")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()



#split grade 3 by cluster
meta$cluster_grade <- paste(meta$cluster_1000_new, meta$grading2021_new, sep = "_")
cols_to_keep <- rownames(meta)[meta$grading2021_new == "3"]
beta_top1000_hyper <- beta_top1000[, cols_to_keep, drop = FALSE]
cluster_grades <- meta[cols_to_keep, "cluster_grade"]
row_means_wide <- sapply(
  split(colnames(beta_top1000_hyper), cluster_grades),
  function(sample_set) {
    rowMeans(beta_top1000_hyper[, sample_set, drop = FALSE])
  }
)

row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      grade_cluster = g,
      row.names = NULL
    )
  })
)
row_means$grade_cluster = factor(row_means$grade_cluster, levels = c("METHlow_3", "METHhigh_3"))

split_cols <- do.call(rbind, strsplit(as.character(row_means$grade_cluster), "_"))
split_cols <- data.frame(grade = split_cols[,1], 
                         cluster = split_cols[,2],
                         stringsAsFactors = FALSE)
row_means <- cbind(row_means, split_cols)
row_means$grade = factor(row_means$grade, levels = c("METHlow", "METHhigh"))
row_means$cluster = factor(row_means$cluster, levels = c("3"))


cairo_pdf(filename = "Top1000_average_grade3-cluster_new.pdf", width = 2, height = 7.5)
ggplot(row_means, 
       aes(cluster, average, fill=interaction(cluster,grade), dodge=grade))+
  geom_point(position = position_jitterdodge(dodge.width = 0.85), alpha=0.95, color="grey60")+
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95))+
  scale_fill_manual(values=c("#810F7C66","#810F7C")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()






#groups
group <- meta2$MCconsensus
row_means_wide <- sapply(
  split(colnames(beta_top1000), group),
  function(sample_set) {
    rowMeans(beta_top1000[, sample_set, drop = FALSE])
  }
)
row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      group = g,
      row.names = NULL
    )
  })
)
row_means$group = factor(row_means$group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

#make plot for averages per WHO grade
cairo_pdf(filename = "Top1000_avg_group.pdf", width = 4, height = 7.5)
ggplot(row_means, 
       aes(group, average, fill=group))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("royalblue2","red3","forestgreen","darkorange2")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#split hypermetabolic by cluster
meta$cluster_group <- paste(meta$cluster_1000_new, meta$mol_group, sep = "_")
cols_to_keep <- rownames(meta)[meta$mol_group == "hypermetabolic"]
beta_top1000_hyper <- beta_top1000[, cols_to_keep, drop = FALSE]
cluster_groups <- meta[cols_to_keep, "cluster_group"]
row_means_wide <- sapply(
  split(colnames(beta_top1000_hyper), cluster_groups),
  function(sample_set) {
    rowMeans(beta_top1000_hyper[, sample_set, drop = FALSE])
  }
)

row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      group_cluster = g,
      row.names = NULL
    )
  })
)
row_means$group_cluster = factor(row_means$group_cluster, levels = c("METHlow_hypermetabolic", "METHhigh_hypermetabolic"))

split_cols <- do.call(rbind, strsplit(as.character(row_means$group_cluster), "_"))
split_cols <- data.frame(group = split_cols[,1], 
                         cluster = split_cols[,2],
                         stringsAsFactors = FALSE)
row_means <- cbind(row_means, split_cols)
row_means$group = factor(row_means$group, levels = c("METHlow", "METHhigh"))
row_means$cluster = factor(row_means$cluster, levels = c("hypermetabolic"))


cairo_pdf(filename = "Top1000_average_group-cluster.pdf", width = 2, height = 7.5)
ggplot(row_means, 
       aes(cluster, average, fill=interaction(cluster,group), dodge=group))+
  geom_point(position = position_jitterdodge(dodge.width = 0.95), alpha=0.95, color="grey60")+
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95))+
  scale_fill_manual(values=c("#228B2266","#228B22")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()






#split proliferative by cluster
meta$cluster_group <- paste(meta$cluster_1000_new, meta$mol_group, sep = "_")
cols_to_keep <- rownames(meta)[meta$mol_group == "proliferative"]
beta_top1000_hyper <- beta_top1000[, cols_to_keep, drop = FALSE]
cluster_groups <- meta[cols_to_keep, "cluster_group"]
row_means_wide <- sapply(
  split(colnames(beta_top1000_hyper), cluster_groups),
  function(sample_set) {
    rowMeans(beta_top1000_hyper[, sample_set, drop = FALSE])
  }
)

row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      group_cluster = g,
      row.names = NULL
    )
  })
)
row_means$group_cluster = factor(row_means$group_cluster, levels = c("METHlow_proliferative", "METHhigh_proliferative"))

split_cols <- do.call(rbind, strsplit(as.character(row_means$group_cluster), "_"))
split_cols <- data.frame(group = split_cols[,1], 
                         cluster = split_cols[,2],
                         stringsAsFactors = FALSE)
row_means <- cbind(row_means, split_cols)
row_means$group = factor(row_means$group, levels = c("METHlow", "METHhigh"))
row_means$cluster = factor(row_means$cluster, levels = c("proliferative"))


cairo_pdf(filename = "Top1000_average_group-cluster.pdf", width = 2, height = 7.5)
ggplot(row_means, 
       aes(cluster, average, fill=interaction(cluster,group), dodge=group))+
  geom_point(position = position_jitterdodge(dodge.width = 0.95), alpha=0.95, color="grey60")+
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95))+
  scale_fill_manual(values=c("#ee760066","#ee7600")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()





#risk
group <- meta$risk_score_new
row_means_wide <- sapply(
  split(colnames(beta_top1000), group),
  function(sample_set) {
    rowMeans(beta_top1000[, sample_set, drop = FALSE])
  }
)

row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      group = g,
      row.names = NULL
    )
  })
)

row_means$group = factor(row_means$group, levels = c("low", "intermediate", "high"))

#make plot for averages per risk group
cairo_pdf(filename = "Top1000_avg_risk_new.pdf", width = 3.2, height = 7.5)
ggplot(row_means, 
       aes(group, average, fill=group))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("dodgerblue2","purple1","firebrick1")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#split intermediate by cluster
meta$cluster_risk <- paste(meta$cluster_1000_new, meta$risk_score_new, sep = "_")
cols_to_keep <- rownames(meta)[meta$risk_score_new == "intermediate"]
beta_top1000_hyper <- beta_top1000[, cols_to_keep, drop = FALSE]
cluster_risk <- meta[cols_to_keep, "cluster_risk"]
row_means_wide <- sapply(
  split(colnames(beta_top1000_hyper), cluster_risk),
  function(sample_set) {
    rowMeans(beta_top1000_hyper[, sample_set, drop = FALSE])
  }
)

row_means <- do.call(
  rbind,
  lapply(colnames(row_means_wide), function(g) {
    data.frame(
      average = row_means_wide[, g],
      risk_cluster = g,
      row.names = NULL
    )
  })
)
row_means$risk_cluster = factor(row_means$risk_cluster, levels = c("METHlow_intermediate", "METHhigh_intermediate"))

split_cols <- do.call(rbind, strsplit(as.character(row_means$risk_cluster), "_"))
split_cols <- data.frame(risk = split_cols[,1], 
                         cluster = split_cols[,2],
                         stringsAsFactors = FALSE)
row_means <- cbind(row_means, split_cols)
row_means$risk = factor(row_means$risk, levels = c("METHlow", "METHhigh"))
row_means$cluster = factor(row_means$cluster, levels = c("intermediate"))

cairo_pdf(filename = "Top1000_average_risk-cluster_new.pdf", width = 2, height = 7.5)
ggplot(row_means, 
       aes(cluster, average, fill=interaction(cluster,risk), dodge=risk))+
  geom_point(position = position_jitterdodge(dodge.width = 0.95), alpha=0.95, color="grey60")+
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.95))+
  scale_fill_manual(values=c("#9D00FF66","#9D00FF")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

save.image()

###calculate mean global methylation levels across WHO grades and molecular groups
#discovery cohort
beta_disc = readRDS(file="Betas_discovery_preprocessed_corrected.rds")
sample_means_disc <- colMeans(beta_disc, na.rm = TRUE)
mean_df_disc <- data.frame(
  ID = names(sample_means_disc),
  mean   = as.numeric(sample_means_disc),
  row.names = NULL
)
mean_df_disc <- merge(mean_df_disc, meta[, c("ID", "mol_group", "grading2021_new")],
                 by = "ID", all.x = TRUE)
mean_df_disc$mol_group = factor(mean_df_disc$mol_group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
mean_df_disc$grading2021_new = factor(mean_df_disc$grading2021_new, levels = c("1", "2", "3"))

cairo_pdf(filename = "GlobaL_mean_discovery_grade.pdf", width = 6, height = 5)
ggplot(mean_df_disc, 
       aes(grading2021_new, mean, fill=grading2021_new))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("#9EBCDA","#8C6BB1","#810F7C")) +
  scale_y_continuous(breaks = seq(0.5,0.75,0.1), limits = c(0.5,0.75))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


cairo_pdf(filename = "Global_mean_discovery_molgroup.pdf", width = 6, height = 5)
ggplot(mean_df_disc, 
       aes(mol_group, mean, fill=mol_group))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("royalblue2","red3","forestgreen","darkorange2")) +
  scale_y_continuous(breaks = seq(0.5,0.75,0.1), limits = c(0.5,0.75))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()










###########UCSF

#chisquare for UCSF
#grade
meta_UCSF = read.csv(file="meta_UCSF_current.csv")

length(which(meta_UCSF$grade == "1" & meta_UCSF$cluster == "METHlow"))
length(which(meta_UCSF$grade == "2" & meta_UCSF$cluster == "METHlow"))
length(which(meta_UCSF$grade == "3" & meta_UCSF$cluster == "METHlow"))
length(which(meta_UCSF$grade == "1" & meta_UCSF$cluster == "METHhigh"))
length(which(meta_UCSF$grade == "2" & meta_UCSF$cluster == "METHhigh"))
length(which(meta_UCSF$grade == "3" & meta_UCSF$cluster == "METHhigh"))

dat1000 <- data.frame(
  c(366, 22),
  c(100, 42),
  c(22,13),
  row.names = c("METHlow", "METHhigh"),
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
df1000$clustering = factor(df1000$clustering, levels = c("METHlow", "METHhigh"))

test1000 <- fisher.test(table(df1000))

tiff("Chisquare_UCSF_grade.tiff", res=300,width=1000,height=2048)
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



length(which(meta_UCSF$mol_group == "Merlinintact" & meta_UCSF$cluster == "METHlow"))
length(which(meta_UCSF$mol_group == "Immuneenriched" & meta_UCSF$cluster == "METHlow"))
length(which(meta_UCSF$mol_group == "hypermetabolic" & meta_UCSF$cluster == "METHlow"))
length(which(meta_UCSF$mol_group == "proliferative" & meta_UCSF$cluster == "METHlow"))
length(which(meta_UCSF$mol_group == "Merlinintact" & meta_UCSF$cluster == "METHhigh"))
length(which(meta_UCSF$mol_group == "Immuneenriched" & meta_UCSF$cluster == "METHhigh"))
length(which(meta_UCSF$mol_group == "hypermetabolic" & meta_UCSF$cluster == "METHhigh"))
length(which(meta_UCSF$mol_group == "proliferative" & meta_UCSF$cluster == "METHhigh"))

dat1000 <- data.frame(
  c(187, 0),
  c(111, 0),
  c(156,16),
  c(34, 61),
  row.names = c("METHlow", "METHhigh"),
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
df1000$clustering = factor(df1000$clustering, levels = c("METHlow", "METHhigh"))
df1000$MCconsensus = factor(df1000$MCconsensus, levels =  c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

test1000 <- chisq.test(table(df1000))

tiff("Chisquare_UCSF_group.tiff", res=300,width=1000,height=2048)
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






#visualize average methylation of top 1000 PC1 probes across grades and molecular groups, also split by cluster

#get EPIC asociated probes for top1000 PC1 probes from EPICv2
meta_UCSF= read.csv(file="meta_UCSF_current.csv")
top1000_EPIC = as.vector(rownames(beta_top1000))
top1000_EPIC = mLiftOver(top1000_EPIC, "EPIC")
UCSF_beta = readRDS(file = "Betas_UCSF_preprocessed.rds")
UCSF_beta_top1000 = UCSF_beta[rownames(UCSF_beta) %in% top1000_EPIC,]
all(colnames(UCSF_beta_top1000) == meta_UCSF$ID)
meta_UCSF = meta_UCSF[match(colnames(UCSF_beta_top1000), meta_UCSF$ID),]
#grades
groups <- split(meta_UCSF$ID, meta_UCSF$grade)
group_means <- lapply(names(groups), function(g) {
  rowMeans(UCSF_beta_top1000[, groups[[g]], drop = FALSE], na.rm = TRUE)
})

group_labels <- c("one", "two", "three")
out_df <- do.call(rbind, lapply(seq_along(group_means), function(i) {
  data.frame(
    mean  = unname(group_means[[i]]),  
    group = group_labels[i]            
  )
}))
out_df$group = factor(out_df$group, levels = c("one", "two", "three"))
str(out_df)

cairo_pdf(filename = "Top1000_avg_grade_UCSF.pdf", width = 3.2, height = 7.5)
ggplot(out_df, 
       aes(group, mean, fill=group))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=1)+
  scale_fill_manual(values=c("#9EBCDA","#8C6BB1","#810F7C")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

##split by cluster for all grades
meta_UCSF$grade_cluster <- paste(meta_UCSF$grade, meta_UCSF$cluster, sep = "_")
groups <- split(meta_UCSF$ID, meta_UCSF$grade_cluster)
group_means <- lapply(groups, function(samples) {
  rowMeans(UCSF_beta_top1000[, samples, drop = FALSE], na.rm = TRUE)
})
out_df <- do.call(rbind, lapply(names(group_means), function(g) {
  data.frame(
    mean = unname(group_means[[g]]),          
    group_condition = g
  )
}))
out_df$group_condition = factor(out_df$group_condition, levels = c("1_METHlow", "1_METHhigh", "2_METHlow", "2_METHhigh","3_METHlow", "3_METHhigh"))
split_cols <- do.call(rbind, strsplit(as.character(out_df$group_condition), "_"))
split_cols <- data.frame(group = split_cols[,1], 
                         condition = split_cols[,2],
                         stringsAsFactors = FALSE)
out_df <- cbind(out_df, split_cols)
out_df$group = factor(out_df$group, levels = c("1","2","3"))
out_df$condition = factor(out_df$condition, levels = c("METHlow", "METHhigh"))


cairo_pdf(filename = "Top1000_avg_UCSF_grade_cluster.pdf", width = 6, height = 7.5)
ggplot(out_df, 
       aes(group, mean, fill=interaction(group,condition), dodge=condition))+
  geom_point(position = position_jitterdodge(0.4), alpha=0.95, color="grey60")+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#9EBCDA66","#8C6BB166","#810F7C66","#9EBCDA","#8C6BB1", "#810F7C")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()




#molecular groups
groups <- split(meta_UCSF$ID, meta_UCSF$mol_group)
str(meta_UCSF)
group_means <- lapply(names(groups), function(g) {
  rowMeans(UCSF_beta_top1000[, groups[[g]], drop = FALSE], na.rm = TRUE)
})

group_labels <- c("hypermetabolic", "Immuneenriched", "Merlinintact", "proliferative")
out_df <- do.call(rbind, lapply(seq_along(group_means), function(i) {
  data.frame(
    mean  = unname(group_means[[i]]),  
    group = group_labels[i]            
  )
}))
str(out_df)
out_df$group = factor(out_df$group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
str(out_df)

cairo_pdf(filename = "Top1000_avg_group_UCSF.pdf", width = 4, height = 7.5)
ggplot(out_df, 
       aes(group, mean, fill=group))+
  geom_point(position = position_jitter(0.2), color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=1)+
  scale_fill_manual(values=c("royalblue2","red3","forestgreen","darkorange2")) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()



##split by cluster for all grades
meta_UCSF$group_cluster <- paste(meta_UCSF$mol_group, meta_UCSF$cluster, sep = "_")
groups <- split(meta_UCSF$ID, meta_UCSF$group_cluster)
group_means <- lapply(groups, function(samples) {
  rowMeans(UCSF_beta_top1000[, samples, drop = FALSE], na.rm = TRUE)
})
out_df <- do.call(rbind, lapply(names(group_means), function(g) {
  data.frame(
    mean = unname(group_means[[g]]),          
    group_condition = g
  )
}))
out_df$group_condition = factor(out_df$group_condition, levels = c("Merlinintact_METHlow", "Merlinintact_METHhigh",
                                                                   "Immuneenriched_METHlow", "Immuneenriched_METHhigh",
                                                                   "hypermetabolic_METHlow", "hypermetabolic_METHhigh",
                                                                   "proliferative_METHlow", "proliferative_METHhigh"))
split_cols <- do.call(rbind, strsplit(as.character(out_df$group_condition), "_"))
split_cols <- data.frame(group = split_cols[,1], 
                         condition = split_cols[,2],
                         stringsAsFactors = FALSE)
out_df <- cbind(out_df, split_cols)
out_df$group = factor(out_df$group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
out_df$condition = factor(out_df$condition, levels = c("METHlow", "METHhigh"))

new_row <- data.frame(
  mean = 0.8,
  group_condition = "Immuneenriched_METHhigh",
  group = "Immuneenriched",
  condition = "METHhigh",
  stringsAsFactors = FALSE
)

out_df = rbind(out_df, new_row)

cairo_pdf(filename = "Top1000_avg_UCSF_group_cluster.pdf", width = 8, height = 7.5)
ggplot(out_df, 
       aes(group, mean, fill=interaction(group,condition), dodge=condition))+
  geom_point(position = position_jitterdodge(0.4), alpha=0.95, color="grey60")+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#436EEE66","#CD000066","#228B2266","#EE760066","#436EEE","#CD0000","#228B22","#EE7600" )) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits=c(0,1))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()






#Sankey plot after methylation cluster and molecular group prediction for UCSF cohort

UCSF_Sankey = read.csv(file="UCSF_For_Sankey.csv")

nodes <- data.frame(
  name = unique(c(
    UCSF_Sankey$WHO_grade,
    UCSF_Sankey$Molecular_group,
    UCSF_Sankey$Methylation_cluster
  )),
  stringsAsFactors = FALSE
)
nodes$label <- ""

links_1 <- UCSF_Sankey %>%
  group_by(WHO_grade, Molecular_group) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(WHO_grade, nodes$name) - 1,
    target = match(Molecular_group, nodes$name) - 1
  ) %>%
  select(source, target, value)

links_2 <- UCSF_Sankey %>%
  group_by(Molecular_group, Methylation_cluster) %>%
  summarise(value = sum(n), .groups = "drop") %>%
  mutate(
    source = match(Molecular_group, nodes$name) - 1,
    target = match(Methylation_cluster, nodes$name) - 1
  ) %>%
  select(source, target, value)

links <- bind_rows(links_1, links_2)

nodes$group <- NA_character_

# WHO grades
nodes$group[nodes$name == "1"] <- "Grade1"
nodes$group[nodes$name == "2"] <- "Grade2"
nodes$group[nodes$name == "3"] <- "Grade3"

# Molecular groups
nodes$group[nodes$name == "Merlinintact"]     <- "Merlin"
nodes$group[nodes$name == "Immuneenriched"]   <- "Immune"
nodes$group[nodes$name == "Hypermetabolic"]   <- "Hyper"
nodes$group[nodes$name == "Proliferative"]    <- "Prolif"

# Methylation clusters
nodes$group[nodes$name == "METHlow"]  <- "METHlow"
nodes$group[nodes$name == "METHhigh"] <- "METHhigh"

stopifnot(!any(is.na(nodes$group)))



my_colour_scale <- '
d3.scaleOrdinal()
  .domain([
    "Grade1","Grade2","Grade3",
    "Merlin","Immune","Hyper","Prolif",
    "METHlow","METHhigh"
  ])
  .range([
    "#9EBCDA", "#8C6BB1", "#810F7C",   // grades
    "	#436eee", "#cd0000", "#228b22", "#ee7600", // molecular
    "#008b8b", "#8b2252"              // methylation
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
  nodeWidth = 30,
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
  nodeWidth = 30,
  sinksRight = FALSE
)

saveWidget(sankey, "sankey_UCSF.html", selfcontained = TRUE)
webshot(
  "sankey_UCSF.html",
  file = "sankey_UCSF.pdf",
  vwidth = 1600,
  vheight = 1000
)



###global mean for UCSF cohort
#UCSF
meta_UCSF= read.csv(file="meta_UCSF_current.csv")
beta_UCSF = readRDS(file="Betas_UCSF_preprocessed.rds")

sample_means_UCSF <- colMeans(beta_UCSF, na.rm = TRUE)
mean_df_UCSF <- data.frame(
  ID = names(sample_means_UCSF),
  mean   = as.numeric(sample_means_UCSF),
  row.names = NULL
)
mean_df_UCSF <- merge(mean_df_UCSF, meta_UCSF[, c("ID", "mol_group", "grade")],
                      by = "ID", all.x = TRUE)
mean_df_UCSF$mol_group = factor(mean_df_UCSF$mol_group, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
mean_df_UCSF$grade = factor(mean_df_UCSF$grade, levels = c("1", "2", "3"))

cairo_pdf(filename = "GlobaL_mean_UCSF_grade.pdf", width = 6, height = 5)
ggplot(mean_df_UCSF, 
       aes(grade, mean, fill=grade))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("#9EBCDA","#8C6BB1","#810F7C")) +
  scale_y_continuous(breaks = seq(0.5,0.75,0.1), limits = c(0.5,0.75))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


cairo_pdf(filename = "Global_mean_UCSF_molgroup.pdf", width = 6, height = 5)
ggplot(mean_df_UCSF, 
       aes(mol_group, mean, fill=mol_group))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("royalblue2","red3","forestgreen","darkorange2")) +
  scale_y_continuous(breaks = seq(0.5,0.75,0.1), limits = c(0.5,0.75))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


####check correlation of average methylation of top1000 probes used for clustering with Ki67 fraction of corresponding tumor tissue

Ki67_cor = beta_top1000
Ki67_cor["average", ] <- colMeans(Ki67_cor, na.rm = TRUE)
Ki67_cor = as.data.frame(t(Ki67_cor))
Ki67_cor = Ki67_cor[,-(1:1000), drop = FALSE]

meta = read.csv(file="meta_discovery_current.csv", header = T)
rownames(meta) = meta$ID
all(rownames(Ki67_cor) == rownames(meta))
meta = meta[match(rownames(Ki67_cor),rownames(meta)),]
all(rownames(Ki67_cor) == rownames(meta))

Ki67_cor$Ki67 = meta$Ki67



#make the plot
ggscatter(
  Ki67_cor,
  x = colnames(Ki67_cor)[1],
  y = colnames(Ki67_cor)[2],
  add = "reg.line",
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson"
)





save.image()



















