library(sesame)
library(caret)
library(MLeval)
library(ggpubr)
library(factoextra)
library(pROC)

#####get data from discovery, longitudinal, UCSF, and Nassiri data and subset to EPIC-matched probes present in all datasets after masking
###for predictions, use only probes that are present in all datasets
###this ensures consistency in predicting both methylation clusters and molecular groups across different cohorts

###get UCSF data and pre-process, keep only complete case probes
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/Raleigh_cohort"
list.files(idat_dir)
meta_Raleigh = read.csv(file="meta_Raleigh.csv")

beta.UCSF = openSesame(idat_dir)
all(colnames(beta.UCSF)==meta_Raleigh$Basename)
beta.UCSF = beta.UCSF[,match(meta_Raleigh$Basename, colnames(beta.UCSF))]
all(colnames(beta.UCSF)==meta_Raleigh$Basename)
colnames(beta.UCSF) = meta_Raleigh$ID
beta.UCSF = beta.UCSF[complete.cases(beta.UCSF),]
saveRDS(beta.UCSF, file="Betas_UCSF_preprocessed.rds")

###get Nassiri data and pre-process, keep only complete case probes
idat_dir = "/Volumes/MAC_backup/Nassiri"
list.files(idat_dir)
meta_Nassiri = read.csv(file="Nassiri_targets.csv")

beta.Nassiri = openSesame(idat_dir)
all(colnames(beta.Nassiri)==meta_Nassiri$Basename)
colnames(beta.Nassiri) = meta_Nassiri$ID
beta.Nassiri = beta.Nassiri[complete.cases(beta.Nassiri),]
saveRDS(beta.Nassiri, file="Betas_Nassiri_preprocessed.rds")

###get discovery data that is already pre-processed, only complete case probes, that math EPIC probes
beta.discovery = readRDS(file="Betas_discovery_preprocessed_corrected_EPIC_match.RDS")

meta_discovery = read.csv(file="meta_discovery_current.csv")


###get longitudinal data and pre-process, keep only complete case probes
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/T26_longitudinal_only"
list.files(idat_dir)
meta_long = read.csv(file="targets_longitudinal.csv")

beta.long = openSesame(idat_dir)
all(colnames(beta.long)==meta_long$Basename)
beta.long = beta.long[,match(meta_long$Basename, colnames(beta.long))]
colnames(beta.long) = meta_long$ID
dim(beta.long)
beta.long = mLiftOver(beta.long, "EPIC")
beta.long = beta.long[complete.cases(beta.long),]
saveRDS(beta.long, file="Betas_longitudinal_preprocessed.rds")


####use PCA to identify components associated with condition of interest, and get ranking of probes for indicated components
###for methylation cluster, use discovery cohort and add clustering information from hierarchical clustering
###for molecular groups, use Nassiri cohort and check associations

#discovery for METH cluster
disc_10k = as.data.frame(beta.discovery)
str(disc_10k)
disc_10k$var = apply(disc_10k,1,var)
disc_10k <- disc_10k[order(disc_10k$var, decreasing = TRUE),]
disc_10k = disc_10k[,-232]
disc_10k = disc_10k[1:10000,]
disc_10k_pca = prcomp(t(disc_10k), scale. = T)

df_disc <- cbind(meta_discovery, disc_10k_pca$x[,1:4])
df_disc$cluster_1000_new = factor(df_disc$cluster_1000_new)

tiff("PCA_METHcluster_discovery_PC1-2.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_disc, x = "PC1", y = "PC2",
  color = "cluster_1000_new", size = 5, alpha = 0.8,
  palette = c("violetred4","cyan4"),
  margin.params = list(fill = "cluster_1000_new", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_METHcluster_discovery_PC3-4.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_disc, x = "PC3", y = "PC4",
  color = "cluster_1000_new", size = 5, alpha = 0.8,
  palette = c("violetred4","cyan4"),
  margin.params = list(fill = "cluster_1000_new", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

###get weigthed pairwise mean overlap of cluster assignment overlap for the first 4 PCs (run iteratively)
#for METHlow and METHhigh
pc_scores <- disc_10k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta_discovery$cluster_1000_new
)
pc1_groups <- split(pc1_df$PC1, pc1_df$Group)
g1 = pc1_groups[["METHlow"]]
g2 = pc1_groups[["METHhigh"]]

#get weigthed pairwise mean overlap 
groups <- list(G1 = g1, G2 = g2)
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

#get probes ranked by PC1
#get names for top contributing probes for PC1
var_10k = get_pca_var(disc_10k_pca)
head(var_10k$contrib,10)
var_contrib_10k = as.data.frame(var_10k$contrib)
var_contrib_sortPC1 = var_contrib_10k[order(var_contrib_10k$Dim.1, decreasing = TRUE),]
METHcluster_PC1_rank = data.frame(PC1=rownames(var_contrib_sortPC1))




#Nassiri for molecular groups
nassiri_10k = as.data.frame(beta.Nassiri)
nassiri_10k$var = apply(nassiri_10k,1,var)
nassiri_10k <- nassiri_10k[order(nassiri_10k$var, decreasing = TRUE),]
nassiri_10k = nassiri_10k[,-122]
nassiri_10k = nassiri_10k[1:10000,]
nassiri_10k_pca = prcomp(t(nassiri_10k), scale. = T)

df_nassiri <- cbind(meta_Nassiri, nassiri_10k_pca$x[,1:4])
df_nassiri$MGgroup = factor(df_nassiri$MGgroup, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))


#check number of meaningful components
fviz_screeplot(nassiri_10k_pca, addlabel=T)

#make PCA plot
tiff("PCA_MolGroup_Nassiri_PC1-2.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_nassiri, x = "PC1", y = "PC2",
  color = "MGgroup", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MGgroup", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()

tiff("PCA_MolGroup_Nassiri_PC3-4.tiff", res=300,width=2460,height=2040)
ggscatterhist(
  df_nassiri, x = "PC3", y = "PC4",
  color = "MGgroup", size = 5, alpha = 0.8,
  palette = c("royalblue2","red3","forestgreen","darkorange2"),
  margin.params = list(fill = "MGgroup", color = "black", size = 0.3),
  margin.plot.size = 0.6
)
dev.off()


###get weigthed pairwise mean overlap of molecular group assignment overlap for the first 4 PCs (run iteratively)
#for METHlow and METHhigh
pc_scores <- nassiri_10k_pca$x
pc1_df <- data.frame(
  Sample = rownames(pc_scores),
  PC1 = pc_scores[, "PC4"],
  Group = meta_Nassiri$MGgroup
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

#get probes ranked by PC1 and PC3
#get names for top contributing probes for PC1
var_10k = get_pca_var(nassiri_10k_pca)
head(var_10k$contrib,10)
var_contrib_10k = as.data.frame(var_10k$contrib)
var_contrib_sortPC1 = var_contrib_10k[order(var_contrib_10k$Dim.3, decreasing = TRUE),]
MolGroup_PC1_rank = data.frame(PC1=rownames(var_contrib_sortPC1))
MolGroup_PC2_rank = data.frame(PC2=rownames(var_contrib_sortPC1))
MolGroup_PC3_rank = data.frame(PC3=rownames(var_contrib_sortPC1))



###subset to probes which are present in all datasets
common_rows <- Reduce(intersect, list(rownames(beta.discovery), rownames(beta.long),rownames(beta.UCSF), rownames(beta.Nassiri)))

#subset all matrices to common rows
beta.discovery.sub = beta.discovery[common_rows, , drop = FALSE]
beta.long.sub = beta.long[common_rows, , drop = FALSE]
beta.UCSF.sub = beta.UCSF[common_rows, , drop = FALSE]
beta.Nassiri.sub = beta.Nassiri[common_rows, , drop = FALSE]



############model METH clustering in discovery data and predict meningioma cluster in UCSF and longitudinal data
#subset discovery by the first 1000 probes associating to PC1, which separates METH clusters in the discovery cohort
pc1_disc_probes <- METHcluster_PC1_rank$PC1
available <- pc1_disc_probes[pc1_disc_probes %in% rownames(beta.discovery.sub)]
selected_1000 <- head(available, 1000)
beta.discovery.sub.sub = beta.discovery.sub[rownames(beta.discovery.sub) %in% selected_1000,]

#transform dataset and include meth group data
dataset <- t(beta.discovery.sub.sub)
dataset <- as.data.frame(dataset)
all(rownames(dataset)==meta_discovery$ID)
condition=factor(meta_discovery$cluster_1000_new, levels = c("METHlow", "METHhigh"))
dataset <- cbind(dataset,condition)

#set up classifier 
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

set.seed(7)
fit.rf <- train(condition~., data=dataset, method="rf", metric=metric, trControl=control)
set.seed(7)
fit.svm <- train(condition~., data=dataset, method="svmRadial", metric=metric, trControl=control)
set.seed(7)
fit.knn <- train(condition~., data=dataset, method="knn", metric=metric, trControl=control)
set.seed(7)
fit.PAM <- train(condition~., data=dataset, method="pam", metric=metric, trControl=control)
set.seed(7)
fit.gbm <- train(condition~., data=dataset, method="gbm", metric=metric, trControl=control)

results <- resamples(list(knn=fit.knn, PAM=fit.PAM, svm=fit.svm, rf=fit.rf, gbm=fit.gbm))
dotplot(results)
#highest accuracy and Kappa for svm

###generate ROC curve for model evaluation to predict cluster assignment in discovery
control <- trainControl(method="cv", number=10, summaryFunction = twoClassSummary, classProbs = T, savePredictions = T)
fit.disc.clust.svm.roc <- train(condition~., data=dataset, method="svmRadial", metric="ROC", trControl=control, verbose=FALSE)
fit.disc.clust.rf.roc <- train(condition~., data=dataset, method="rf", metric="ROC", trControl=control, verbose=FALSE)
fit.disc.clust.gbm.roc <- train(condition~., data=dataset, method="gbm", metric="ROC", trControl=control, verbose=FALSE)

cairo_pdf(filename = "ROC_discovery_cluster_prediction_severalmodels_crossvalidation.pdf", width = 8, height = 7)
ROC_disccovery_cluster_prediction_several = evalm(list(fit.disc.clust.svm.roc, fit.disc.clust.rf.roc,fit.disc.clust.gbm.roc),
                                                  gnames=c("svm","rf","gbm"),
                                                  cols=c("red3", "green4", "blue"),plots = "r", rlinethick = 1.5, fsize = 8)
dev.off()


#predict methylation cluster cluster for UCSF data
beta.UCSF.sub = t(beta.UCSF.sub)
predict.UCSF <- predict(fit.svm, beta.UCSF.sub)
all(rownames(beta.UCSF.sub)==meta_Raleigh$ID)
meta_Raleigh$cluster = predict.UCSF
write.csv(meta_Raleigh, file="meta_UCSF_current.csv")
meta_UCSF = read.csv(file="meta_UCSF_current.csv")

#predict methylation cluster for longitudinal data
beta.long.sub = t(beta.long.sub)
predict.long <- predict(fit.svm, beta.long.sub)
all(rownames(beta.long.sub)==meta_long$ID)
meta_long$cluster = predict.long
write.csv(meta_long, file="meta_longitudinal_current.csv")







####predict molecular group from Nassiri data in discovery, longitudinal, and UCSF cohort
#select top 1000 probes from both PC1 and PC3 to capture differences in molecular groups

get_top_n_present <- function(ranked_probes, available_probes, n = 1000) {
  present <- ranked_probes[ranked_probes %in% available_probes]
  head(present, n)
}
available <- rownames(beta.Nassiri.sub)

pc1_1000 <- get_top_n_present(MolGroup_PC1_rank$PC1, available, 1000)
pc3_1000 <- get_top_n_present(MolGroup_PC3_rank$PC3, available, 1000)
selected_probes <- unique(c(pc1_1000, pc3_1000))
beta.Nassiri.sub.sub <- beta.Nassiri.sub[rownames(beta.Nassiri.sub) %in% selected_probes, ]

#transform dataset and include meth group data
dataset <- t(beta.Nassiri.sub.sub)
dataset <- as.data.frame(dataset)
all(rownames(dataset)==meta_Nassiri$ID)
condition=factor(meta_Nassiri$MGgroup, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
dataset <- cbind(dataset,condition)

#set up classifier 
control <- trainControl(method="cv", number=10)
metric <- "Accuracy"

set.seed(7)
fit.rf <- train(condition~., data=dataset, method="rf", metric=metric, trControl=control)
set.seed(7)
fit.svm <- train(condition~., data=dataset, method="svmRadial", metric=metric, trControl=control)
set.seed(7)
fit.knn <- train(condition~., data=dataset, method="knn", metric=metric, trControl=control)
set.seed(7)
fit.PAM <- train(condition~., data=dataset, method="pam", metric=metric, trControl=control)
set.seed(7)
fit.gbm <- train(condition~., data=dataset, method="gbm", metric=metric, trControl=control)

results <- resamples(list(knn=fit.knn, PAM=fit.PAM, svm=fit.svm, rf=fit.rf, gbm=fit.gbm))
dotplot(results)
#highest accuracy for svm


#ROC curve for multi-class classifier
#micro-averaged ROC curves

control <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  savePredictions = "final"
)

#for svm
set.seed(7)
model <- train(condition~., data=dataset, method="svmRadial", metric=metric, trControl=control)
probs <- predict(model, type = "prob")
y_true <- model$trainingData$.outcome
classes <- colnames(probs)
y_bin <- sapply(classes, function(cl) {
  as.numeric(y_true == cl)
})
y_micro <- as.vector(y_bin)
p_micro <- as.vector(as.matrix(probs))
roc_micro_A <- roc(
  response = y_micro,
  predictor = p_micro
)


#for rf
set.seed(7)
model <- train(condition~., data=dataset, method="rf", metric=metric, trControl=control)
probs <- predict(model, type = "prob")
y_true <- model$trainingData$.outcome
classes <- colnames(probs)
y_bin <- sapply(classes, function(cl) {
  as.numeric(y_true == cl)
})
y_micro <- as.vector(y_bin)
p_micro <- as.vector(as.matrix(probs))
roc_micro_B <- roc(
  response = y_micro,
  predictor = p_micro
)

#for gbm
set.seed(7)
model <- train(condition~., data=dataset, method="gbm", metric=metric, trControl=control)
probs <- predict(model, type = "prob")
y_true <- model$trainingData$.outcome
classes <- colnames(probs)
y_bin <- sapply(classes, function(cl) {
  as.numeric(y_true == cl)
})
y_micro <- as.vector(y_bin)
p_micro <- as.vector(as.matrix(probs))
roc_micro_C <- roc(
  response = y_micro,
  predictor = p_micro
)


roc_list <- list(
  Model_A = roc_micro_A,
  Model_B = roc_micro_B,
  Model_C = roc_micro_C
)

cols <- c("red3", "green4", "blue")


cairo_pdf(filename = "Microaveraged_ROC_curves_Nassirir_molecular_groups.pdf", width = 8, height = 7)
plot(roc_list[[1]], col = cols[1], legacy.axes = TRUE, xlim = c(1, 0), xaxs = "i")
lines(roc_list[[2]], col = cols[2])
lines(roc_list[[3]], col = cols[3])
legend(
  "bottomright",
  legend = paste0(
    names(roc_list),
    " (AUC = ",
    round(sapply(roc_list, auc), 3),
    ")"
  ),
  col = cols,
  lwd = 2,
  bty = "n"
)
dev.off()



#predict mol groups in discovery using svm
beta.discovery.sub = t(beta.discovery.sub)
predict.discovery <- predict(fit.svm, beta.discovery.sub)
all(rownames(beta.discovery.sub)==meta_discovery$ID)
meta_discovery$mol_group = predict.discovery
write.csv(meta_discovery, file="meta_discovery_current.csv")

#predict mol groups in longitudinal data using svm
predict.long <- predict(fit.svm, beta.long.sub)
all(rownames(beta.long.sub)==meta_long$ID)
meta_long$mol_group = predict.long
write.csv(meta_long, file="meta_long_current.csv")

#predict mol groups in UCSF data using svm

predict.UCSF <- predict(fit.svm, beta.UCSF.sub)
all(rownames(beta.UCSF.sub)==meta_UCSF$ID)
meta_UCSF$mol_group = predict.UCSF
write.csv(meta_UCSF, file="meta_UCSF_current.csv")

save.image()



