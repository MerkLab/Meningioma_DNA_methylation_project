library(sesame)
library(ggplot2)
library(minfi)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(vegan)
library(ggplot2)
library(pals)
library(factoextra)
library(limma)

####sample-level QC
#get data from all preselected idat files
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/T26_preselect"
list.files(idat_dir)

#minfi library only for unprocessed data, comparison methlyated to unmethylated intensities
#read in meta data
targets <- read.metharray.sheet(idat_dir, pattern="meta_preselect.csv")
targets
rgSet <- read.metharray.exp(base = idat_dir, targets=targets)
rgSet
sampleNames(rgSet) <- targets$ID
rgSet
annotation(rgSet) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
rgSet

#make methylset, no normalization
MSet <- preprocessRaw(rgSet)

#get QC
qc <- getQC(MSet)
rownames(qc)
plotQC(qc,badSampleCutoff = 10.1)
#samples suggested to be of inferior quality: T26_0225, T26_0285, T26_0282

###use sesame pipeline for preprocessing
#path to idats and meta data
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/T26_preselect"
list.files(idat_dir)
setwd("/Volumes/MAC_backup/Tübingen_MNG_datasets/T26_preselect")
meta = read.csv((file="meta_preselect.csv"))
rownames(meta) = meta$X
str(meta)

#get beta values
beta.preselect = openSesame(idat_dir)
all(colnames(beta.preselect) %in% meta$Basename)
beta.preselect <- beta.preselect[, meta$Basename]
all(colnames(beta.preselect) == meta$Basename)
colnames(beta.preselect) = meta$ID
head(beta.preselect)

#delete inferior samples from beta and meta based on intensities
cols_to_remove = c("T26_0225", "T26_0285", "T26_0282")
beta.preselect.flt = beta.preselect[, !colnames(beta.preselect) %in% cols_to_remove]
dim(beta.preselect.flt)

rows_to_remove = c("T26_0225", "T26_0285", "T26_0282")
meta.flt = meta[!rownames(meta) %in% rows_to_remove,]

#check beta value distribution

samples <- colnames(beta.preselect.flt)
densityPlot(beta.preselect[, samples, drop = FALSE])
densityPlot(beta.preselect[, "T26_0065", drop = FALSE])
densityPlot(beta.preselect[, "T26_0005", drop = FALSE])
#exclude T26_0065 for uncommon beta value distribution

#delete inferior samples from beta and meta based beta distribution
cols_to_remove = c("T26_0065")
beta.preselect.flt = beta.preselect.flt[, !colnames(beta.preselect.flt) %in% cols_to_remove]
dim(beta.preselect.flt)

rows_to_remove = c("T26_0065")
meta.flt = meta.flt[!rownames(meta.flt) %in% rows_to_remove,]


#check sample mean of detP

#get detP values from sesame pOOBAH
pvalspOOBAH = openSesame(idat_dir, func = pOOBAH, return.pval = TRUE)
all(colnames(pvalspOOBAH) %in% meta$Basename)
pvalspOOBAH <- pvalspOOBAH[, meta$Basename]
all(colnames(pvalspOOBAH) == meta$Basename)
colnames(pvalspOOBAH) = meta$ID
head(pvalspOOBAH)
dim(pvalspOOBAH)

#delete previously excluded samples
cols_to_remove = c("T26_0225", "T26_0285", "T26_0282","T26_0065")
pvalspOOBAH.flt = pvalspOOBAH[, !colnames(pvalspOOBAH) %in% cols_to_remove]
dim(pvalspOOBAH.flt)
all(colnames(pvalspOOBAH.flt) == meta.flt$ID)

#get sample mean detP and graph
colMean = sort(colMeans(pvalspOOBAH.flt, na.rm = T))
barplot(colMean, ylim = c(0, 0.06), col="darkblue")
abline(h=0.05, col="red")
#remove 2 samples based on mean detP pOOBAH > 0.05 (T26_0184 + T26_0201)
#final discovery cohort has n=231 meningioma samples of good quality




######show effect of probe filtering

###sex chromosome associated probes effect on clustering
#read in data (n=231 MNG)
idat_dir = "/Users/lab/Desktop/Meningioma/data/T26_discovery"
targets <- read.metharray.sheet(idat_dir, pattern="T26_discovery_targets.csv")

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

mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

#make MethylSet a GenomicRatioSet
mSetSqFlt = mapToGenome(mSetSqFlt)
mSetSqFlt = ratioConvert(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# tag sex chromosome probes for removal
data("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

keep <- !(featureNames(mSetSqFlt) %in% annoEPICv2$Name[annoEPICv2$chr %in% 
                                                         c("chrX","chrY")])
table(keep)

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


#Probe-level QC for final discovery
#path to idats
idat_dir = "/Volumes/MAC_backup/Tübingen_MNG_datasets/T26_discovery"
list.files(idat_dir)
meta_discovery = read.csv((file="meta_discovery.csv"))

#Sesame QC pipeline for deteciton success
qcs = openSesame(idat_dir, prep="", func=sesameQC_calcStats)
qc_df = do.call(rbind, lapply(qcs, as.data.frame))
all(rownames(qc_df) == meta_discovery$Basename)

rownames(qc_df) = meta_discovery$ID
rownames(qc_df)

#make graph
det_suc = qc_df$num_dt
names(det_suc) = rownames(qc_df)
det_suc = sort(det_suc)
det_suc

df_det_suc <- data.frame(
  id    = seq_along(det_suc),   
  value = det_suc
)

ggplot(df_det_suc, aes(x = id, y = value)) +
  geom_col(fill = "black") +
  scale_x_continuous(
    breaks = seq(0, 250, 50),  
    expand = c(0, 0)
  ) +
  labs(x = "Sample index (sorted)", y = "Probes detection success") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_classic() 


##determine number of probes with detP > 0.05 per sample
pvalspOOBAH = openSesame(idat_dir, func = pOOBAH, return.pval = TRUE)
head(pvalspOOBAH)
all(colnames(pvalspOOBAH) == meta_discovery$Basename)
colnames(pvalspOOBAH) = meta_discovery$ID

probes_detP = colSums(pvalspOOBAH > 0.05,na.rm = TRUE)
probes_detP
write.csv(probes_detP, file = "Probes_excluded_detP.csv")

#make plot
probes_detP_sorted <- sort(probes_detP)

df_detP <- data.frame(
  id    = seq_along(probes_detP_sorted),   
  value = probes_detP_sorted
)

ggplot(df_detP, aes(x = id, y = value)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(
    breaks = seq(0, 250, 50),  
    expand = c(0, 0)
  ) +
  labs(x = "Sample index (sorted)", y = "Probes detP > 0.05") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_classic() 


#determine total number of probes masked across cohort
beta.discovery = openSesame(idat_dir)
colnames(beta.discovery)
all(colnames(beta.discovery) == meta_discovery$Basename)
colnames(beta.discovery) = meta_discovery$ID

rows_with_na <- rowSums(is.na(beta.discovery)) > 0
sum(rows_with_na) #401.311 probes masked in total

#determine number of probes masked per sample
masked_per_sample = colSums(is.na(beta.discovery))
masked_per_sample

#make plot
probes_mask_sorted <- sort(masked_per_sample)

df_mask <- data.frame(
  id    = seq_along(probes_mask_sorted),   
  value = probes_mask_sorted
)

ggplot(df_mask, aes(x = id, y = value)) +
  geom_col(fill = "darkorchid2") +
  scale_x_continuous(
    breaks = seq(0, 250, 50),  
    expand = c(0, 0)
  ) +
  labs(x = "Sample index (sorted)", y = "Probes detP > 0.05") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_classic() 


####merged density distribution
dim(beta.discovery) #937.690 probes in total

#remove sex chromosome probes
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

keep <- !(rownames(beta.discovery) %in% annoEPICv2$Name[annoEPICv2$chr %in% 
                                                          c("chrX","chrY")])
table(keep)
beta.discovery.sub = beta.discovery[keep,]
dim(beta.discovery.sub) 

#remove all probes which are masked in at least one case
beta.discovery.sub = beta.discovery.sub[complete.cases(beta.discovery.sub),]
dim(beta.discovery.sub) #527.352

#no prep beta values as contrast, but also remove sex probes
beta.noPrep = openSesame(idat_dir, prep = "", mask= FALSE)
dim(beta.noPrep)

beta.noPrep.sub= beta.noPrep[keep,]
dim(beta.noPrep.sub)

beta.noPrep.sub = beta.noPrep.sub[complete.cases(beta.noPrep.sub),] #no NAs in the data
dim(beta.noPrep.sub) #913.178 probes in total


#make density plots
densityPlot(beta.noPrep.sub)
densityPlot(beta.discovery.sub)

save.image()



#####check effects of technical factors (batches of DNA isolation, methylation arrays, tumor content, and age of FFPE block)

#get data
df = beta.discovery.sub

#perform principal component analysis
pca_df = prcomp(t(df), scale. = T)
meta_pca = meta_discovery[,c("ID", "DNA_batch", "Array_batch", "Tumor_content", "FFPE_age")]

pca_final = cbind(meta_pca, pca_df$x[,1:10])
str(pca_final)
pca_final$DNA_batch <- factor(
  pca_final$DNA_batch,
  levels = paste0("DNA_", sort(as.numeric(sub("DNA_", "", levels(pca_final$DNA_batch)))))
)
pca_final$DNA_batch <- as.factor(pca_final$DNA_batch)
pca_final$Array_batch <- as.factor(pca_final$Array_batch)
pca_final$Tumor_content <- as.factor(pca_final$Tumor_content)



#####use permutational test for homogeneity of multivariate dispersions (PERMDISP)
##if <0.05, unequal dispersions, permanova significance might be inflated, else permanova reliable
pc_mat <- pca_final[, paste0("PC", 1:10)]

# Euclidean distances in PC space
dist_mat <- dist(pc_mat, method = "euclidean")

# PERMDISP for DNA_batch
disp_dna <- betadisper(dist_mat, pca_final$DNA_batch)
permutest(disp_dna, permutations = 999) #dispersion homogenous, PERMANOVA reliable

# PERMDISP for Array_batch
disp_array <- betadisper(dist_mat, pca_final$Array_batch)
permutest(disp_array, permutations = 999) #dispersion homogenous, PERMANOVA reliable

# PERMDISP for Tumor_content
disp_content <- betadisper(dist_mat, pca_final$Tumor_content)
permutest(disp_content, permutations = 999) #slighlty significant (0.034), permanova might be driven by heterogeneity

# PERMDISP for age does not work, continous factor
#calculate distance of each sample to centroid and check correlation with age
centroid <- colMeans(pc_mat)
dist_to_centroid <- apply(pc_mat, 1, function(x) sqrt(sum((x - centroid)^2)))
cor.test(dist_to_centroid, pca_final$FFPE_age)


#show dispersion
boxplot(disp_dna)
boxplot(disp_array)
boxplot(disp_content)


#PERMANOVA test

# DNA batch
adonis2(dist_mat ~ DNA_batch, data = pca_final, permutations = 999, method = "euclidean")

# Array batch
adonis2(dist_mat ~ Array_batch, data = pca_final, permutations = 999, method = "euclidean")

# Tumor content
adonis2(dist_mat ~ Tumor_content, data = pca_final, permutations = 999, method = "euclidean")

# FFPE age
adonis2(dist_mat ~ FFPE_age, data = pca_final, permutations = 999, method = "euclidean")

#results: DNA-isolation, Array_run and age of tissue have a significant effect on multivariate structure!!!



###Scree plot for percent variances for uncorrected data

tiff("Scree_all_probes.tiff", res=300,width=2560,height=2048)
fviz_screeplot(pca_df, addlabels = T)
dev.off()

tiff("Scree_all_probes_clean.tiff", res=300,width=2560,height=2048)
fviz_screeplot(pca_df)
dev.off()


##visualize significant factors in PCA plots before correction
#DNA_batch
my_colors <- polychrome()
# Subset to the number of levels for DNA
my_colors_DNA <- polychrome()[1:length(levels(pca_final$DNA_batch))]
length(my_colors_DNA)
names(my_colors_DNA) <- levels(pca_final$DNA_batch)

tiff("PCA_DNA_batch_before_correction_legend.tiff", res=300,width=2460,height=1648)
ggplot(pca_final, aes(x=PC1, y=PC2, color = DNA_batch))+
  geom_jitter(height = 0.1, width = 0.1, size=3,show.legend = FALSE) +
  scale_color_manual(values=my_colors_DNA)+
  theme_classic()
dev.off()


#Array_batch
# Subset to the number of levels for DNA
my_colors_Array <- polychrome()[1:length(levels(pca_final$Array_batch))]
length(my_colors_Array)
names(my_colors_Array) <- levels(pca_final$Array_batch)

tiff("PCA_Array_batch_before_correction_legend.tiff", res=300,width=2460,height=1648)
ggplot(pca_final, aes(x=PC1, y=PC2, color = Array_batch))+
  geom_jitter(height = 0.1, width = 0.1, size=3,show.legend = FALSE) +
  scale_color_manual(values=my_colors_Array)+
  theme_classic()
dev.off()


#FFPE_Age
tiff("PCA_Age_batch_before_correction_legend.tiff", res=300,width=2460,height=1648)
ggplot(pca_final, aes(x=PC1, y=PC2, color = FFPE_age))+
  geom_jitter(height = 0.1, width = 0.1, size=3,show.legend = FALSE) +
  theme_classic()
dev.off()


######use limma removeBatchEffect to account for batch effects
##use subset of high quality probe beta values
#Ensure metadata rows match beta matrix columns
beta.discovery.sub
dim(beta.discovery.sub)
rownames(meta_discovery) = meta_discovery$ID
all(rownames(meta_discovery) == colnames(beta.discovery.sub))


#####for correction, convert to M values
#functions to convert B to M, and M to B
beta_to_m <- function(beta) {
  beta[beta <= 0] <- 1e-6
  beta[beta >= 1] <- 1 - 1e-6
  log2(beta / (1 - beta))
}

m_to_beta <- function(M) { 2^M / (1 + 2^M) }

M <- beta_to_m(beta.discovery.sub)


m_to_beta <- function(M) { 2^M / (1 + 2^M) }

#make meta data categorical factors
str(meta_discovery)
meta_discovery$gender = factor(meta_discovery$gender)
meta_discovery$status = factor(meta_discovery$status)
meta_discovery$grading2021 = factor(meta_discovery$grading2021)
meta_discovery$risk_score = factor(meta_discovery$risk_score)
meta_discovery$MCconsensus = factor(meta_discovery$MCconsensus)
meta_discovery$priorRT = factor(meta_discovery$priorRT)
meta_discovery$adjuvantRT = factor(meta_discovery$adjuvantRT)

#define design matrix to define biological factors to keep
design <- model.matrix(~ gender + status + grading2021 + risk_score + MCconsensus + priorRT + adjuvantRT, data = meta_discovery)

#define batches to remove
batch_factors <- meta_discovery[, c("DNA_batch", "Array_batch", "FFPE_age")]
batch_factors$DNA_batch = factor(batch_factors$DNA_batch)
batch_factors$Array_batch = factor(batch_factors$Array_batch)
str(batch_factors)

#use limma to remove batch effects

M_corrected <- removeBatchEffect(
  M,
  batch  = batch_factors$DNA_batch,
  batch2 = batch_factors$Array_batch,
  covariates = batch_factors$FFPE_age,
  design = design
)

#convert back to beta values
beta.discovery.corrected = m_to_beta(M_corrected)
dim(beta.discovery.corrected)


#######re-run PCA and PERMANOVA
df2 = beta.discovery.corrected

#perform principal component analysis
pca_df2 = prcomp(t(df2), scale. = T)
pca_2_final = cbind(meta_pca, pca_df2$x[,1:10])
str(pca_2_final)

pca_2_final$DNA_batch = as.factor(pca_2_final$DNA_batch)
pca_2_final$Array_batch = as.factor(pca_2_final$Array_batch)
pca_2_final$Tumor_content = as.factor(pca_2_final$Tumor_content)

#####use permutational test for homogeneity of multivariate dispersions (PERMDISP) again on corrected data
##if <0.05, unequal dispersions, permanova significance might be inflated, else permanova reliable
pc_2_mat <- pca_2_final[, paste0("PC", 1:10)]

# Euclidean distances in PC space
dist_2_mat <- dist(pc_2_mat, method = "euclidean")

#PERMANOVA test

# DNA batch
adonis2(dist_2_mat ~ DNA_batch, data = pca_2_final, permutations = 999, method = "euclidean")

# Array batch
adonis2(dist_2_mat ~ Array_batch, data = pca_2_final, permutations = 999, method = "euclidean")

# Tumor content
adonis2(dist_2_mat ~ Tumor_content, data = pca_2_final, permutations = 999, method = "euclidean")

# FFPE age
adonis2(dist_2_mat ~ FFPE_age, data = pca_2_final, permutations = 999, method = "euclidean")

#results: none of the batches has significant effect on multivariate structure!!!

save.image()



###Scree plot for percent variances for corrected data

tiff("Scree_all_probes_correction.tiff", res=300,width=2560,height=2048)
fviz_screeplot(pca_df2, addlabels = T)
dev.off()

tiff("Scree_all_probes_correction_clean.tiff", res=300,width=2560,height=2048)
fviz_screeplot(pca_df2)
dev.off()


##visualize significant factors in PCA plots before correction
#DNA_batch
my_colors <- polychrome()
# Subset to the number of levels for DNA
my_colors_DNA <- polychrome()[1:length(levels(pca_final$DNA_batch))]
length(my_colors_DNA)
names(my_colors_DNA) <- levels(pca_final$DNA_batch)

tiff("PCA_DNA_batch_after_correction_legend.tiff", res=300,width=2460,height=1648)
ggplot(pca_2_final, aes(x=PC1, y=PC2, color = DNA_batch))+
  geom_jitter(height = 0.1, width = 0.1, size=3,show.legend = FALSE) +
  scale_color_manual(values=my_colors_DNA)+
  theme_classic()
dev.off()


#Array_batch
# Subset to the number of levels for DNA
my_colors_Array <- polychrome()[1:length(levels(pca_final$Array_batch))]
length(my_colors_Array)
names(my_colors_Array) <- levels(pca_final$Array_batch)

tiff("PCA_Array_batch_after_correction_legend.tiff", res=300,width=2460,height=1648)
ggplot(pca_2_final, aes(x=PC1, y=PC2, color = Array_batch))+
  geom_jitter(height = 0.1, width = 0.1, size=3,show.legend = FALSE) +
  scale_color_manual(values=my_colors_Array)+
  theme_classic()
dev.off()


#FFPE_Age
tiff("PCA_Age_batch_after_correction_legend.tiff", res=300,width=2460,height=1648)
ggplot(pca_2_final, aes(x=PC1, y=PC2, color = FFPE_age))+
  geom_jitter(height = 0.1, width = 0.1, size=3,show.legend = FALSE) +
  theme_classic()
dev.off()


# save corrected beta values to .RDS file for later use
saveRDS(beta.discovery.corrected, file = "Betas_discovery_preprocessed_corrected.RDS")

save.image()



