devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")

library(conumee2)
library(minfi)
library(ggplot2)
library(stringr)
library(lme4)
library(lmerTest) 
library(emmeans)
library(sesame)
library(data.table)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(here)

###clone figshare repository including the given folder structure to your local directory
###populate raw data subfolders with data from either GEO or figshare itself
###refer to the README.txt file at https://github.com/MerkLab/Meningioma_DNA_methylation_project for further information

#read in Pri_no_rec data
idat_dir = here("data", "GSE304094")
meta_dir = here("data", "general_datasets")
targets <- read.metharray.sheet(meta_dir, pattern="meta_discovery_current.csv")
setwd(here("data", "general_datasets"))
keep.samples =  read.csv(file = "meta_validation.csv", header = T)
keep.samples = keep.samples[keep.samples$setting_detail=="Primary_no_rec",]
keep.samples = keep.samples$Basename

RGSet_Pri_no_rec = read.metharray.exp(idat_dir, targets = targets)
RGSet_Pri_no_rec = RGSet_Pri_no_rec[, sampleNames(RGSet_Pri_no_rec) %in% keep.samples]
RGSet_Pri_no_rec@colData@rownames
sampleNames(RGSet_Pri_no_rec)
sampleNames(RGSet_Pri_no_rec) <- keep.samples
annotation(RGSet_Pri_no_rec) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSetPri_no_rec = preprocessNoob(RGSet_Pri_no_rec)

#change probe names to EPIC names, and reduce to probes on EPIC
MSetPri_no_rec@NAMES
MSetPrinorecnames = as.vector(MSetPri_no_rec@NAMES)
MSetPrinorecnamesnew = str_sub(MSetPrinorecnames, end = -6)
MSetPri_no_rec@NAMES = MSetPrinorecnamesnew
setwd(here("data", "10_longitudinal_CNV"))
EPICoverlap = read.csv(file = "EPICv2_EPIC_overlap.csv", header = F)
EPICoverlap = as.vector(EPICoverlap$V1)
EPICshareunique = EPICoverlap[ave(EPICoverlap, EPICoverlap, FUN = length) == 1]
MSet_Pri_no_rec_EPIC = subset(MSetPri_no_rec, rownames(MSetPri_no_rec) %in% EPICshareunique)
MSet_Pri_no_rec_EPIC
#dataset contains 718,114 probes identical in EPICv2 and EPIC

#read in control data
idat_dir_cont = here("data", "GSE109381")
cont = read.metharray.exp(base = idat_dir_cont)
cont
cont = preprocessNoob(cont)

#generate annotation shared for 450k and EPIC/EPICv2
data("exclude_regions")
anno = CNV.create_anno(array_type = "450k",bin_minprobes = 100,
                       chrXY = FALSE)

MSet_Pri_no_rec_EPIC = subset(MSet_Pri_no_rec_EPIC, rownames(MSet_Pri_no_rec_EPIC) %in% names(anno@probes))
MSet_Pri_no_rec_EPIC
anno@probes = subset(anno@probes, names(anno@probes) %in% rownames(MSet_Pri_no_rec_EPIC))
cont <- subset(cont, rownames(cont) %in% names(anno@probes))

#perform CNV analysis
t_Pri_no_rec = CNV.load(MSet_Pri_no_rec_EPIC)
c = CNV.load(cont)
save.image()

Pri_no_rec_fit = conumee2::CNV.fit(t_Pri_no_rec, c, anno)
Pri_no_rec_bin = CNV.bin(Pri_no_rec_fit)
Pri_no_rec_detail = CNV.detail(Pri_no_rec_bin)
Pri_no_rec_seg = CNV.segment(Pri_no_rec_detail)

setwd(here("results"))
tiff("Pri_no_rec_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Pri_no_rec_seg, threshold = 0.1)
dev.off()

tiff("Pri_no_rec_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Pri_no_rec_seg,zlim = c(-1, 1))
dev.off()


#read in Pri_rec data
idat_dir = here("data", "GSE304096")
meta_dir = here("data", "general_datasets")
targets = read.metharray.sheet(meta_dir, pattern="meta_validation.csv")
targets_Pri_rec = subset(targets, setting_detail == "Primary_rec")
RGSet_Pri_rec = read.metharray.exp(idat_dir, targets = targets)
keep.samples = targets_Pri_rec$Basename
RGSet_Pri_rec = RGSet_Pri_rec[, sampleNames(RGSet_Pri_rec) %in% keep.samples]
RGSet_Pri_rec@colData@rownames
sampleNames(RGSet_Pri_rec)
all(sampleNames(RGSet_Pri_rec) %in% targets_Pri_rec$Basename)
all(sampleNames(RGSet_Pri_rec) == targets_Pri_rec$Basename)
sampleNames(RGSet_Pri_rec) <- targets_Pri_rec$ID
annotation(RGSet_Pri_rec) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSetPri_rec = preprocessNoob(RGSet_Pri_rec)

#change probe names to EPIC names, and reduce to probes on EPIC
MSetPri_rec@NAMES
MSetPrirecnames = as.vector(MSetPri_rec@NAMES)
MSetPrirecnamesnew = str_sub(MSetPrirecnames, end = -6)
MSetPri_rec@NAMES = MSetPrirecnamesnew
MSet_Pri_rec_EPIC = subset(MSetPri_rec, rownames(MSetPri_rec) %in% EPICshareunique)
MSet_Pri_rec_EPIC
#dataset contains 718,114 probes identical in EPICv2 and EPIC

#generate annotation shared for 450k and EPIC/EPICv2

MSet_Pri_rec_EPIC = subset(MSet_Pri_rec_EPIC, rownames(MSet_Pri_rec_EPIC) %in% names(anno@probes))
MSet_Pri_rec_EPIC

#perform CNV analysis
t_Pri_rec = CNV.load(MSet_Pri_rec_EPIC)

Pri_rec_fit = CNV.fit(t_Pri_rec, c, anno)
Pri_rec_bin = CNV.bin(Pri_rec_fit)
Pri_rec_detail = CNV.detail(Pri_rec_bin)
Pri_rec_seg = CNV.segment(Pri_rec_detail)

setwd(here("results"))
tiff("Pri_rec_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Pri_rec_seg, threshold = 0.1)
dev.off()

tiff("Pri_rec_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Pri_rec_seg,zlim = c(-1, 1))
dev.off()


#read in Rec1 data
idat_dir = here("data", "GSE304096")
meta_dir = here("data", "general_datasets")
targets = read.metharray.sheet(meta_dir, pattern="meta_validation.csv")
targets_Rec1 = subset(targets, setting_detail == "Rec_1")
RGSet_Rec1 = read.metharray.exp(idat_dir, targets = targets)
keep.samples = targets_Rec1$Basename
RGSet_Rec1 = RGSet_Rec1[, sampleNames(RGSet_Rec1) %in% keep.samples]
RGSet_Rec1@colData@rownames
sampleNames(RGSet_Rec1)
all(sampleNames(RGSet_Rec1) %in% targets_Rec1$Basename)
all(sampleNames(RGSet_Rec1) == targets_Rec1$Basename)
sampleNames(RGSet_Rec1) <- targets_Rec1$ID
annotation(RGSet_Rec1) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSet_Rec1 = preprocessNoob(RGSet_Rec1)

#change probe names to EPIC names, and reduce to probes on EPIC
MSet_Rec1@NAMES
MSetRec1names = as.vector(MSet_Rec1@NAMES)
MSetRec1namesnew = str_sub(MSetRec1names, end = -6)
MSet_Rec1@NAMES = MSetRec1namesnew
MSet_Rec1_EPIC = subset(MSet_Rec1, rownames(MSet_Rec1) %in% EPICshareunique)
MSet_Rec1_EPIC
#dataset contains 718,114 probes identical in EPICv2 and EPIC

#generate annotation shared for 450k and EPIC/EPICv2

MSet_Rec1_EPIC = subset(MSet_Rec1_EPIC, rownames(MSet_Rec1_EPIC) %in% names(anno@probes))
MSet_Rec1_EPIC

#perform CNV analysis
t_Rec1 = CNV.load(MSet_Rec1_EPIC)

Rec1_fit = CNV.fit(t_Rec1, c, anno)
Rec1_bin = CNV.bin(Rec1_fit)
Rec1_detail = CNV.detail(Rec1_bin)
Rec1_seg = CNV.segment(Rec1_detail)

setwd(here("results"))
tiff("Rec1_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Rec1_seg, threshold = 0.1)
dev.off()

tiff("Rec1_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Rec1_seg,zlim = c(-1, 1))
dev.off()



#read in Rec2 data
idat_dir = here("data", "GSE304096")
meta_dir = here("data", "general_datasets")
targets = read.metharray.sheet(meta_dir, pattern="meta_validation.csv")
targets_Rec2 = subset(targets, setting_detail == "Rec_2")
RGSet_Rec2 = read.metharray.exp(idat_dir, targets = targets)
keep.samples = targets_Rec2$Basename
RGSet_Rec2 = RGSet_Rec2[, sampleNames(RGSet_Rec2) %in% keep.samples]
RGSet_Rec2@colData@rownames
sampleNames(RGSet_Rec2)
all(sampleNames(RGSet_Rec2) %in% targets_Rec2$Basename)
all(sampleNames(RGSet_Rec2) == targets_Rec2$Basename)
sampleNames(RGSet_Rec2) <- targets_Rec2$ID
annotation(RGSet_Rec2) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSet_Rec2 = preprocessNoob(RGSet_Rec2)

#change probe names to EPIC names, and reduce to probes on EPIC
MSet_Rec2@NAMES
MSetRec2names = as.vector(MSet_Rec2@NAMES)
MSetRec2namesnew = str_sub(MSetRec2names, end = -6)
MSet_Rec2@NAMES = MSetRec2namesnew
MSet_Rec2_EPIC = subset(MSet_Rec2, rownames(MSet_Rec2) %in% EPICshareunique)
MSet_Rec2_EPIC
#dataset contains 718,114 probes identical in EPICv2 and EPIC

#generate annotation shared for 450k and EPIC/EPICv2

MSet_Rec2_EPIC = subset(MSet_Rec2_EPIC, rownames(MSet_Rec2_EPIC) %in% names(anno@probes))
MSet_Rec2_EPIC

#perform CNV analysis
t_Rec2 = CNV.load(MSet_Rec2_EPIC)

Rec2_fit = CNV.fit(t_Rec2, c, anno)
Rec2_bin = CNV.bin(Rec2_fit)
Rec2_detail = CNV.detail(Rec2_bin)
Rec2_seg = CNV.segment(Rec2_detail)

setwd(here("results"))
tiff("Rec2_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Rec2_seg, threshold = 0.1)
dev.off()

tiff("Rec2_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Rec2_seg,zlim = c(-1, 1))
dev.off()


#########make CNV plots and get CNV information
##get entire validation info
#read in files
idat_dir = here("data", "GSE304096")
RGSet = read.metharray.exp(base = idat_dir)
annotation(RGSet) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
MSet = preprocessNoob(RGSet)

#change probe names of EPICv2 to match EPIC, and subset MSet of EPICv2 to only include probes shared by EPICv2 and EPIC
MSetnames = as.vector(MSet@NAMES)
MSetnamesnew = str_sub(MSetnames, end = -6)
MSet@NAMES = MSetnamesnew

#get info on EPICv2 and EPIC overlap and subset MSets by overlap
setwd(here("data", "10_longitudinal_CNV"))
EPICoverlap = read.csv(file = "EPICv2_EPIC_overlap.csv", header = F)
EPICoverlap = as.vector(EPICoverlap$V1)
EPICshareunique = EPICoverlap[ave(EPICoverlap, EPICoverlap, FUN = length) == 1]
MSet = subset(MSet, rownames(MSet) %in% EPICshareunique)
MSet

#read in control reference files
idat_dir_cont = here("data", "GSE109381")
cont = read.metharray.exp(base = idat_dir_cont)
cont
cont = preprocessNoob(cont)

#annotation with "overlap" to compare 450K and EPIC
data("exclude_regions")
anno = CNV.create_anno(array_type = "overlap", bin_minprobes = 50,
                       chrXY = FALSE, exclude_regions=exclude_regions,detail_regions = "Detail_NF2_CDKN2AB.bed")

#reduced MSet and cont to probes on anno
MSet = subset(MSet, rownames(MSet) %in% names(anno@probes))
MSet
anno@probes = subset(anno@probes, names(anno@probes) %in% rownames(MSet))
cont = subset(cont, rownames(cont) %in% names(anno@probes))
cont

#make useful names
setwd(here("data", "general_datasets"))
targets <- read.metharray.sheet(idat_dir, pattern="targets_validation_with_predictions.csv")
all(colnames(MSet) %in% targets$Basename)
all(colnames(MSet) == targets$Basename)
colnames(MSet) = targets$ID
colnames(MSet)


#bring into format
c = CNV.load(cont)
t = CNV.load(MSet)

#make the plots and get text files for CNVs (important use conumee, not conumee2, does not work)
setwd(here("results"))
for(x in names(t)){
  tmp <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(t[x], c,anno))))
  tiff(paste0(x,".tiff"),compression="lzw",res=300,width=5120,height=2048);
  CNV.genomeplot(tmp, ylim = c(-1.1, 1.1), cols = c("blue4", "blue2", "lightgrey", "darkorange2", "darkorange4"));dev.off()
  write.table(CNV.write(tmp,what="segments"),sep="\t",quote=F,row.names=F,file=paste0(x,".tsv"))
}



#######get CNV data for all samples
#coerce to one dataframe to work with

# get locations of all cn files
files <- Sys.glob(file.path(here("results"),"*.tsv"))

# create function to read in and format data
a <- function(x){
  # read data and set column names
  data <- read.delim(x, header=TRUE)
  colnames(data) <- c("ID", "chromosome", "start", "end", "probes","bstat", "pavl", "segmean", "segmedian")
  
  # return the data
  return(data)
}

# run the anonymous function defined above
CNVData <- lapply(files, a)

# turn the list of data frames into a single data frame
CNVData <- do.call("rbind", CNVData)
setwd(here("results"))
write.csv(CNVData, file= "CNVData_per_sample.csv")

##calculate genome instability for each sample as percentage of genome affected by CNVs
######get CNVs with cutoff 0.3!
setwd(here("data", "10_longitudinal_CNV"))
CNV_sizes = read.csv(file = "All_samples_CNV_sizes.csv", header = T)
str(CNV_sizes)
CNV_sizes = aggregate(CNV_sizes[-1], by = list(CNV_sizes$ID), FUN=sum)
setwd(here("results"))
write.csv(CNV_sizes, file = "CNV_sizes_aggregate.csv")

#make genome instability plot
setwd(here("data", "10_longitudinal_CNV"))
GeIn = read.csv(file="Percentage_genome_disrupted.csv", header = T)

setwd(here("results"))
cairo_pdf(filename = "Boxplot_genome_instability_MNG_settings.pdf", width = 5, height = 8)
ggplot(GeIn, 
       aes(condition, Percentage))+
  geom_point(size=4,position = position_jitter(width = 0.15),alpha=0.95, aes(color = condition))+
  scale_color_manual(values = c("#1CFDB2", "darkgreen", "#FFCB1B", "#D95F02"))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  scale_y_continuous(transform = "sqrt")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#test statistics
# Fit mixed-effects model
model <- lmer(
  Percentage ~ condition + (1 | Patient_ID),
  data = GeIn
)

# Overall test for condition
anova(model)

#post-hoc pairwise comparisions
emm <- emmeans(model, ~ condition)
pairs(emm, adjust = "tukey")




##make plot for unfavorable CNVs in all conditions
#get data
setwd(here("data", "10_longitudinal_CNV"))
CNVunfavor_MNG_settings = read.csv("CNV_unfavor_condition.csv", header = T)
str(CNVunfavor_MNG_settings)

setwd(here("results"))
cairo_pdf(filename = "Boxplot_CNVs_unfavorable_validation.pdf", width = 4.9, height = 8)
ggplot(CNVunfavor_MNG_settings, 
       aes(condition, CNV_unfavor))+
  geom_point(size=4,position = position_jitter(width = 0.2, height = 0),alpha=0.95, aes(color = condition))+
  scale_color_manual(values = c("#1CFDB2", "darkgreen", "#FFCB1B", "#D95F02"))+
  geom_boxplot(outlier.shape = NA, alpha=0.5)+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

# Fit mixed-effects model
model <- lmer(
  CNV_unfavor ~ condition + (1 | Patient_ID),
  data = CNVunfavor_MNG_settings
)

# Overall test for condition
anova(model)

#post-hoc pairwise comparisions
emm <- emmeans(model, ~ condition)
pairs(emm, adjust = "tukey")



#######get bin mean and make correlation
#########make CNV plots and get CNV information
##get entire validation info
#read in files
idat_dir = here("data", "GSE304094")
RGSet = read.metharray.exp(base = idat_dir)
annotation(RGSet) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
MSet = preprocessNoob(RGSet)

#change probe names to EPIC names, and reduce to probes on EPIC
MSet@NAMES
MSetnames = as.vector(MSet@NAMES)
MSetnamesnew = str_sub(MSetnames, end = -6)
MSet@NAMES = MSetnamesnew
MSet_EPIC = subset(MSet, rownames(MSet) %in% EPICshareunique)
MSet_EPIC
#dataset contains 718,114 probes identical in EPICv2 and EPIC

#read in control data
idat_dir_cont = here("data", "GSE109381")
cont = read.metharray.exp(base = idat_dir_cont)
cont
cont = preprocessNoob(cont)

#generate new anno file for fixed bin size
anno_bin= CNV.create_anno(bin_minprobes = 15, bin_minsize = 300000,
                bin_maxsize = 300000, array_type = "450k", chrXY = FALSE,
                exclude_regions = NULL, detail_regions = NULL)

#reduce validation, Cont, an anno_bin to same eprobes
#only 354,842 probes in validaiton set after overlap with anno_bin
MSet_EPIC = subset(MSet_EPIC, rownames(MSet_EPIC) %in% names(anno_bin@probes))
MSet_EPIC
anno_bin@probes = subset(anno_bin@probes, names(anno_bin@probes) %in% rownames(MSet_EPIC))
cont = subset(cont, rownames(cont) %in% names(anno_bin@probes))
#now all validation set, cont set, and anno have to same probes (n=354,842)

#perform CNV analysis up to bins
c = CNV.load(cont)
Val_load = conumee2::CNV.load(MSet_EPIC)
Val_fit = conumee2::CNV.fit(Val_load, c, anno_bin)
Val_bin = conumee2::CNV.bin(Val_fit)

# turn the list of data frames into a single data frame
Val_bin_ratios <- do.call("rbind", Val_bin@bin$ratio)
Val_bin_ratios = Val_bin_ratios[,order(colnames(Val_bin_ratios))]
Val_bin_ratios = as.data.frame(Val_bin_ratios)
all(rownames(Val_bin_ratios) %in% targets$Basename)
all(rownames(Val_bin_ratios) == targets$Basename)
rownames(Val_bin_ratios) = targets$ID
Val_bin_ratios = Val_bin_ratios[,-c(5552:5559)]
Val_bin_ratios = t(Val_bin_ratios)

#make correlation matrix/heatmap
Val_CNV_cor = cor(Val_bin_ratios)
anno_Val_CNV = targets[,c(1,4)]
rownames(anno_Val_CNV) = anno_Val_CNV$ID
anno_Val_CNV = anno_Val_CNV[,-1]
anno_Val_CNV = as.data.frame(anno_Val_CNV)
rownames(anno_Val_CNV) = targets$ID
colnames(anno_Val_CNV) = "setting_detail"
str(anno_Val_CNV)

ann_colors = list(setting_detail = c(Primary_no_rec="#1CFDB2",
                                      Primary_rec="darkgreen",
                                      Rec_1 = "#FFCB1B",
                                      Rec_2 = "#D95F02"))

pal=colorRampPalette(c("white", "#27408B"))


setwd(here("results"))
cairo_pdf(filename = "Correlation_heatmap_MNG_settings.pdf", width = 10, height = 8)
pheatmap(Val_CNV_cor,
         annotation = anno_Val_CNV, annotation_colors = ann_colors, border_color = NA, fontsize = 8, color=pal(1000))
dev.off()


#get selected boxplots for that correlation heatmap
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}

Cor_mat_flat = flattenCorrMatrix(Val_CNV_cor)
setwd(here("results"))
write.csv(Cor_mat_flat,file="Cor_matrix_CNV_Val_flat.csv")

#get sample IDs for Pri_no_rec, Pri_rec, and Recurrences
setwd(here("data", "10_longitudinal_CNV"))
Pri_No_rec_IDs = read.csv(file = "Pri_no_rec_IDs.csv", header = T)
Pri_No_rec_IDs = Pri_No_rec_IDs$ID

Pri_rec_IDs = read.csv(file = "Pri_rec_IDs.csv", header = T)
Pri_rec_IDs = Pri_rec_IDs$ID

Recurrences_IDs = read.csv(file = "Recurrences_IDs.csv", header = T)
Recurrences_IDs = Recurrences_IDs$ID

#create all pairwise combinations for those
pairs_Pri_no_rec = as.data.frame(t(combn(Pri_No_rec_IDs,2)))
colnames(pairs_Pri_no_rec) = c("row", "column")

pairs_Pri_rec = as.data.frame(t(combn(Pri_rec_IDs,2)))
colnames(pairs_Pri_rec) = c("row", "column")

pairs_recurrences = as.data.frame(t(combn(Recurrences_IDs,2)))
colnames(pairs_recurrences) = c("row", "column")



#match those pairs to flattened matrix
pairwise_cor_Pri_no_rec = Cor_mat_flat %>%
  filter(
    paste(pmin(row, column), pmax(row, column)) %in%
      paste(pmin(pairs_Pri_no_rec$row, pairs_Pri_no_rec$column), pmax(pairs_Pri_no_rec$row, pairs_Pri_no_rec$column))
  )

pairwise_cor_Pri_rec = Cor_mat_flat %>%
  filter(
    paste(pmin(row, column), pmax(row, column)) %in%
      paste(pmin(pairs_Pri_rec$row, pairs_Pri_rec$column), pmax(pairs_Pri_rec$row, pairs_Pri_rec$column))
  )

pairwise_cor_recurrences = Cor_mat_flat %>%
  filter(
    paste(pmin(row, column), pmax(row, column)) %in%
      paste(pmin(pairs_recurrences$row, pairs_recurrences$column), pmax(pairs_recurrences$row, pairs_recurrences$column))
  )

#also get within patient correlations (which is only primaries and recurrences from the same patient)
sample_info = data.frame(
  Sample = targets_rec$ID,
  Patient = targets_rec$patient
)

within_patient_pairs = sample_info %>%
  group_by(Patient) %>%
  summarise(Pairs = list(as.data.frame(t(combn(Sample,2)))), .groups = "drop") %>%
  tidyr::unnest(Pairs)
colnames(within_patient_pairs) = c("Patient", "row", "column")

pairwise_cor_within_patient = Cor_mat_flat %>%
  filter(
    paste(pmin(row, column), pmax(row, column)) %in%
      paste(pmin(within_patient_pairs$row, within_patient_pairs$column), pmax(within_patient_pairs$row, within_patient_pairs$column))
  )

setwd(here("results"))
write.csv(pairwise_cor_Pri_no_rec, file="Pairwise_cor_Pri_no_rec.csv")
write.csv(pairwise_cor_Pri_rec, file="Pairwise_cor_Pri_rec.csv")
write.csv(pairwise_cor_recurrences, file="Pairwise_cor_Recurrences.csv")
write.csv(pairwise_cor_within_patient, file="Pairwise_cor_within_patient.csv")

#get data for boxplot
setwd(here("data", "10_longitudinal_CNV"))
intra_cor_all = read.csv(file="Selected_intra_correlations.csv", header = T)
intra_cor_all$condition = factor(intra_cor_all$condition, levels = c("intra_pri_no_rec","intra_pri_rec", "intra_rec", "intra_patient"))

pal(1000)[(0.560739502-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#8290BC"
pal(1000)[(0.312164753-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#BAC2DA"
pal(1000)[(0.325717396-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#B7BFD8"
pal(1000)[(0.854554462-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#405698"

setwd(here("results"))
cairo_pdf(filename = "Boxplot_CNV_correlations.pdf", width = 3, height = 8)
ggplot(intra_cor_all, aes(x=condition, y=correlation, fill=condition))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#8290BC","#BAC2DA","#B7BFD8","#405698"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()






#####make custom genome plots, using sesame pipeline
##get entire validation info
#read in files
idat_dir = here("data", "GSE304096")
sdfs.q <- openSesame(idat_dir, prep = "QCDPB", func = NULL)

#read in control reference files
idat_dir_cont = here("data", "GSE109381")
sdfs.c <- openSesame(idat_dir_cont, prep = "QCDPB", func = NULL)

# create CNV data object from list of combined intensitiy values (SeSAMe)

data.q <- CNV.load(do.call(cbind, lapply(sdfs.q, totalIntensities)))
data.c <- CNV.load(do.call(cbind, lapply(sdfs.c, totalIntensities)))
data.q
data.c

#annotation with "overlap" to compare 450K and EPIC, define hg38 genome for EPICv2
data("exclude_regions")
data(detail_regions.hg38)
detail_regions

#make annotation object
anno <- conumee2::CNV.create_anno(array_type = c("450k", "EPICv2"), genome= "hg38", bin_minprobes = 45,exclude_regions = exclude_regions, detail_regions = "detail_genes_4.bed")

#make useful names
meta_dir = here("data", "general_datasets")
targets <- read.metharray.sheet(meta_dir, pattern="targets_validation_with_predictions.csv")
all(colnames(data.q@intensity) %in% targets$Basename)
all(colnames(data.q@intensity) == targets$Basename)
colnames(data.q@intensity) = targets$ID
all(colnames(data.q@intensity) == targets$ID)
colnames(data.q@intensity)

#perform CNV analysis
data(package="conumee2")
data(cancer_genes_hg38)

x <- conumee2::CNV.fit(data.q, data.c, anno)
x <- conumee2::CNV.bin(x)
x <- conumee2::CNV.detail(x)
x <- conumee2::CNV.segment(x)


#####add custom genomeblot function with no detail text
###based on the CNV.genomplot function from conumee2, this generates a custom function that produces a genomeplot with indicated colors
###based on detail genes from anno file, this will highlight selected genes as red dots

.cumsum0 <- function(x, left = TRUE, right = FALSE, n = NULL) {
  xx <- c(0, cumsum(as.numeric(x)))
  if (!left) 
    xx <- xx[-1]
  if (!right) 
    xx <- head(xx, -1)
  names(xx) <- n
  xx
}

setGeneric("CNV.genomeplotcustom", function(object, ...) {
  standardGeneric("CNV.genomeplotcustom")
})

#' @rdname CNV.genomeplotcustom
setMethod("CNV.genomeplotcustom", signature(object = "CNV.analysis"), function(object, chr = "all", centromere = TRUE, detail = TRUE,
                                                                         main = NULL, sig_cgenes = FALSE, nsig_cgenes = 3, output = "local", directory = getwd(), ylim = c(-1.25, 1.25),
                                                                         bins_cex = 0.75, set_par = TRUE,
                                                                         width = 12, height = 6, res = 720, cols = c("darkblue","darkblue", "lightgrey", "#F16729", "#F16729")){
  
  # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
  if (length(object@bin) == 0)
    stop("bin unavailable, run CNV.bin")
  # if(length(object@detail) == 0) stop('bin unavailable, run CNV.detail')
  if (length(object@seg) == 0)
    stop("bin unavailable, run CNV.seg")
  if (nrow(object@fit$ratio) < 300000) {
    centromere = FALSE
  }
  
  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
  }
  
  if (is.null(main)) {
    main = colnames(object@fit$ratio)
  }
  
  if (!is.null(main) & length(main) != ncol(object@fit$ratio)) {
    stop("please provide names for every sample")
  }
  
  for (i in 1:ncol(object@fit$ratio)) {
    
    message(main[i])
    
    if(output == "pdf"){
      p_names <- paste(directory,"/", main,"_genomeplot",".pdf",sep="")
      pdf(p_names[i], width = width, height = height)
      par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
    }
    
    if(output == "png"){
      p_names <- paste(directory,"/", main[i],"_genomeplot",".png",sep="")
      png(p_names[i], units = "in", width = width, height = height, res = res)
      par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
    }
    
    if (chr[1] == "all") {
      chr <- object@anno@genome$chr
    } else {
      chr <- intersect(chr, object@anno@genome$chr)
    }
    
    chr.cumsum0 <- .cumsum0(object@anno@genome[chr, "size"], n = chr)
    
    plot(NA, xlim = c(0, sum(as.numeric(object@anno@genome[chr, "size"])) -
                        0), ylim = ylim, xaxs = "i", xaxt = "n", yaxt = "n", xlab = NA,
         ylab = NA, main = main[i])
    abline(v = .cumsum0(object@anno@genome[chr, "size"], right = TRUE),
           col = "grey")
    if (centromere) {
      abline(v = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr,
                                                                                "pq"], col = "grey", lty = 2)
    }
    
    axis(1, at = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr,
                                                                                "size"]/2, labels = object@anno@genome[chr, "chr"], las = 2)
    if (all(ylim == c(-1.25, 1.25))) {
      axis(2, at = round(seq(-1.2, 1.2, 0.4), 1), las = 2)
    } else {
      axis(2, las = 2)
    }
    
    # ratio
    bin.ratio <- object@bin$ratio[[i]] - object@bin$shift[i]
    bin.ratio[bin.ratio < ylim[1]] <- ylim[1]
    bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
    
    p_size <- 1/object@bin$variance[[i]][names(object@anno@bins)]
    
    if(bins_cex == "standardized") {
      p_size[p_size <15] <- 0.2
      p_size[p_size >= 15 & p_size <22.5] <- 0.3
      p_size[p_size >= 22.5 & p_size <30] <- 0.4
      p_size[p_size >= 30 & p_size <37.5] <- 0.5
      p_size[p_size >= 37.5 & p_size <45] <- 0.6
      p_size[p_size >= 45 & p_size <52.5] <- 0.7
      p_size[p_size >= 52.5 & p_size <60] <- 0.8
      p_size[p_size > 60] <- 0.9
    } else if(bins_cex == "sample_level") {
      b <- boxplot.stats(p_size)
      outliers <- names(b$out)
      p_size[outliers] <- as.numeric(b$stats[5])
      p_size <- round(0.7*((p_size - min(p_size))/(max(p_size) - min(p_size)))+ 0.2, digits = 2) #scaling from 0.1:0.8 for cex using predefined bins to enable comparability between plots
    } else {
      p_size <- bins_cex
    }
    
    bin.ratio.cols <- apply(colorRamp(cols)((bin.ratio + max(abs(ylim)))/(2 *max(abs(ylim)))),
                            1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
    
    lines(chr.cumsum0[as.vector(seqnames(object@anno@bins))] + values(object@anno@bins)$midpoint,
          bin.ratio, type = "p", pch = 16, cex = p_size, col = bin.ratio.cols)
    
    
    for (l in seq(length(object@seg$summary[[i]]$seg.median))) {
      lines(c(object@seg$summary[[i]]$loc.start[l] + chr.cumsum0[object@seg$summary[[i]]$chrom[l]],
              object@seg$summary[[i]]$loc.end[l] + chr.cumsum0[object@seg$summary[[i]]$chrom[l]]),
            rep(min(ylim[2], max(ylim[1], object@seg$summary[[i]]$seg.median[l])),
                2) - object@bin$shift[i], col = "darkblue", lwd = 2)
    }
    
    # detail
    
    if (detail & length(object@detail) > 0 & ncol(object@anno@genome) == 2) {        #mouse array
      
      detail.ratio <- object@detail$ratio[[i]] - object@bin$shift[i]
      detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
      detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
      detail.ratio.above <- (detail.ratio > 0 & detail.ratio < 0.85) |
        detail.ratio < -0.85
      
      lines(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
            + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
            detail.ratio, type = "p", pch = 16, col = "red")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, detail.ratio, NA), labels = paste("  ", values(object@anno@detail)$name, sep = ""), adj = c(0,0.5),srt = 90, col = "#FFFFFF00")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, NA, detail.ratio), labels = paste(values(object@anno@detail)$name, "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "#FFFFFF00")
      
    } else if(ncol(object@anno@genome) != 2 & detail & length(object@detail) > 0) { #human array
      
      detail.ratio <- object@detail$ratio[[i]] - object@bin$shift[i]
      detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
      detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
      detail.ratio.above <- (detail.ratio > 0 & detail.ratio < 0.85) |
        detail.ratio < -0.85
      
      lines(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
            + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
            detail.ratio, type = "p", pch = 16, col = "red")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, detail.ratio, NA), labels = paste("  ", values(object@anno@detail)$name, sep = ""), adj = c(0,0.5),srt = 90, col = "#FFFFFF00")
      text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
           + chr.cumsum0[as.vector(seqnames(object@anno@detail))],
           ifelse(detail.ratio.above, NA, detail.ratio), labels = paste(values(object@anno@detail)$name, "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "#FFFFFF00")
      
      
      if(!length(object@detail$amp.detail.regions[[i]]) == 0 || !length(object@detail$del.detail.regions[[i]]) == 0){ #CNV focal was used
        
        sig.detail.regions.ratio <- c(object@detail$amp.detail.regions[[i]], object@detail$del.detail.regions[[i]])
        
        sig.detail.regions.ratio[sig.detail.regions.ratio < ylim[1]] <- ylim[1]
        sig.detail.regions.ratio[sig.detail.regions.ratio > ylim[2]] <- ylim[2]
        sig.detail.regions.ratio.above <- (sig.detail.regions.ratio > 0 & sig.detail.regions.ratio < 0.85) |
          sig.detail.regions.ratio < -0.85
        
        sig.detail.regions <- object@anno@detail[object@anno@detail$name %in% names(sig.detail.regions.ratio)]
        names(sig.detail.regions) <- sig.detail.regions$name
        sig.detail.regions.ratio <- sig.detail.regions.ratio[names(sig.detail.regions)]
        
        lines(start(sig.detail.regions) + (end(sig.detail.regions) - start(sig.detail.regions)) /2
              + chr.cumsum0[as.vector(seqnames(sig.detail.regions))],
              sig.detail.regions.ratio, type = "p", pch = 16, col = "red")
        text(start(sig.detail.regions) + (end(sig.detail.regions) - start(sig.detail.regions)) /2
             + chr.cumsum0[as.vector(seqnames(sig.detail.regions))],
             ifelse(sig.detail.regions.ratio.above, sig.detail.regions.ratio, NA), labels = paste("  ", names(sig.detail.regions), sep = ""), adj = c(0,0.5), srt = 90, col = "red")
        text(start(sig.detail.regions) + (end(sig.detail.regions) - start(sig.detail.regions)) /2
             + chr.cumsum0[as.vector(seqnames(sig.detail.regions))],
             ifelse(sig.detail.regions.ratio.above, NA, sig.detail.regions.ratio), labels = paste(names(sig.detail.regions), "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "red")
        
        if(sig_cgenes) {
          
          if(object@anno@args$genome == "hg38"){
            data("consensus_cancer_genes_hg38")
          } else {
            data("consensus_cancer_genes_hg19")
          }
          
          sig.cancer.genes.sorted <- names(sort(c(object@detail$amp.cancer.genes[[i]], abs(object@detail$del.cancer.genes[[i]])), decreasing = TRUE))
          
          if(nsig_cgenes > length(sig.cancer.genes.sorted)){
            nsig_cgenes <- length(sig.cancer.genes.sorted)
          }
          
          sig.cancer.genes.ratio <- c(object@detail$amp.cancer.genes[[i]], object@detail$del.cancer.genes[[i]])[sig.cancer.genes.sorted[1:nsig_cgenes]]
          
          sig.cancer.genes.ratio[sig.cancer.genes.ratio < ylim[1]] <- ylim[1]
          sig.cancer.genes.ratio[sig.cancer.genes.ratio > ylim[2]] <- ylim[2]
          sig.cancer.genes.ratio.above <- (sig.cancer.genes.ratio > 0 & sig.cancer.genes.ratio < 0.85) |
            sig.cancer.genes.ratio < -0.85
          
          sig.cancer.genes <- cancer_genes[names(sig.cancer.genes.ratio)]
          
          lines(start(sig.cancer.genes) + (end(sig.cancer.genes) - start(sig.cancer.genes)) /2
                + chr.cumsum0[as.vector(seqnames(sig.cancer.genes))],
                sig.cancer.genes.ratio, type = "p", pch = 16, col = "red")
          text(start(sig.cancer.genes) + (end(sig.cancer.genes) - start(sig.cancer.genes)) /2
               + chr.cumsum0[as.vector(seqnames(sig.cancer.genes))],
               ifelse(sig.cancer.genes.ratio.above, sig.cancer.genes.ratio, NA), labels = paste("  ", names(sig.cancer.genes), sep = ""), adj = c(0,0.5), srt = 90, col = "red")
          text(start(sig.cancer.genes) + (end(sig.cancer.genes) - start(sig.cancer.genes)) /2
               + chr.cumsum0[as.vector(seqnames(sig.cancer.genes))],
               ifelse(sig.cancer.genes.ratio.above, NA, sig.cancer.genes.ratio), labels = paste(names(sig.cancer.genes), "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "red")
          
        }}}
    
    if(is.element(output, c("pdf", "png"))){
      dev.off()
    }
  }
  
  if(is.element(output, c("pdf", "png"))){
    message(paste(ncol(object@fit$ratio)," files were created.", sep = ""))
  }
  
  if (set_par)
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
  
})


CNV.genomeplotcustom(x[25], sig_cgenes=FALSE, nsig_cgenes=3,bins_cex = 1.5)



#detect focal CNVs
data(package = "conumee2")
data(cancer_genes_hg38)
x = CNV.focal(x, conf=0.99)

name_list <- lapply(x@detail$amp.detail.regions, names)
max_len <- max(lengths(name_list))
padded <- lapply(name_list, function(x) {
  length(x) <- max_len
  return(x)
})
CancerAmpGenes <- as.data.frame(padded, stringsAsFactors = FALSE)

setwd(here("results"))
write.csv(CancerAmpGenes, file = "Cancer_AMP_genes.csv")


name_list <- lapply(x@detail$del.detail.regions, names)
max_len <- max(lengths(name_list))
padded <- lapply(name_list, function(x) {
  length(x) <- max_len
  return(x)
})
CancerDelGenes <- as.data.frame(padded, stringsAsFactors = FALSE)
write.csv(CancerDelGenes, file = "Cancer_DEL_genes.csv")

name_list <- lapply(focal@detail$del.detail.regions, names)
max_len <- max(lengths(name_list))
padded <- lapply(name_list, function(x) {
  length(x) <- max_len
  return(x)
})
DetailDelGenes <- as.data.frame(padded, stringsAsFactors = FALSE)
write.csv(DetailDelGenes, file = "Detail_DEL_genes.csv")

name_list <- lapply(focal@detail$amp.detail.regions, names)
max_len <- max(lengths(name_list))
padded <- lapply(name_list, function(x) {
  length(x) <- max_len
  return(x)
})
DetailAMPGenes <- as.data.frame(padded, stringsAsFactors = FALSE)
write.csv(DetailAMPGenes, file = "Detail_AMP_genes.csv")


####make barplots of depleted/amplified cancer gene census genes
#CancerDelGenes
setwd(here("data", "10_longitudinal_CNV"))
Plot_CancerDelGenes = read.csv(file="CancerDelGenes_barplot.csv", header = T)
str(Plot_CancerDelGenes)
Plot_CancerDelGenes$condition = factor(Plot_CancerDelGenes$condition, 
                                       levels = c("shared_rec1_rec2","shared_pri_rec2","shared_pri_rec1","unique_rec2","unique_rec1","unique_pri","common"))

setwd(here("results"))
cairo_pdf(filename = "Boxplot_Gene_CNVs_DEL_per_patient.pdf", width = 8, height = 5)
ggplot(Plot_CancerDelGenes, aes(fill=condition, y=value, x=patient)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("dodgerblue4","dodgerblue1", "deepskyblue", "mediumorchid4", "mediumorchid3", "violet","darkblue" ))+
  theme_classic()
dev.off()


tiff("CancerDelGenes_barplot.tiff", res=300,width=3060,height=2048)
ggplot(Plot_CancerDelGenes, aes(fill=condition, y=value, x=patient)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("dodgerblue4","dodgerblue1", "deepskyblue", "mediumorchid4", "mediumorchid3", "violet","darkblue" ))+
  theme_classic()
dev.off()


#CancerAmpGenes
setwd(here("data", "10_longitudinal_CNV"))
Plot_CancerAmpGenes = read.csv(file="CancerAmpGenes_barplot.csv", header = T)
str(Plot_CancerAmpGenes)
Plot_CancerAmpGenes$condition = factor(Plot_CancerAmpGenes$condition, 
                                       levels = c("shared_rec1_rec2","shared_pri_rec2","shared_pri_rec1","unique_rec2","unique_rec1","unique_pri","common"))

setwd(here("results"))
cairo_pdf(filename = "Boxplot_Gene_CNVs_AMP_per_patient.pdf", width = 8, height = 5)
ggplot(Plot_CancerAmpGenes, aes(fill=condition, y=value, x=patient)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("dodgerblue4","dodgerblue1", "deepskyblue", "mediumorchid4", "mediumorchid3", "violet","darkblue" ))+
  theme_classic()
dev.off()

tiff("CancerAmpGenes_barplot.tiff", res=300,width=3060,height=2048)
ggplot(Plot_CancerAmpGenes, aes(fill=condition, y=value, x=patient)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = c("dodgerblue4","dodgerblue1", "deepskyblue", "mediumorchid4", "mediumorchid3", "violet","darkblue" ))+
  theme_classic()
dev.off()


#check statistics on AMP and DEL genes across MNG settings (cohort size too small for meaningful statistics!)
#get validation cohort gene CNV matrix
setwd(here("data", "10_longitudinal_CNV"))
cnv_mat = as.matrix(read.csv(file="Gene_CNV_mat.csv"))
rownames(cnv_mat) = cnv_mat[,1]
cnv_mat = cnv_mat[,-1]

#get group data
groups <- c(rep("Pri_no_rec", 13), rep("Pri_rec", 18), rep("Recurrences", 24))
groups <- c(rep("Pri_no_rec", 13), rep("Recurrence-associated", 42))
groups <- c(rep("Primaries", 31), rep("Recurrences", 24))
group <- factor(groups)
group

#test_enrichment function
test_enrichment <- function(cnv_mat, group,
                            method = c("fisher", "chisq")) {
  method <- match.arg(method)
  stopifnot(ncol(cnv_mat) == length(group))
  
  p_fun <- switch(method,
                  fisher = function(x, g) fisher.test(table(g, x))$p.value,
                  chisq  = function(x, g) chisq.test(table(g, x),
                                                     correct = FALSE)$p.value)
  
  pvals <- apply(cnv_mat, 1, p_fun, g = group)
  
  tibble::tibble(
    event = rownames(cnv_mat),
    pval  = pvals,
    padj  = p.adjust(pvals, method = "BH")
  ) |>
    dplyr::arrange(padj)
}

# Run the test
res <- test_enrichment(cnv_mat, group, method = "fisher")

# View results
print(res)
setwd(here("results"))
write.csv(res, file="Gene_CNV_enrichment_MNG_settings.csv")



#######make oncoprint for selected CNV genes
#read in pri_no_rec, pri-rec, and rec individually and get optimized order
setwd(here("data", "10_longitudinal_CNV"))
mat_pri_no_rec = read.table(file = "CNVs_validated_pri_no_rec.txt", header = TRUE, stringsAsFactors = FALSE, 
                            sep = "\t",na.strings = character(0),fill = TRUE,colClasses = "character")
rownames(mat_pri_no_rec) = mat_pri_no_rec[, 1]
mat_pri_no_rec = mat_pri_no_rec[, -1]
mat_pri_no_rec = as.matrix(mat_pri_no_rec)

setwd(here("data", "10_longitudinal_CNV"))
mat_pri_rec = read.table(file = "CNVs_validated_pri_rec.txt", header = TRUE, stringsAsFactors = FALSE, 
                            sep = "\t",na.strings = character(0),fill = TRUE,colClasses = "character")
rownames(mat_pri_rec) = mat_pri_rec[, 1]
mat_pri_rec = mat_pri_rec[, -1]
mat_pri_rec = as.matrix(mat_pri_rec)

setwd(here("data", "10_longitudinal_CNV"))
mat_rec = read.table(file = "CNVs_validated_rec.txt", header = TRUE, stringsAsFactors = FALSE, 
                         sep = "\t",na.strings = character(0),fill = TRUE,colClasses = "character")
rownames(mat_rec) = mat_rec[, 1]
mat_rec = mat_rec[, -1]
mat_rec = as.matrix(mat_rec)


col_disc = c("DEL" = "blue4", "AMP" = "darkorange2")

alter_fun_disc = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  DEL = alter_graphic("rect", fill = col_disc["DEL"]),
  AMP = alter_graphic("rect", fill = col_disc["AMP"])
)

heatmap_legend_param_disc = list(title = "CNVs", at = c("DEL", "AMP"), 
                                 labels = c("Deletion", "Amplification"))


x=oncoPrint(mat_pri_no_rec,
          alter_fun = alter_fun_disc, col = col_disc, 
          heatmap_legend_param = heatmap_legend_param_disc, 
          alter_fun_is_vectorized = FALSE)

sampleorder_pri_no_rec = x@column_order
sampleorder_pri_no_rec
setwd(here("results"))
write.csv(sampleorder_pri_no_rec, file ="sampleorder_pri_no_rec.csv")

x=oncoPrint(mat_pri_rec,
            alter_fun = alter_fun_disc, col = col_disc, 
            heatmap_legend_param = heatmap_legend_param_disc, 
            alter_fun_is_vectorized = FALSE)

sampleorder_pri_rec = x@column_order
sampleorder_pri_rec
setwd(here("results"))
write.csv(sampleorder_pri_rec, file ="sampleorder_pri_rec.csv")

x=oncoPrint(mat_rec,
            alter_fun = alter_fun_disc, col = col_disc, 
            heatmap_legend_param = heatmap_legend_param_disc, 
            alter_fun_is_vectorized = FALSE)

sampleorder_rec = x@column_order
sampleorder_rec
setwd(here("results"))
write.csv(sampleorder_rec, file ="sampleorder_rec.csv")


#re-order full matrix
setwd(here("data", "10_longitudinal_CNV"))
sampleorder_all = read.csv(file ="sampleorder_all.csv", header = T)
sampleorder_all = sampleorder_all$order

setwd(here("data", "10_longitudinal_CNV"))
mat_all = read.table(file = "CNVs_validated_extended_all.txt", header = TRUE, stringsAsFactors = FALSE, 
                     sep = "\t",na.strings = character(0),fill = TRUE,colClasses = "character")
rownames(mat_all) = mat_all[, 1]
mat_all = mat_all[, -1]
mat_all = as.matrix(mat_all)
mat_all_order = mat_all[,sampleorder_all]


#get oncoprint

row_order = c("CCND1", "MDM2", "GATA3", "PIM1", "VHL", "CDKN2A/B", "NF2")

setwd(here("data", "10_longitudinal_CNV"))
oncoprint_anno = read.csv(file ="oncoprint_anno.csv", header = T)
oncoprint_anno=oncoprint_anno[,-1]

top_anno =HeatmapAnnotation(df = oncoprint_anno,
                            col = list(setting = c(Primary_no_rec="#1CFDB2",
                                                          Primary_rec="darkgreen",
                                                          Rec_1 = "#FFCB1B",
                                                          Rec_2 = "#D95F02"),
                                       grade = c("1" = "#9EBCDA", "2" = "#8C6BB1", "3" = "#810F7C"),
                                       group = c("Merlinintact" = "royalblue2","Immuneenriched" = "red3",
                                                 "hypermetabolic" = "forestgreen","proliferative" = "darkorange2")))


setwd(here("results"))
tiff("Oncoprint_full.tiff", res=300,width=5120,height=2048)
oncoPrint(mat_all_order,alter_fun = alter_fun_disc, col = col_disc,alter_fun_is_vectorized = FALSE,
          column_order = colnames(mat_all_order),row_order = row_order,
          heatmap_legend_param = heatmap_legend_param_disc,top_annotation = top_anno)
dev.off()
