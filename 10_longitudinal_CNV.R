library(conumee2)
library(minfi)
library(ggplot2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("conumee")
library(conumee)

remove.packages("conumee2")

devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")

library(conumee2)


#read in Pri_no_rec data
idat_Pri_no_rec = "/Users/lab/Desktop/Meningioma/data/T26_validation_Pri_no_rec"
targets_Pri_no_rec <- read.metharray.sheet(idat_Pri_no_rec, pattern="targets_Pri_no_rec.csv")

RGSet_Pri_no_rec = read.metharray.exp(idat_Pri_no_rec, targets = targets_Pri_no_rec)
RGSet_Pri_no_rec
RGSet_Pri_no_rec@colData@rownames
sampleNames(RGSet_Pri_no_rec)
all(sampleNames(RGSet_Pri_no_rec) %in% targets_Pri_no_rec$Basename)
all(sampleNames(RGSet_Pri_no_rec) == targets_Pri_no_rec$Basename)
sampleNames(RGSet_Pri_no_rec) <- targets_Pri_no_rec$ID
annotation(RGSet_Pri_no_rec) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSetPri_no_rec = preprocessNoob(RGSet_Pri_no_rec)

#change probe names to EPIC names, and reduce to probes on EPIC
library(stringr)
MSetPri_no_rec@NAMES
MSetPrinorecnames = as.vector(MSetPri_no_rec@NAMES)
MSetPrinorecnamesnew = str_sub(MSetPrinorecnames, end = -6)
MSetPri_no_rec@NAMES = MSetPrinorecnamesnew
EPICoverlap = read.csv(file = "EPICv2_EPIC_overlap.csv", header = F)
EPICoverlap = as.vector(EPICoverlap$V1)
EPICshareunique = EPICoverlap[ave(EPICoverlap, EPICoverlap, FUN = length) == 1]
MSet_Pri_no_rec_EPIC = subset(MSetPri_no_rec, rownames(MSetPri_no_rec) %in% EPICshareunique)
MSet_Pri_no_rec_EPIC
#dataset contains 718,114 probes identical in EPICv2 and EPIC

#read in control data
cont = read.metharray.exp("/Users/lab/Desktop/Meningioma/data/Capper_cont/")
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

tiff("Pri_no_rec_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Pri_no_rec_seg, threshold = 0.1)
dev.off()

tiff("Pri_no_rec_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Pri_no_rec_seg,zlim = c(-1, 1))
dev.off()

save.image()




#read in Pri_rec data
idat_Pri_rec = "/Users/lab/Desktop/Meningioma/data/T26_validation_Pri_rec"
targets_Pri_rec <- read.metharray.sheet(idat_Pri_rec, pattern="targets_Pri_rec.csv")

RGSet_Pri_rec = read.metharray.exp(idat_Pri_rec, targets = targets_Pri_rec)
RGSet_Pri_rec
RGSet_Pri_rec@colData@rownames
sampleNames(RGSet_Pri_rec)
all(sampleNames(RGSet_Pri_rec) %in% targets_Pri_rec$Basename)
all(sampleNames(RGSet_Pri_rec) == targets_Pri_rec$Basename)
sampleNames(RGSet_Pri_rec) <- targets_Pri_rec$ID
annotation(RGSet_Pri_rec) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSetPri_rec = preprocessNoob(RGSet_Pri_rec)

#change probe names to EPIC names, and reduce to probes on EPIC
library(stringr)
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

tiff("Pri_rec_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Pri_rec_seg, threshold = 0.1)
dev.off()

tiff("Pri_rec_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Pri_rec_seg,zlim = c(-1, 1))
dev.off()

save.image()


#read in Rec1 data
idat_Rec1 = "/Users/lab/Desktop/Meningioma/data/T26_validation_Rec1"
targets_Rec1 <- read.metharray.sheet(idat_Rec1, pattern="targets_Rec1.csv")

RGSet_Rec1 = read.metharray.exp(idat_Rec1, targets = targets_Rec1)
RGSet_Rec1
RGSet_Rec1@colData@rownames
sampleNames(RGSet_Rec1)
all(sampleNames(RGSet_Rec1) %in% targets_Rec1$Basename)
all(sampleNames(RGSet_Rec1) == targets_Rec1$Basename)
sampleNames(RGSet_Rec1) <- targets_Rec1$ID
annotation(RGSet_Rec1) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSet_Rec1 = preprocessNoob(RGSet_Rec1)

#change probe names to EPIC names, and reduce to probes on EPIC
library(stringr)
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

tiff("Rec1_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Rec1_seg, threshold = 0.1)
dev.off()

tiff("Rec1_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Rec1_seg,zlim = c(-1, 1))
dev.off()



#read in Rec2 data
idat_Rec2 = "/Users/lab/Desktop/Meningioma/data/T26_validation_Rec2"
targets_Rec2 <- read.metharray.sheet(idat_Rec2, pattern="targets_Rec2.csv")

RGSet_Rec2 = read.metharray.exp(idat_Rec2, targets = targets_Rec2)
RGSet_Rec2
RGSet_Rec2@colData@rownames
sampleNames(RGSet_Rec2)
all(sampleNames(RGSet_Rec2) %in% targets_Rec2$Basename)
all(sampleNames(RGSet_Rec2) == targets_Rec2$Basename)
sampleNames(RGSet_Rec2) <- targets_Rec2$ID
annotation(RGSet_Rec2) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

#normalize
MSet_Rec2 = preprocessNoob(RGSet_Rec2)

#change probe names to EPIC names, and reduce to probes on EPIC
library(stringr)
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

tiff("Rec2_CNV_summaryplot.tiff", res=300,width=5120,height=2048)
CNV.summaryplot(Rec2_seg, threshold = 0.1)
dev.off()

tiff("Rec2_CNV_heatmap.tiff", res=300, width=5120,height=4000)
CNV.heatmap(Rec2_seg,zlim = c(-1, 1))
dev.off()

save.image()



#########make CNV plots and get CNV information
##get entire validation info
#read in files
getwd()
idat_dir = "/Users/lab/Desktop/Meningioma/data/T26_validation"
RGSet = read.metharray.exp(base = idat_dir)
annotation(RGSet) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
MSet = preprocessNoob(RGSet)

#change probe names of EPICv2 to match EPIC, and subset MSet of EPICv2 to only include probes shared by EPICv2 and EPIC
library(stringr)
MSetnames = as.vector(MSet@NAMES)
MSetnamesnew = str_sub(MSetnames, end = -6)
MSet@NAMES = MSetnamesnew

#get info on EPICv2 and EPIC overlap and subset MSets by overlap
EPICoverlap = read.csv(file = "EPICv2_EPIC_overlap.csv", header = F)
EPICoverlap = as.vector(EPICoverlap$V1)
EPICshareunique = EPICoverlap[ave(EPICoverlap, EPICoverlap, FUN = length) == 1]
MSet = subset(MSet, rownames(MSet) %in% EPICshareunique)
MSet

#read in control reference files
idat_dir_cont = "/Users/lab/Desktop/Meningioma/data/Capper_cont"
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
targets <- read.metharray.sheet(idat_dir, pattern="targets_validation_with_predictions.csv")
all(colnames(MSet) %in% targets$Basename)
all(colnames(MSet) == targets$Basename)
colnames(MSet) = targets$ID
colnames(MSet)

#get custom function

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

setMethod("CNV.genomeplotcustom", signature(object = "CNV.analysis"), function(object, 
                                                                               chr = "all", chrX = TRUE, chrY = TRUE, centromere = TRUE, detail = TRUE, 
                                                                               main = NULL, ylim = c(-1.25, 1.25), set_par = TRUE, cols = c("red", 
                                                                                                                                            "red", "lightgrey", "green", "green")) {
  # if(length(object@fit) == 0) stop('fit unavailable, run CNV.fit')
  if (length(object@bin) == 0) 
    stop("bin unavailable, run CNV.bin")
  # if(length(object@detail) == 0) stop('bin unavailable, run
  # CNV.detail')
  if (length(object@seg) == 0) 
    stop("bin unavailable, run CNV.seg")
  
  if (set_par) {
    mfrow_original <- par()$mfrow
    mar_original <- par()$mar
    oma_original <- par()$oma
    par(mfrow = c(1, 1), mar = c(4, 4, 4, 4), oma = c(0, 0, 0, 0))
  }
  
  if (is.null(main)) 
    main <- object@name
  if (chr[1] == "all") {
    chr <- object@anno@genome$chr
  } else {
    chr <- intersect(chr, object@anno@genome$chr)
  }
  chr.cumsum0 <- .cumsum0(object@anno@genome[chr, "size"], n = chr)
  if (!chrX & is.element("chrX", names(chr.cumsum0))) 
    chr.cumsum0["chrX"] <- NA
  if (!chrY & is.element("chrY", names(chr.cumsum0))) 
    chr.cumsum0["chrY"] <- NA
  
  plot(NA, xlim = c(0, sum(as.numeric(object@anno@genome[chr, "size"])) - 
                      0), ylim = ylim, xaxs = "i", xaxt = "n", yaxt = "n", xlab = NA, 
       ylab = NA, main = main)
  abline(v = .cumsum0(object@anno@genome[chr, "size"], right = TRUE), 
         col = "grey")
  if (centromere) 
    abline(v = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr, 
                                                                              "pq"], col = "grey", lty = 2)
  axis(1, at = .cumsum0(object@anno@genome[chr, "size"]) + object@anno@genome[chr, 
                                                                              "size"]/2, labels = object@anno@genome[chr, "chr"], las = 2)
  if (all(ylim == c(-1.25, 1.25))) {
    axis(2, at = round(seq(-1.2, 1.2, 0.4), 1), las = 2)
  } else {
    axis(2, las = 2)
  }
  
  # ratio
  bin.ratio <- object@bin$ratio - object@bin$shift
  bin.ratio[bin.ratio < ylim[1]] <- ylim[1]
  bin.ratio[bin.ratio > ylim[2]] <- ylim[2]
  bin.ratio.cols <- apply(colorRamp(cols)((bin.ratio + max(abs(ylim)))/(2 * 
                                                                          max(abs(ylim)))), 1, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  
  lines(chr.cumsum0[as.vector(seqnames(object@anno@bins))] + values(object@anno@bins)$midpoint, 
        bin.ratio, type = "p", pch = 16, cex = 0.75, col = bin.ratio.cols)
  
  for (i in seq(length(object@seg$summary$seg.median))) {
    lines(c(object@seg$summary$loc.start[i] + chr.cumsum0[object@seg$summary$chrom[i]], 
            object@seg$summary$loc.end[i] + chr.cumsum0[object@seg$summary$chrom[i]]), 
          rep(min(ylim[2], max(ylim[1], object@seg$summary$seg.median[i])), 
              2) - object@bin$shift, col = "darkblue", lwd = 2)
  }
  
  # detail
  if (detail & length(object@detail) > 0) {
    detail.ratio <- object@detail$ratio - object@bin$shift
    detail.ratio[detail.ratio < ylim[1]] <- ylim[1]
    detail.ratio[detail.ratio > ylim[2]] <- ylim[2]
    detail.ratio.above <- (detail.ratio > 0 & detail.ratio < 0.85) | 
      detail.ratio < -0.85
    
    lines(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
          + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
          detail.ratio, type = "p", pch = 19, col = "red")
    text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
         + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
         ifelse(detail.ratio.above, detail.ratio, NA), labels = paste("  ", 
                                                                      values(object@anno@detail)$name, sep = ""), adj = c(0, 
                                                                                                                          0.5), srt = 90, col = "#FFFFFF00")
    text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
         + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
         ifelse(detail.ratio.above, NA, detail.ratio), labels = paste(values(object@anno@detail)$name, 
                                                                      "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "#FFFFFF00")
  }
  
  if (set_par) 
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})

#bring into format
c = CNV.load(cont)
t = CNV.load(MSet)

library(conumee)
#make the plots and get text files for CNVs (important use conumee, not conumee2, does not work)
for(x in names(t)){
  tmp <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(t[x], c,anno))))
  tiff(paste0(x,".tiff"),compression="lzw",res=300,width=5120,height=2048);
  CNV.genomeplotcustom(tmp, ylim = c(-1.1, 1.1), cols = c("blue4", "blue2", "lightgrey", "darkorange2", "darkorange4"));dev.off()
  write.table(CNV.write(tmp,what="segments"),sep="\t",quote=F,row.names=F,file=paste0(x,".tsv"))
}

#detail

for(x in names(t)){
  tmp <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(t[x], c,anno))))
  tiff(paste0(x,".tiff"),compression="lzw",res=300,width=5120,height=2048);
  CNV.genomeplotcustom(tmp, ylim = c(-1.1, 1.1), cols = c("blue4", "blue2", "lightgrey", "darkorange2", "darkorange4"));dev.off()
}


#######make new plots for 
##get entire validation info
#try with old conumee package


#this is a try with the new package
library(sesame)
#read in files
idat_dir = "/Users/lab/Desktop/Meningioma/data/T26_validation"
sdfs.q <- openSesame(idat_dir, prep = "QCDPB", func = NULL)

#read in control reference files
idat_dir_cont = "/Users/lab/Desktop/Meningioma/data/Capper_cont"
sdfs.c <- openSesame(idat_dir_cont, prep = "QCDPB", func = NULL)

# create CNV data object from list of combined intensitiy values (SeSAMe)

data.q <- CNV.load(do.call(cbind, lapply(sdfs.q, totalIntensities)))
data.c <- CNV.load(do.call(cbind, lapply(sdfs.c, totalIntensities)))
data.q
data.c

#make useful names
targets <- read.metharray.sheet(idat_dir, pattern="targets_validation_with_predictions.csv")
all(colnames(data.q@intensity) %in% targets$Basename)
all(colnames(data.q@intensity) == targets$Basename)
colnames(data.q@intensity) = targets$ID
all(colnames(data.q@intensity) == targets$ID)
colnames(data.q@intensity)

#new annotation
anno <- CNV.create_anno(array_type = c("450k", "EPICv2"), genome= "hg38", exclude_regions = exclude_regions, detail_regions = "Detail_gene_CNVs.bed")

x <- CNV.fit(data.q, data.c, anno)
x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)




#######get CNV data for all samples
#coerce to one dataframe to work with
library("stringr")

# get locations of all cn files
files <- Sys.glob("/Users/lab/Desktop/Meningioma/T26_validation_CNV/results/CNV_plots_txt/*.tsv")

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
write.csv(CNVData, file= "CNVData_per_sample.csv")

##calculate genome instability for each sample as percentage of genome affected by CNVs
######get CNVs with cutoff 0.3!
CNV_sizes = read.csv(file = "All_samples_CNV_sizes.csv", header = T)
str(CNV_sizes)
CNV_sizes = aggregate(CNV_sizes[-1], by = list(CNV_sizes$ID), FUN=sum)
write.csv(CNV_sizes, file = "CNV_sizes_aggregate.csv")

#make genome instability plot
GeIn = read.csv(file="Percentage_genome_disrupted.csv", header = T)

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
library(lme4)
library(lmerTest) 
library(emmeans)

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
CNVunfavor_MNG_settings = read.csv("CNV_unfavor_condition.csv", header = T)
str(CNVunfavor_MNG_settings)

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
getwd()
idat_dir = "/Users/lab/Desktop/Meningioma/data/T26_validation"
RGSet = read.metharray.exp(base = idat_dir)
annotation(RGSet) = c(array= "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
MSet = preprocessNoob(RGSet)

#change probe names to EPIC names, and reduce to probes on EPIC
library(stringr)
MSet@NAMES
MSetnames = as.vector(MSet@NAMES)
MSetnamesnew = str_sub(MSetnames, end = -6)
MSet@NAMES = MSetnamesnew
MSet_EPIC = subset(MSet, rownames(MSet) %in% EPICshareunique)
MSet_EPIC
#dataset contains 718,114 probes identical in EPICv2 and EPIC

#read in control data
cont = read.metharray.exp("/Users/lab/Desktop/Meningioma/data/Capper_cont/")
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
library(pheatmap)
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


cairo_pdf(filename = "Correlation_heatmap_MNG_settings.pdf", width = 10, height = 8)
pheatmap(Val_CNV_cor,
         annotation = anno_Val_CNV, annotation_colors = ann_colors, border_color = NA, fontsize = 8, color=pal(1000))
dev.off()

save.image()


#get selected boxplots for that correlation heatmap
library(data.table)
flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}

Cor_mat_flat = flattenCorrMatrix(Val_CNV_cor)
write.csv(Cor_mat_flat,file="Cor_matrix_CNV_Val_flat.csv")

#get sample IDs for Pri_no_rec, Pri_rec, and Recurrences
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
library(dplyr)
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

write.csv(pairwise_cor_Pri_no_rec, file="Pairwise_cor_Pri_no_rec.csv")
write.csv(pairwise_cor_Pri_rec, file="Pairwise_cor_Pri_rec.csv")
write.csv(pairwise_cor_recurrences, file="Pairwise_cor_Recurrences.csv")
write.csv(pairwise_cor_within_patient, file="Pairwise_cor_within_patient.csv")

#get data for boxplot
intra_cor_all = read.csv(file="Selected_intra_correlations.csv", header = T)
intra_cor_all$condition = factor(intra_cor_all$condition, levels = c("intra_pri_no_rec","intra_pri_rec", "intra_rec", "intra_patient"))


pal(1000)[(0.560739502-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#8290BC"
pal(1000)[(0.312164753-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#BAC2DA"
pal(1000)[(0.325717396-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#B7BFD8"
pal(1000)[(0.854554462-min(intra_cor_all$correlation))/((max(intra_cor_all$correlation)-min(intra_cor_all$correlation))/1000)] #"#405698"

cairo_pdf(filename = "Boxplot_CNV_correlations.pdf", width = 3, height = 8)
ggplot(intra_cor_all, aes(x=condition, y=correlation, fill=condition))+
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=c("#8290BC","#BAC2DA","#B7BFD8","#405698"))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()






#####check focal CNVs and which of them are cancer-associated genes, use sesame pipeline as minfi generates error at segment step
library(conumee2)
library(sesame)

##get entire validation info
#read in files
idat_dir = "/Users/lab/Desktop/Meningioma/data/T26_validation"
sdfs.q <- openSesame(idat_dir, prep = "QCDPB", func = NULL)

#read in control reference files
idat_dir_cont = "/Users/lab/Desktop/Meningioma/data/Capper_cont"
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

#try using CNV.define_detail for cancer census genes (should however be covered by sig_cgenes below, but somehow is not)

cancer_genes_custom = cancer_genes@elementMetadata@listData$SYMBOL
cancer_genes_custom = as.data.frame(cancer_genes_custom)
cancer_genes_custom <- paste0('"', cancer_genes_custom$cancer_genes_custom, '"', collapse = ", ")
cat(cancer_genes_custom)
detail_regions_custom_1 = CNV.define_detail(array_type = "EPICv2", gene = "predefined")

#does not seem to work, problem with seqnames which are not in seqinfo, means that by default anno function uses seqinfo other than what is given in seqnames
library(biomaRt)
genes <- c("SKI", "TNFRSF14", "PRDM16", "RPL22", "CAMTA1", "MTOR", "PRDM2", "CASP9", "SPEN", "SDHB",
           "ARHGEF10L", "PAX7", "ID3", "MDS2", "ARID1A", "LCK", "SFPQ", "THRAP3", "CSF3R", "MYCL", "MPL", "MUTYH",
           "PIK3R3", "TAL1", "STIL", "CDKN2C", "EPS15", "JUN", "JAK1", "FUBP1", "BCL10", "RPL5", "RBM15", "TRIM33",
           "NRAS", "ATP1A1", "TENT5C", "NOTCH2", "BCL9", "PDE4DIP", "ARNT", "SETDB1", "MLLT11", "S100A7", "TPM3",
           "MUC1", "LMNA", "PRCC", "NTRK1", "FCRL4", "SDHC", "FCGR2B", "DDR2", "PBX1", "PRRX1", "ABL2", "TPR",
           "CDC73", "ASPM", "PTPRC", "ELF3", "BTG2", "MDM4", "ELK4", "SLC45A3", "RGS7", "FH", "AKT3",
           "MYCN", "WDCP", "NCOA1", "DNMT3A", "ASXL2", "ALK", "BIRC6", "STRN", "EML4", "SIX2", "EPAS1", "MSH2",
           "MSH6", "FBXO11", "BCL11A", "REL", "XPO1", "PCBP1", "DCTN1", "CTNNA2", "TMEM127", "AFF3", "RGPD3",
           "RANBP2", "PAX8", "ERCC3", "CXCR4", "LRP1B", "ACVR2A", "ACVR1", "HOXD13", "HOXD11", "NFE2L2", "ITGAV",
           "COL3A1", "PMS1", "SF3B1", "CASP8", "CD28", "CREB1", "IDH1", "ERBB4", "BARD1", "ATIC", "FEV", "PAX3", "ACSL3",
           "CUL3", "ACKR3", "SRGAP3", "FANCD2", "VHL", "PPARG", "RAF1", "FBLN2", "XPC", "TGFBR2", "CCR4", "MLH1", "MYD88",
           "CTNNB1", "SETD2", "NCKIPSD", "RHOA", "BAP1", "PBRM1", "CACNA1D", "FHIT", "MITF", "FOXP1", "ROBO2", "EPHA3",
           "TFG", "CBLB", "GSK3B", "POLQ", "GATA2", "RPN1", "CNBP", "STAG1", "PIK3CB", "FOXL2", "ATR", "WWTR1", "GMPS",
           "MLF1", "MECOM", "TBL1XR1", "PIK3CA", "SOX2", "MAP3K13", "IGF2BP2", "ETV5", "EIF4A2", "BCL6", "LPP", "TP63",
           "MB21D2", "MUC4", "TFRC", "FGFR3", "NSD2", "SLC34A2", "N4BP2", "RHOH", "PHOX2B", "TEC", "FIP1L1", "CHIC2",
           "PDGFRA", "KIT", "KDR", "PTPN13", "AFF1", "RAP1GDS1", "TET2", "LEF1", "IL2", "FAT4", "FBXW7", "CASP3", "FAT1",
           "SDHA", "TERT", "CTNND2", "CDH10", "DROSHA", "GOLPH3", "SUB1", "IL7R", "LIFR", "IL6ST", "MAP3K1", "PIK3R1",
           "RAD17", "APC", "ACSL6", "RAD50", "AFF4", "CTNNA1", "ARHGAP26", "CSF1R", "PDGFRB", "CD74", "ITK", "EBF1",
           "PWWP2A", "TLX3", "NPM1", "FGFR4", "NSD1", "FLT4", "IRF4", "DEK", "TRIM27", "HLA-A",
           "POU5F1", "DAXX", "HMGA1", "FANCE", "SRSF3", "CDKN1A", "PIM1", "TFEB", "CCND3", "HSP90AB1", "NFKBIE", "BMP5",
           "EPHA7", "CCNC", "PRDM1", "FOXO3", "ROS1", "GOPC", "RSPO3", "PTPRK", "SGK1", "MYB", "BCLAF1", "TNFAIP3",
           "ECT2L", "LATS1", "ESR1", "ARID1B", "EZR", "QKI", "AFDN", "CARD11", "TNRC18", "PMS2", "RAC1",
           "ETV1", "MACC1", "HNRNPA2B1", "HOXA9", "HOXA11", "HOXA13", "JAZF1", "FKBP9", "SFRP4", "IKZF1", "EGFR", "ZNF479",
           "SBDS", "ELN", "HIP1", "HGF", "GRM3", "AKAP9", "CDK6", "TRRAP", "CUX1", "MET", "POT1", "SND1", "SMO", "CREB3L2",
           "TRIM24", "KIAA1549", "BRAF", "FAM131B", "CNTNAP2", "EZH2", "KMT2C", "MNX1", "ARHGEF10", "PCM1", "LEPROTL1",
           "WRN", "NRG1", "NSD3", "FGFR1", "ANK1", "KAT6A", "IKBKB", "HOOK3", "TCEA1", "LYN", "PLAG1", "CHCHD7", "PREX2",
           "NCOA2", "HEY1", "CNBD1", "NBN", "RUNX1T1", "CDH17", "COX6C", "PABPC1", "UBR5", "RSPO2", "EIF3E", "CSMD3",
           "RAD21", "EXT1", "MYC", "NDRG1", "FAM135B", "PLEC", "RECQL4", "JAK2", "CD274", "PDCD1LG2", "PTPRD", "NFIB",
           "PSIP1", "MLLT3", "CDKN2A", "FANCG", "PAX5", "GNAQ", "NTRK2", "SYK", "OMD", "WNK2", "FANCC", "PTCH1", "XPA",
           "NR4A3", "TAL2", "KLF4", "TNC", "CNTRL", "PPP6C", "SET", "FNBP1", "ABL1", "NUP214", "TSC1", "RALGDS", "BRD3",
           "RXRA", "NOTCH1", "LARP4B", "KLF6", "GATA3", "MLLT10", "ABI1", "ZEB1", "KIF5B", "RET", "NCOA4", "A1CF", "CCDC6",
           "TET1", "PRF1", "KAT6B", "BMPR1A", "NUTM2D", "PTEN", "FAS", "CPEB3", "CYP2C8", "TLX1", "NFKB2",
           "SUFU", "NT5C2", "VTI1A", "TCF7L2", "SHTN1", "FGFR2", "MGMT", "HRAS", "MUC6", "NUP98", "LMO1",
           "RRAS2", "MYOD1", "FANCF", "WT1", "LMO2", "EXT2", "CREB3L1", "DDB2", "CLP1", "CTNND1", "SDHAF2", "FEN1",
           "MEN1", "MALAT1", "CCND1", "FADD", "NUMA1", "PICALM", "EED", "FAT3", "MAML2", "BIRC3", "ATM", "DDX10",
           "POU2AF1", "SDHD", "ZBTB16", "PAFAH1B2", "KMT2A", "DDX6", "BCL9L", "FOXR1", "CBL", "ARHGEF12", "FLI1",
           "KCNJ5", "KDM5A", "ERC1", "CCND2", "CHD4", "ZNF384", "PTPN6", "ETV6", "CDKN1B", "ETNK1", "KRAS", "PPFIBP1",
           "ARID2", "COL2A1", "KMT2D", "PRPF40B", "SMARCD1", "ATF1", "ACVR1B", "HOXC13", "HOXC11", "ERBB3", "NACA",
           "NAB2", "STAT6", "GLI1", "DDIT3", "CDK4", "LRIG3", "WIF1", "HMGA2", "RAP1B", "MDM2", "PTPRB", "BTG1",
           "USP44", "CHST11", "SH2B3", "ALDH2", "PTPN11", "TBX3", "HNF1A", "SETD1B", "BCL7A", "CLIP1", "ZCCHC8", "NCOR2",
           "POLE", "ZMYM2", "LATS2", "CDX2", "FLT3", "BRCA2", "NBEA", "LHFPL6", "FOXO1", "LCP1", "RB1", "CYSLTR2", "GPC5",
           "SOX21", "ERCC5", "CCNB1IP1", "PRKD1", "ARHGAP5", "BAZ1A", "NKX2-1", "FOXA1", "NIN", "KTN1", "SIX1", "HIF1A",
           "MAX", "GPHN", "RAD51B", "TSHR", "TRIP11", "GOLGA5", "DICER1", "TCL1A", "BCL11B", "HSP90AA1", "AKT1", "NUTM1",
           "BUB1B", "KNSTRN", "KNL1", "B2M", "HMGN2P46", "USP8", "MYO5A", "C15orf65", "TCF12", "MAP2K1", "SMAD3", "PML",
           "NTRK3", "POLG", "IDH2", "CRTC3", "BLM", "FES", "CHD2", "NR2F2", "AXIN1", "NTHL1", "TSC2", "TRAF7", "CREBBP",
           "GRIN2A", "CIITA", "SOCS1", "RMI2", "TNFRSF17", "SNX29", "ERCC4", "MYH11", "PALB2", "PRKCB", "IL21R", "FUS",
           "CYLD", "HERPUD1", "CDH11", "CBFB", "CTCF", "CDH1", "ZFHX3", "RFWD3", "MAF", "CBFA2T3", "FANCA", "YWHAE",
           "USP6", "RABEP1", "POLR2A", "TP53", "PER1", "GAS7", "MAP2K4", "NCOR1", "FLCN", "SPECC1", "NF1", "SUZ12",
           "TAF15", "MLLT6", "LASP1", "CDK12", "ERBB2", "IKZF3", "RARA", "CCR7", "SMARCE1", "STAT5B", "STAT3", "BRCA1",
           "ETV4", "SPOP", "KAT7", "COL1A1", "HLF", "MSI2", "RNF43", "CLTC", "PPM1D", "BRIP1", "CD79B", "DDX5", "AXIN2",
           "PRKAR1A", "SRSF2", "CANT1", "RNF213", "ASPSCR1", "ZNF521", "SS18", "SETBP1", "SMAD2",
           "SMAD4", "DCC", "MALT1", "BCL2", "KDSR", "FSTL3", "STK11", "TCF3", "GNA11", "MAP2K2", "SH3GL1", "MLLT1",
           "VAV1", "CD209", "MUC16", "KEAP1", "DNM2", "SMARCA4", "CALR", "LYL1", "PRKACA", "DNAJB1", "BRD4", "TPM4",
           "JAK3", "ELL", "CRTC1", "ZNF429", "CCNE1", "CEP89", "CEBPA", "LSM14A", "AKT2", "CD79A", "CIC", "BCL3",
           "CBLC", "ERCC2", "ARHGAP35", "BAX", "BCL2L12", "POLD1", "KLK2", "PPP2R1A", "ZNF331", "TFPT", "CNOT3",
           "SIRPA", "CRNKL1", "ASXL1", "SRC", "MAFB", "TOP1", "PLCG1", "PTPRT", "SDC4", "NFATC2", "SALL4", "GNAS",
           "SS18L1", "PTK6", "OLIG2", "RUNX1", "ERG", "TMPRSS2", "U2AF1", "CLTCL1", "DGCR8", "LZTR1", "MAPK1",
           "BCR", "SMARCB1", "MN1", "CHEK2", "ZNRF3", "EWSR1", "NF2", "PATZ1", "ISX", "MYH9", "APOBEC3B", "PDGFB",
           "MRTFA", "EP300",  "TMSB4X", "ZRSR2", "EIF1AX", "FAM47C", "BCOR", "USP9X", "DDX3X",
           "KDM6A", "RBM10", "ARAF", "WAS", "GATA1", "TFE3", "KDM5C", "SMC1A", "AMER1",
           "MSN", "AR", "FOXO4", "ZMYM3", "NONO", "ATRX", "BTK", "IRS4", "STAG2", "DCAF12L2",
           "BCORL1", "ELF4", "GPC3", "PHF6", "ATP2B3", "FLNA", "RPL10", "MTCP1")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_coords <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
                     filters = "hgnc_symbol",
                     values = genes,
                     mart = mart)
gene_coords <- gene_coords[gene_coords$chromosome_name %in% c(1:22, "X", "Y"), ]

valid_coords <- gene_coords[
  gene_coords$chromosome_name %in% c(1:22) &
    !is.na(gene_coords$start_position), 
]

cancer_genes_custom <- paste0('"', valid_coords$hgnc_symbol, '"', collapse = ", ")
cat(cancer_genes_custom)

detail_regions_custom_3 = CNV.define_detail(array_type = "EPICv2", gene = c("A1CF", "ABI1", "ABL1", "ABL2", "ACKR3", "ACSL3", "ACSL6", "ACVR1", "ACVR1B", "ACVR2A", "AFDN", "AFF1", "AFF3", "AFF4",
                                                                            "AKAP9", "AKT1", "AKT2", "AKT3", "ALDH2", "ALK", "ANK1", "APC", "APOBEC3B", "ARHGAP26", "ARHGAP35", "ARHGAP5", "ARHGEF10",
                                                                            "ARHGEF10L", "ARHGEF12", "ARID1A", "ARID1B", "ARID2", "ARNT", "ASPM", "ASPSCR1", "ASXL1", "ASXL2", "ATF1", "ATIC", "ATM", "ATP1A1",
                                                                            "ATR", "AXIN1", "AXIN2", "B2M", "BAP1", "BARD1", "BAX", "BAZ1A", "BCL10", "BCL11A", "BCL11B", "BCL2", "BCL2L12", "BCL3", "BCL6",
                                                                            "BCL7A", "BCL9", "BCL9L", "BCLAF1", "BCR", "BIRC3", "BIRC6", "BLM", "BMP5", "BMPR1A", "BRAF", "BRCA1", "BRCA2", "BRD3", "BRD4",
                                                                            "BRIP1", "BTG1", "BTG2", "BUB1B", "CACNA1D", "CALR", "CAMTA1", "CANT1", "CARD11", "CASP3", "CASP8", "CASP9", "CBFA2T3", "CBFB", "CBL",
                                                                            "CBLB", "CBLC", "CCDC6", "CCNB1IP1", "CCNC", "CCND1", "CCND2", "CCND3", "CCNE1", "CCR4", "CCR7", "CD209", "CD274", "CD28", "CD74",
                                                                            "CD79A", "CD79B", "CDC73", "CDH1", "CDH10", "CDH11", "CDH17", "CDK12", "CDK4", "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2C", "CDX2",
                                                                            "CEBPA", "CEP89", "CHCHD7", "CHD2", "CHD4", "CHEK2", "CHIC2", "CHST11", "CIC", "CIITA", "CLIP1", "CLP1", "CLTC", "CLTCL1", "CNBD1",
                                                                            "CNBP", "CNOT3", "CNTNAP2", "CNTRL", "COL1A1", "COL2A1", "COL3A1", "COX6C", "CPEB3", "CREB1", "CREB3L1", "CREB3L2", "CREBBP", "CRNKL1", 
                                                                            "CRTC1", "CRTC3", "CSF1R", "CSF3R", "CSMD3", "CTCF", "CTNNA1", "CTNNA2", "CTNNB1", "CTNND1", "CTNND2", "CUL3", "CUX1", "CXCR4", "CYLD",
                                                                            "CYP2C8", "CYSLTR2", "DAXX", "DCC", "DCTN1", "DDB2", "DDIT3", "DDR2", "DDX10", "DDX5", "DDX6", "DEK", "DGCR8", "DICER1", "DNAJB1", "DNM2",
                                                                            "DNMT3A", "DROSHA", "EBF1", "ECT2L", "EED", "EGFR", "EIF3E", "EIF4A2", "ELF3", "ELK4", "ELL", "ELN", "EML4", "EP300", "EPAS1", "EPHA3",
                                                                            "EPHA7", "EPS15", "ERBB2", "ERBB3", "ERBB4", "ERC1", "ERCC2", "ERCC3", "ERCC4", "ERCC5", "ERG", "ESR1", "ETNK1", "ETV1", "ETV4", "ETV5",
                                                                            "ETV6", "EWSR1", "EXT1", "EXT2", "EZH2", "EZR", "FADD", "FAM131B", "FAM135B", "FANCA", "FANCC", "FANCD2", "FANCE", "FANCF", "FANCG", "FAS", 
                                                                            "FAT1", "FAT3", "FAT4", "FBLN2", "FBXO11", "FBXW7", "FCGR2B", "FCRL4", "FEN1", "FES", "FEV", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FH", "FHIT",
                                                                            "FIP1L1", "FKBP9", "FLCN", "FLI1", "FLT3", "FLT4", "FNBP1", "FOXA1", "FOXL2", "FOXO1", "FOXO3", "FOXP1", "FOXR1", "FSTL3", "FUBP1", "FUS", "GAS7",
                                                                            "GATA2", "GATA3", "GLI1", "GMPS", "GNA11", "GNAQ", "GNAS", "GOLGA5", "GOLPH3", "GOPC", "GPC5", "GPHN", "GRIN2A", "GRM3", "GSK3B", "HERPUD1",
                                                                            "HEY1", "HGF", "HIF1A", "HIP1", "HLA-A", "HLF", "HMGA1", "HMGA2", "HMGN2P46", "HNF1A", "HNRNPA2B1", "HOOK3", "HOXA11", "HOXA13", "HOXA9", 
                                                                            "HOXC11", "HOXC13", "HOXD11", "HOXD13", "HRAS", "HSP90AA1", "HSP90AB1", "ID3", "IDH1", "IDH2", "IGF2BP2", "IKBKB", "IKZF1", "IKZF3", "IL2",
                                                                            "IL21R", "IL6ST", "IL7R", "IRF4", "ISX", "ITGAV", "ITK", "JAK1", "JAK2", "JAK3", "JAZF1", "JUN", "KAT6A", "KAT6B", "KAT7", "KCNJ5", "KDM5A",
                                                                            "KDR", "KDSR", "KEAP1", "KIAA1549", "KIF5B", "KIT", "KLF4", "KLF6", "KLK2", "KMT2A", "KMT2C", "KMT2D", "KNL1", "KNSTRN", "KRAS", "KTN1",
                                                                            "LARP4B", "LASP1", "LATS1", "LATS2", "LCK", "LCP1", "LEF1", "LEPROTL1", "LHFPL6", "LIFR", "LMNA", "LMO1", "LMO2", "LPP", "LRIG3", "LRP1B", 
                                                                            "LSM14A", "LYL1", "LYN", "LZTR1", "MACC1", "MAF", "MAFB", "MALAT1", "MALT1", "MAML2", "MAP2K1", "MAP2K2", "MAP2K4", "MAP3K1", "MAP3K13", 
                                                                            "MAPK1", "MAX", "MB21D2", "MDM2", "MDM4", "MDS2", "MECOM", "MEN1", "MET", "MGMT", "MITF", "MLF1", "MLH1", "MLLT1", "MLLT10", "MLLT11",
                                                                            "MLLT3", "MLLT6", "MN1", "MNX1", "MPL", "MRTFA", "MSH2", "MSH6", "MSI2", "MTOR", "MUC1", "MUC16", "MUC4", "MUC6", "MUTYH", "MYB", "MYC",
                                                                            "MYCL", "MYCN", "MYD88", "MYH11", "MYH9", "MYO5A", "MYOD1", "N4BP2", "NAB2", "NACA", "NBEA", "NBN", "NCKIPSD", "NCOA1", "NCOA2", "NCOA4", 
                                                                            "NCOR1", "NCOR2", "NDRG1", "NF1", "NF2", "NFATC2", "NFE2L2", "NFIB", "NFKB2", "NFKBIE", "NIN", "NKX2-1", "NOTCH1", "NOTCH2", "NPM1", "NR2F2",
                                                                            "NR4A3", "NRAS", "NRG1", "NSD1", "NSD2", "NSD3", "NT5C2", "NTHL1", "NTRK1", "NTRK2", "NTRK3", "NUMA1", "NUP214", "NUP98", "NUTM1", "NUTM2D",
                                                                            "OLIG2", "OMD", "PABPC1", "PAFAH1B2", "PALB2", "PATZ1", "PAX3", "PAX5", "PAX7", "PAX8", "PBRM1", "PBX1", "PCBP1", "PCM1", "PDCD1LG2", "PDE4DIP",
                                                                            "PDGFB", "PDGFRA", "PDGFRB", "PER1", "PHOX2B", "PICALM", "PIK3CA", "PIK3CB", "PIK3R1", "PIK3R3", "PIM1", "PLAG1", "PLCG1", "PLEC", "PML", "PMS1",
                                                                            "PMS2", "POLD1", "POLE", "POLG", "POLQ", "POLR2A", "POT1", "POU2AF1", "POU5F1", "PPARG", "PPFIBP1", "PPM1D", "PPP2R1A", "PPP6C", "PRCC",
                                                                            "PRDM1", "PRDM16", "PRDM2", "PREX2", "PRF1", "PRKACA", "PRKAR1A", "PRKCB", "PRKD1", "PRPF40B", "PRRX1", "PSIP1", "PTCH1", "PTEN", "PTK6", 
                                                                            "PTPN11", "PTPN13", "PTPN6", "PTPRB", "PTPRC", "PTPRD", "PTPRK", "PTPRT", "PWWP2A", "QKI", "RABEP1", "RAC1", "RAD17", "RAD21", "RAD50",
                                                                            "RAD51B", "RAF1", "RALGDS", "RANBP2", "RAP1B", "RAP1GDS1", "RARA", "RB1", "RBM15", "RECQL4", "REL", "RET", "RFWD3", "RGPD3", "RGS7", "RHOA",
                                                                            "RHOH", "RMI2", "RNF213", "RNF43", "ROBO2", "ROS1", "RPL22", "RPL5", "RPN1", "RRAS2", "RSPO2", "RSPO3", "RUNX1", "RUNX1T1", "RXRA", "S100A7", 
                                                                            "SALL4", "SBDS", "SDC4", "SDHA", "SDHAF2", "SDHB", "SDHC", "SDHD", "SET", "SETBP1", "SETD1B", "SETD2", "SETDB1", "SF3B1", "SFPQ", "SFRP4",
                                                                            "SGK1", "SH2B3", "SH3GL1", "SHTN1", "SIRPA", "SIX1", "SIX2", "SKI", "SLC34A2", "SLC45A3", "SMAD2", "SMAD3", "SMAD4", "SMARCA4", "SMARCB1",
                                                                            "SMARCD1", "SMARCE1", "SMO", "SND1", "SNX29", "SOCS1", "SOX2", "SOX21", "SPECC1", "SPEN", "SPOP", "SRC", "SRGAP3", "SRSF2", "SRSF3", "SS18",
                                                                            "SS18L1", "STAG1", "STAT3", "STAT5B", "STAT6", "STIL", "STK11", "STRN", "SUB1", "SUFU", "SUZ12", "SYK", "TAF15", "TAL1", "TAL2", "TBL1XR1",
                                                                            "TBX3", "TCEA1", "TCF12", "TCF3", "TCF7L2", "TCL1A", "TEC", "TENT5C", "TERT", "TET1", "TET2", "TFEB", "TFG", "TFPT", "TFRC", "TGFBR2", "THRAP3", 
                                                                            "TLX1", "TLX3", "TMEM127", "TMPRSS2", "TNC", "TNFAIP3", "TNFRSF14", "TNFRSF17", "TNRC18", "TOP1", "TP53", "TP63", "TPM3", "TPM4", "TPR", 
                                                                            "TRAF7", "TRIM24", "TRIM27", "TRIM33", "TRIP11", "TRRAP", "TSC1", "TSC2", "TSHR", "U2AF1", "UBR5", "USP44", "USP6", "USP8", "VAV1", "VHL", 
                                                                            "VTI1A", "WDCP", "WIF1", "WNK2", "WRN", "WT1", "WWTR1", "XPA", "XPC", "XPO1", "YWHAE", "ZBTB16", "ZCCHC8", "ZEB1", "ZFHX3", "ZMYM2", "ZNF331",
                                                                            "ZNF384", "ZNF429", "ZNF479", "ZNF521", "ZNRF3"))

seqlevels(detail_regions_custom_2)
anno_dummy <- CNV.create_anno(array_type = c("450k", "EPICv2"), genome = "hg38")
seqlevels(anno_dummy@detail)

anno <- conumee2::CNV.create_anno(array_type = c("450k", "EPICv2"), genome= "hg38", bin_minprobes = 45,exclude_regions = exclude_regions, detail_regions = "detail_genes_4.bed")

#make useful names
targets <- read.metharray.sheet(idat_dir, pattern="targets_validation_with_predictions.csv")
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

x <- conumee2::CNV.focal(x,sig_cgenes = TRUE)


#####add custom genomeblot function with no detail text

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



names(x)
conumee2::CNV.genomeplot(x[39], sig_cgenes=FALSE, nsig_cgenes=3)

CNV.genomeplotcustom(x[25], sig_cgenes=FALSE, nsig_cgenes=3,bins_cex = 1.5)


CNV.plotly(x[9])


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

######################re-do for selected MNG genes
detail_regions_custom_4 = CNV.define_detail(array_type = "EPICv2", gene = c("CDKN2A", "CDKN2B", "NF2")) 

anno <- CNV.create_anno(array_type = c("450k", "EPICv2"), genome= "hg38", exclude_regions = exclude_regions, detail_regions = detail_regions_custom_4)

#perform CNV analysis

x <- CNV.fit(data.q, data.c, anno)
x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)

#detect focal CNVs
data(package = "conumee2")
data(cancer_genes_hg38)
x = CNV.focal(x, conf=0.99)




####make barplots of depleted/amplified cancer gene census genes
library(ggplot2)

#CancerDelGenes

Plot_CancerDelGenes = read.csv(file="CancerDelGenes_barplot.csv", header = T)
str(Plot_CancerDelGenes)
Plot_CancerDelGenes$condition = factor(Plot_CancerDelGenes$condition, 
                                       levels = c("shared_rec1_rec2","shared_pri_rec2","shared_pri_rec1","unique_rec2","unique_rec1","unique_pri","common"))

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

Plot_CancerAmpGenes = read.csv(file="CancerAmpGenes_barplot.csv", header = T)
str(Plot_CancerAmpGenes)
Plot_CancerAmpGenes$condition = factor(Plot_CancerAmpGenes$condition, 
                                       levels = c("shared_rec1_rec2","shared_pri_rec2","shared_pri_rec1","unique_rec2","unique_rec1","unique_pri","common"))

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

write.csv(res, file="Gene_CNV_enrichment_MNG_settings.csv")





#######make oncoprint for selected CNV genes
library(ComplexHeatmap)
#read in pri_no_rec, pri-rec, and rec individually and get optimized order
mat_pri_no_rec = read.table(file = "CNVs_validated_pri_no_rec.txt", header = TRUE, stringsAsFactors = FALSE, 
                            sep = "\t",na.strings = character(0),fill = TRUE,colClasses = "character")
rownames(mat_pri_no_rec) = mat_pri_no_rec[, 1]
mat_pri_no_rec = mat_pri_no_rec[, -1]
mat_pri_no_rec = as.matrix(mat_pri_no_rec)

mat_pri_rec = read.table(file = "CNVs_validated_pri_rec.txt", header = TRUE, stringsAsFactors = FALSE, 
                            sep = "\t",na.strings = character(0),fill = TRUE,colClasses = "character")
rownames(mat_pri_rec) = mat_pri_rec[, 1]
mat_pri_rec = mat_pri_rec[, -1]
mat_pri_rec = as.matrix(mat_pri_rec)

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
write.csv(sampleorder_pri_no_rec, file ="sampleorder_pri_no_rec.csv")

x=oncoPrint(mat_pri_rec,
            alter_fun = alter_fun_disc, col = col_disc, 
            heatmap_legend_param = heatmap_legend_param_disc, 
            alter_fun_is_vectorized = FALSE)

sampleorder_pri_rec = x@column_order
sampleorder_pri_rec
write.csv(sampleorder_pri_rec, file ="sampleorder_pri_rec.csv")

x=oncoPrint(mat_rec,
            alter_fun = alter_fun_disc, col = col_disc, 
            heatmap_legend_param = heatmap_legend_param_disc, 
            alter_fun_is_vectorized = FALSE)

sampleorder_rec = x@column_order
sampleorder_rec
write.csv(sampleorder_rec, file ="sampleorder_rec.csv")


#re-order full matrix
sampleorder_all = read.csv(file ="sampleorder_all.csv", header = T)
sampleorder_all = sampleorder_all$order


mat_all = read.table(file = "CNVs_validated_extended_all.txt", header = TRUE, stringsAsFactors = FALSE, 
                     sep = "\t",na.strings = character(0),fill = TRUE,colClasses = "character")
rownames(mat_all) = mat_all[, 1]
mat_all = mat_all[, -1]
mat_all = as.matrix(mat_all)
mat_all_order = mat_all[,sampleorder_all]


#get oncoprint

row_order = c("CCND1", "MDM2", "GATA3", "PIM1", "VHL", "CDKN2A/B", "NF2")

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



tiff("Oncoprint_full.tiff", res=300,width=5120,height=2048)
oncoPrint(mat_all_order,alter_fun = alter_fun_disc, col = col_disc,alter_fun_is_vectorized = FALSE,
          column_order = colnames(mat_all_order),row_order = row_order,
          heatmap_legend_param = heatmap_legend_param_disc,top_annotation = top_anno)
dev.off()


oncoPrint(mat_all_order,alter_fun = alter_fun_disc, col = col_disc,alter_fun_is_vectorized = FALSE,
          column_order = colnames(mat_disc_order),row_order = rownames(mat_disc_order),
          heatmap_legend_param = heatmap_legend_param_disc,top_annotation = top_anno)


###work on code for CNV plot for app


library(conumee2)
library(sesame)
library(sesameData)

anno = conumee2::CNV.create_anno(array_type = c("450k","EPICv2"), exclude_regions = exclude_regions, detail_regions = detail_regions)

#read in files
idat_dir = "/Users/lab/Desktop/Meningioma/data/idat_EPICv2"
sdfs.q <- openSesame(idat_dir, prep = "QCDPB", func = NULL)
intensity_q = totalIntensities(sdfs.q)
intensity_q <- data.frame(intensity_q)
colnames(intensity_q) <- "query"
intensity_q = subset(intensity_q, rownames(intensity_q) %in% anno@probes$IlmnID)
data.q <- CNV.load(intensity_q)

#read in control reference files
idat_dir_cont = "/Users/lab/Desktop/Meningioma/data/Capper_cont"
sdfs.c <- openSesame(idat_dir_cont, prep = "QCDPB", func = NULL)
intensity_c <- lapply(sdfs.c, totalIntensities)
intensity_c <- do.call(cbind, intensity_c)
intensity_c = subset(intensity_c, rownames(intensity_c) %in% names(anno@probes))
data.c <- CNV.load(intensity_c)

anno@probes = subset(anno@probes, names(anno@probes) %in% rownames())

#annotation with "overlap" to compare 450K and EPIC, define hg38 genome for EPICv2
data("exclude_regions")
data("detail_regions")

x=conumee2::CNV.fit(data.q, data.c, anno)


###new try
anno_450k = conumee2::CNV.create_anno(array_type = "450k")

# Automatically detect array type and prepare query for CNV.fit
library(dplyr)
library(S4Vectors)  # for names() on GRanges
library(sesame)     # openSesame / betas
library(conumee2)

prepare_query_for_cnv <- function(total_intensities) {
  # total_intensities: CNV.data@intensity or output from totalIntensities()
  # Returns a CNV.data object ready for CNV.fit with 450k controls
  
  # -------------------------
  # 1️⃣ Detect array type
  # -------------------------
  num_probes <- nrow(total_intensities)
  
  # Check suffix pattern if available
  probe_names <- rownames(total_intensities)
  if(any(grepl("_TC", probe_names))) {
    array_type <- "EPICv2"
  } else if(any(grepl("_BC", probe_names))) {
    array_type <- "EPIC"
  } else if(num_probes >= 937000) {
    array_type <- "EPICv2"
  } else if(num_probes >= 850000) {
    array_type <- "EPIC"
  } else {
    stop("Cannot determine array type from probe names or probe count.")
  }
  message("Detected array type: ", array_type)
  
  # -------------------------
  # 2️⃣ Load 450k-compatible annotation
  # -------------------------
  # Make sure anno object includes 'IlmID' column mapping EPIC/EPICv2 to 450k Probe_ID
  anno <- CNV.create_anno(array_type = c("450k", array_type))
  
  # -------------------------
  # 3️⃣ Map query probes to 450k-compatible IDs
  # -------------------------
  # Build a mapping vector: query probe -> 450k Probe_ID
  mapping <- setNames(anno@probes$Probe_ID, anno@probes$IlmID)
  
  common_probes <- intersect(rownames(total_intensities), names(mapping))
  if(length(common_probes) == 0) stop("No common probes found with annotation.")
  
  # Subset and rename rownames
  query_intensity <- total_intensities[common_probes, , drop = FALSE]
  rownames(query_intensity) <- mapping[common_probes]
  
  # -------------------------
  # 4️⃣ Replace CNV.data@intensity
  # -------------------------
  query_cnv <- CNV.data()
  query_cnv@intensity <- query_intensity
  
  return(query_cnv)
}

query_total_intensity = intensity_q

query_processed <- prepare_query_for_cnv(query_total_intensity, anno_450k)




