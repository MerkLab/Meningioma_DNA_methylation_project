library(stringr)
library(conumee2)


#read in files
idat_dir = "Directory_discovery"
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

#read in control reference files
idat_dir_cont = "Controls_Capper_paper"
cont = read.metharray.exp(base = idat_dir_cont)
cont
cont = preprocessNoob(cont)

#annotation with "overlap" to compare 450K and EPIC
data("exclude_regions")
anno = CNV.create_anno(array_type = "overlap", bin_minprobes = 50,
                       chrXY = FALSE,exclude_regions = exclude_regions, detail_regions = "Detail_NF2_CDKN2AB.bed")

#reduced MSet and cont to probes on anno
MSet = subset(MSet, rownames(MSet) %in% names(anno@probes))
anno@probes = subset(anno@probes, names(anno@probes) %in% rownames(MSet))
cont = subset(cont, rownames(cont) %in% names(anno@probes))

#make useful names
targets <- read.metharray.sheet(idat_dir, pattern="PC1_cluster_results.csv")
all(colnames(MSet) %in% targets$Basename)
all(colnames(MSet) == targets$Basename)
targets <- targets[match(colnames(MSet),targets$Basename),] 
all(colnames(MSet) == targets$Basename)
colnames(MSet) = targets$ID
colnames(MSet)

#bring into format
c = CNV.load(cont)
t = CNV.load(MSet)

#make CNV
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results/detail_2")

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
          detail.ratio, type = "p", pch = 16, col = "red")
    text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
         + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
         ifelse(detail.ratio.above, detail.ratio, NA), labels = paste("  ", 
                                                                      values(object@anno@detail)$name, sep = ""), adj = c(0, 
                                                                                                                          0.5), srt = 90, col = "red")
    text(start(object@anno@detail) + (end(object@anno@detail) - start(object@anno@detail)) /2
         + chr.cumsum0[as.vector(seqnames(object@anno@detail))], 
         ifelse(detail.ratio.above, NA, detail.ratio), labels = paste(values(object@anno@detail)$name, 
                                                                      "  ", sep = ""), adj = c(1, 0.5), srt = 90, col = "red")
  }
  
  if (set_par) 
    par(mfrow = mfrow_original, mar = mar_original, oma = oma_original)
})



#make the plots and get text files for CNVs
for(x in names(t)){
  tmp <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(t[x], c,anno))))
  tiff(paste0(x,".tiff"),compression="lzw",res=300,width=5120,height=2048);
  CNV.genomeplot(tmp, ylim = c(-1.1, 1.1), cols = c("blue4", "blue2", "lightgrey", "darkorange2", "darkorange4"));dev.off()
  write.table(CNV.write(tmp,what="segments"),sep="\t",quote=F,row.names=F,file=paste0(x,".tsv"))
}

#make CNV plots for CDKN2A/B und NF2
for(x in names(t)){
  tmp <- CNV.segment(CNV.detail(CNV.bin(CNV.fit(t[x], c,anno))))
  tiff(paste0(x,".tiff"),compression="lzw",res=300,width=5120,height=2048);
  CNV.genomeplotcustom(tmp, ylim = c(-1.1, 1.1), cols = c("blue4", "blue2", "lightgrey", "darkorange2", "darkorange4"));dev.off()
}

save.image()



#coerce to one dataframe to work with
library("stringr")

# get locations of all cn files
files <- Sys.glob("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results/*.tsv")

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


##calculate genome instability for each sample as percentage of genome affected by CNVs, cutoff 0.3!
getwd()
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
CNV_sizes = read.csv(file = "All_samples_CNV_sizes.csv", header = T)
str(CNV_sizes)
CNV_sizes = aggregate(CNV_sizes[-1], by = list(CNV_sizes$ID), FUN=sum)
write.csv(CNV_sizes, file = "CNV_sizes_aggregate.csv")

genome_instability = read.csv(file="Percentage_genome_disrupted.csv", header = T)
all(genome_instability$ID %in% sampleorder)
all(genome_instability$ID == sampleorder_T26) #order fits to oncoprint

##########work on generating an oncoprint for most recurrent (>10% of sampels) CNVs

#T26 discovery all
library(ComplexHeatmap)
getwd()
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")

mat_disc = read.table(file = "Broad_CNVs_recurrent_discovery.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
rownames(mat_disc) = mat_disc[, 1]
mat_disc = mat_disc[, -1]
mat_disc = t(as.matrix(mat_disc))


mat_c1 = read.table(file = "Broad_CNVs_cluster1.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
rownames(mat_c1) = mat_c1[, 1]
mat_c1 = mat_c1[, -1]
mat_c1 = t(as.matrix(mat_c1))

mat_c2 = read.table(file = "Broad_CNVs_cluster2.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
rownames(mat_c2) = mat_c2[, 1]
mat_c2 = mat_c2[, -1]
mat_c2 = t(as.matrix(mat_c2))

col_disc = c("DEL" = "blue4", "AMP" = "darkorange2")

alter_fun_disc = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  DEL = alter_graphic("rect", fill = col_disc["DEL"]),
  AMP = alter_graphic("rect", fill = col_disc["AMP"])
)

heatmap_legend_param_disc = list(title = "CNVs", at = c("DEL", "AMP"), 
                                 labels = c("Deletion", "Amplification"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/oncoprint_CNVs")

tiff("Oncoprint_c1.tiff", res=300,width=5120,height=2048)
oncoPrint(mat_c1,
          alter_fun = alter_fun_disc, col = col_disc, 
          column_title = column_title_disc, heatmap_legend_param = heatmap_legend_param_disc, 
          alter_fun_is_vectorized = FALSE)
dev.off()

tiff("Oncoprint_c2.tiff", res=300,width=1837,height=2048)
oncoPrint(mat_c2,
          alter_fun = alter_fun_disc, col = col_disc, 
          column_title = column_title_disc, heatmap_legend_param = heatmap_legend_param_disc, 
          alter_fun_is_vectorized = FALSE)
dev.off()

#get optimzed sample order for c1 and c2

x = oncoPrint(mat_c1,
              alter_fun = alter_fun_disc, col = col_disc, 
              column_title = column_title_disc, heatmap_legend_param = heatmap_legend_param_disc, 
              alter_fun_is_vectorized = FALSE)
y = oncoPrint(mat_c2,
              alter_fun = alter_fun_disc, col = col_disc, 
              column_title = column_title_disc, heatmap_legend_param = heatmap_legend_param_disc, 
              alter_fun_is_vectorized = FALSE)

sampleorder_c1 = x@column_order
sampleorder_c1
write.csv(sampleorder_c1, file ="sampleorder_c1.csv")

sampleorder_c2 = y@column_order
sampleorder_c2
write.csv(sampleorder_c2, file ="sampleorder_c2.csv")

#re-order full matrix
sampleorder_c1 = read.csv(file ="sampleorder_c1.csv", header = F)
sampleorder_c1 = sampleorder_c1$V1
sampleorder_c2 = read.csv(file ="sampleorder_c2.csv", header = F)
sampleorder_c2 = sampleorder_c2$V1
sampleorder_T26 = read.csv(file ="sampleorder_T26.csv", header = F)
sampleorder_T26 = sampleorder_T26$V1

mat_disc_order = mat_disc[,sampleorder_T26]


#make oncoprint, also add several annotations (genome instabilty, cluster assignment, grade, mol group)
#get annotation data
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/oncoprint_CNVs")
anno_genome_instability = genome_instability$Percentage
oncoprint_anno = read.csv(file ="oncoprint_anno.csv", header = T)

top_anno =HeatmapAnnotation(df = oncoprint_anno,
                            points = anno_barplot(anno_genome_instability,gp=gpar(border =NA,fill="deeppink3",
                                                                                  lty="blank"),border = T,height = unit(2, "cm")),
                            col = list(cluster = c("1" = "cyan4", "2" = "violetred4"),
                                       grade = c("1" = "#9EBCDA", "2" = "#8C6BB1", "3" = "#810F7C"),
                                       group = c("Merlinintact" = "royalblue2","Immuneenriched" = "red3",
                                                 "hypermetabolic" = "forestgreen","proliferative" = "darkorange2"),
                                       risk = c("low"="dodgerblue2", "intermediate" ="purple1","high" = "firebrick1")))



increase_incidence = read.csv(file="Incidence_increase_cluster2.csv", header = T)
increase_incidence = increase_incidence$increase
row_anno = rowAnnotation(incidence = anno_barplot(increase_incidence,gp=gpar(border =NA,fill="deepskyblue4",
                                                                             lty="blank"),border = T,height = unit(2, "cm")))


tiff("Oncoprint_full.tiff", res=300,width=5120,height=2048)
oncoPrint(mat_disc_order,alter_fun = alter_fun_disc, col = col_disc,alter_fun_is_vectorized = FALSE,
          column_order = colnames(mat_disc_order),row_order = rownames(mat_disc_order),
          heatmap_legend_param = heatmap_legend_param_disc,top_annotation = top_anno,right_annotation = row_anno)
dev.off()




######check distribution of genome instability across grading, groups, and cluster
#grade
getwd()
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
GeIn_grade_T26 = read.csv("instability_grade.csv", header = T)
str(GeIn_grade_T26)
GeIn_grade_T26$condition = factor(GeIn_grade_T26$condition, levels = c("1", "2", "3"))

library(ggplot2)
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "Instability_T26_grade.pdf", width = 4, height = 9)
ggplot(GeIn_grade_T26, 
       aes(condition, instability, fill=condition))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("#9EBCDA","#8C6BB1","#810F7C")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#groups
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
GeIn_groups_T26 = read.csv("instability_group.csv", header = T)
str(GeIn_groups_T26)
GeIn_groups_T26$condition = factor(GeIn_groups_T26$condition, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))
str(GeIn_groups_T26)

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "Instability_T26_groups.pdf", width = 5, height = 9)
ggplot(GeIn_groups_T26, 
       aes(condition, instability, fill=condition))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("royalblue2","red3","forestgreen","darkorange2")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#cluster
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
GeIn_cluster_T26 = read.csv("instability_cluster.csv", header = T)
str(GeIn_cluster_T26)
GeIn_cluster_T26$condition = factor(GeIn_cluster_T26$condition, levels = c("1", "2"))
str(GeIn_cluster_T26)

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "Instability_T26_cluster.pdf", width = 3, height = 9)
ggplot(GeIn_cluster_T26, 
       aes(condition, instability, fill=condition))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("cyan4","violetred4")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#risk
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
GeIn_risk_T26 = read.csv("instability_risk.csv", header = T)
str(GeIn_risk_T26)
GeIn_risk_T26$condition = factor(GeIn_risk_T26$condition, levels = c("low", "intermediate", "high"))
str(GeIn_risk_T26)

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "Instability_T26_risk.pdf", width = 4, height = 9)
ggplot(GeIn_risk_T26, 
       aes(condition, instability, fill=condition))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("dodgerblue2","purple1","firebrick1")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


####check genome instability in intermediate MNG (grade 2, hypermetabolic, intermediate risk) stratified by cluster assignment
#difference in genome instability in grade 2 MNG stratified by cluster
getwd()
GeIn_grade2_cluster = read.csv("stability_grade2_cluster.csv", header = T)
GeIn_grade2_cluster$condition = factor(GeIn_grade2_cluster$condition)
GeIn_grade2_cluster$cluster = factor(GeIn_grade2_cluster$cluster, levels = c("1", "2"))
str(GeIn_grade2_cluster)

cairo_pdf(filename = "Instability_grade2_cluster.pdf", width = 3, height = 9)
ggplot(GeIn_grade2_cluster, 
       aes(condition, instability, fill=interaction(condition,cluster), dodge=cluster))+
  geom_point(position = position_jitterdodge(0.4), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("#8C6BB166","#8C6BB1")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#difference in genome instability in hypermetabolic MNG stratified by cluster
GeIn_group_cluster = read.csv("Instability_group_cluster.csv", header = T)
GeIn_group_cluster$condition = factor(GeIn_group_cluster$condition)
GeIn_group_cluster$cluster = factor(GeIn_group_cluster$cluster, levels = c("1", "2"))
str(GeIn_group_cluster)

cairo_pdf(filename = "Instability_group_cluster.pdf", width = 3, height = 9)
ggplot(GeIn_group_cluster, 
       aes(condition, instability, fill=interaction(condition,cluster), dodge=cluster))+
  geom_point(position = position_jitterdodge(0.4), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("#228B2266","#228B22")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#difference in genome instability in intermediate risk MNG stratified by cluster
GeIn_risk_cluster = read.csv("Instability_risk_cluster.csv", header = T)
GeIn_risk_cluster$condition = factor(GeIn_risk_cluster$condition)
GeIn_risk_cluster$cluster = factor(GeIn_risk_cluster$cluster, levels = c("1", "2"))
str(GeIn_risk_cluster)

cairo_pdf(filename = "Instability_risk_cluster.pdf", width = 3, height = 9)
ggplot(GeIn_risk_cluster, 
       aes(condition, instability, fill=interaction(condition,cluster), dodge=cluster))+
  geom_point(position = position_jitterdodge(0.4), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("plum3","purple1")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

######check number of unfavorable CNVs across grading, groups, and cluster
#grade
library(ggplot2)

getwd()
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
CNVunfavor_grade_T26 = read.csv("CNV_unfavor_grade.csv", header = T)
str(CNVunfavor_grade_T26)
CNVunfavor_grade_T26$condition = factor(CNVunfavor_grade_T26$condition, levels = c("1", "2", "3"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "CNVs_unfavorable_T26_grade.pdf", width = 4, height = 9)
ggplot(CNVunfavor_grade_T26, 
       aes(condition, CNV_unfavor, fill=condition))+
  geom_jitter(height = 0, width = 0.1, color ="grey50")+
  geom_boxplot(outlier.shape = NA,alpha=0.8)+
  scale_fill_manual(values=c("#9EBCDA","#8C6BB1","#810F7C")) +
  scale_y_continuous(breaks = seq(0,7,1), limits=c(0,7))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()



#group
getwd()
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
CNVunfavor_group_T26 = read.csv("CNV_unfavor_group.csv", header = T)
str(CNVunfavor_group_T26)
CNVunfavor_group_T26$condition = factor(CNVunfavor_group_T26$condition, levels = c("Merlinintact", "Immuneenriched", "hypermetabolic", "proliferative"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "CNVs_unfavorable_T26_group.pdf", width = 5, height = 9)
ggplot(CNVunfavor_group_T26, 
       aes(condition, CNV_unfavor, fill=condition))+
  geom_jitter(height = 0, width = 0.1, color ="grey50")+
  geom_boxplot(outlier.shape = NA,alpha=0.8)+
  scale_fill_manual(values=c("royalblue2","red3","forestgreen","darkorange2")) +
  scale_y_continuous(breaks = seq(0,7,1), limits=c(0,7))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#risk
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
CNVunfavor_risk_T26 = read.csv("CNV_unfavor_risk.csv", header = T)
str(CNVunfavor_risk_T26)
CNVunfavor_risk_T26$condition = factor(CNVunfavor_risk_T26$condition, levels = c("low", "intermediate", "high"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "CNVs_unfavorable_T26_risk.pdf", width = 4, height = 9)
ggplot(CNVunfavor_risk_T26, 
       aes(condition, CNV_unfavor, fill=condition))+
  geom_jitter(height = 0, width = 0.1, color ="grey50")+
  geom_boxplot(outlier.shape = NA,alpha=0.8)+
  scale_fill_manual(values=c("dodgerblue2","purple1","firebrick1")) +
  scale_y_continuous(breaks = seq(0,7,1), limits=c(0,7))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#cluster
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
CNVunfavor_cluster_T26 = read.csv("CNV_unfavor_cluster.csv", header = T)
str(CNVunfavor_cluster_T26)
CNVunfavor_cluster_T26$condition = factor(CNVunfavor_cluster_T26$condition, levels = c("1", "2"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "CNVs_unfavorable_T26_cluster.pdf", width = 3, height = 9)
ggplot(CNVunfavor_cluster_T26, 
       aes(condition, CNV_unfavor, fill=condition))+
  geom_jitter(height = 0, width = 0.1, color ="grey50")+
  geom_boxplot(outlier.shape = NA,alpha=0.8)+
  scale_fill_manual(values=c("cyan4","violetred4")) +
  scale_y_continuous(breaks = seq(0,7,1), limits=c(0,7))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()



###check influence of 1p and 22q deletions on genome stability
##in all T26
#regarding instability
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
x = read.csv("Instability_T26_1p22q.csv", header = T)
str(x)
x$condition = factor(x$condition, levels = c("other", "1p", "22q", "codeletion"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "Instability_T26_1p-22q.pdf", width = 5, height = 9)
ggplot(x, 
       aes(condition, instability, fill=condition))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("gray30","deepskyblue","deepskyblue3","dodgerblue3")) +
  scale_y_continuous(breaks = seq(0,0.5,0.25), limits=c(0,0.5))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#regarding unfavorable CNVs
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
x = read.csv("CNV_unfavor_1p22q.csv", header = T)
str(x)
x$condition = factor(x$condition, levels = c("other", "1p", "22q", "codeletion"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "CNV_unfavorable_T26_1p-22q.pdf", width = 5, height = 9)
ggplot(x, 
       aes(condition, CNV_favor, fill=condition))+
  geom_jitter(height = 0, width = 0.1, color ="grey50")+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("gray30","deepskyblue","deepskyblue3","dodgerblue3")) +
  scale_y_continuous(breaks = seq(0,7,1), limits=c(0,7))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


##in METH low
#regarding instability
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
x = read.csv("Instability_c1_1p22q.csv", header = T)
str(x)
x$condition = factor(x$condition, levels = c("other", "1p", "22q", "codeletion"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "Instability_c1_1p-22q.pdf", width = 5, height = 9)
ggplot(x, 
       aes(condition, instability, fill=condition))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50", size=2)+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("cyan4","cyan4","cyan4","cyan4")) +
  scale_y_continuous(breaks = seq(0,0.4,0.1), limits=c(0,0.43))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#regarding unfavorable CNVs
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
x = read.csv("CNV_unfavorable_c1_1p22q.csv", header = T)
str(x)
x$condition = factor(x$condition, levels = c("other", "1p", "22q", "codeletion"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "CNV_unfavorable_c1_1p-22q.pdf", width = 5, height = 9)
ggplot(x, 
       aes(condition, CNV_favor, fill=condition))+
  geom_jitter(height = 0, width = 0.1, color ="grey50", size=2)+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("cyan4","cyan4","cyan4","cyan4")) +
  scale_y_continuous(breaks = seq(0,7,1), limits=c(0,7))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()



##in METH high
#regarding instability
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
x = read.csv("Instability_c2_1p22q.csv", header = T)
str(x)
x$condition = factor(x$condition, levels = c("other", "1p", "22q", "codeletion"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "Instability_c2_1p-22q.pdf", width = 4, height = 9)
ggplot(x, 
       aes(condition, instability, fill=condition))+
  geom_point(position = position_jitter(0.2), alpha=0.95, color="grey50", size=2)+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("violetred4","violetred4","violetred4","violetred4")) +
  scale_y_continuous(breaks = seq(0,0.4,0.1), limits=c(0,0.43))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

#regarding unfavorable CNVs
setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/data")
x = read.csv("CNV_unfavorable_c2_1p22q.csv", header = T)
str(x)
x$condition = factor(x$condition, levels = c("other", "1p", "22q", "codeletion"))

setwd("/Users/lab/Desktop/Meningioma/T26_CNV_samples/results")
cairo_pdf(filename = "CNV_unfavorable_c2_1p-22q.pdf", width = 4, height = 9)
ggplot(x, 
       aes(condition, CNV_favor, fill=condition))+
  geom_jitter(height = 0, width = 0.1, color ="grey50", size=2)+
  geom_boxplot(outlier.shape = NA, alpha=0.8)+
  scale_fill_manual(values=c("violetred4","violetred4","violetred4","violetred4")) +
  scale_y_continuous(breaks = seq(0,7,1), limits=c(0,7))+
  theme_classic()+
  theme(legend.position = "none")
dev.off()
