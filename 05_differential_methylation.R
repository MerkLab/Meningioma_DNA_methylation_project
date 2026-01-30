library(sesame)
library(ggplot2)
library(SummarizedExperiment)
library(dplyr)
library(plyr)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(circlize)
library(rockchalk)

#get discovery data, preprocessed and batch corrected
beta.disc = readRDS(file = "Betas_discovery_preprocessed_corrected.rds")
dim(beta.disc)
meta.disc = read.csv(file="meta_discovery_current.csv")
rownames(meta.disc) = meta.disc$ID

# differential methylation in discovery cohort
#first, make summarizedExperiment
se.discovery <- SummarizedExperiment(beta.disc, colData = meta.disc)
str(se.discovery)
dim(se.discovery)
colData(se.discovery)
assay(se.discovery)

####determine differentially methylated probes for several predictors of interest (cluster, grading (grade 2 and 3 vs 1), and hypermetabolic and proliferative vs NF2-wildtype/immuneenriched)
#####prepare predictor variables
#turn discrete contrast variables to factor and relevel
colData(se.discovery)$cluster_1000_new <- relevel(factor(colData(se.discovery)$cluster_1000_new), "METHlow")
str(se.discovery$cluster_1000_new)


###model methylation differences
####for methylation cluster
smry.discovery.cluster = DML(se.discovery, ~cluster_1000_new)
lm_results_discovery_cluster = summaryExtractTest(smry.discovery.cluster)
str(lm_results_discovery_cluster)

# adjust p values for status
padj.c2vsc1 = p.adjust(lm_results_discovery_cluster$Pval_cluster_1000_newMETHhigh, method = "BH", n = length(lm_results_discovery_cluster$Pval_cluster_1000_newMETHhigh))
lm_results_discovery_cluster = cbind(lm_results_discovery_cluster, padj.c2vsc1)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in METHhigh vs METHlow
#442 hypomethylated, 6,289 hypermethylated in c2

hypo_c2_vs_c1 = lm_results_discovery_cluster %>% dplyr::filter(Est_cluster_1000_newMETHhigh < -0.2 & padj.c2vsc1 < 0.05)

hyper_c2_vs_c1 = lm_results_discovery_cluster %>% dplyr::filter(Est_cluster_1000_newMETHhigh > 0.2 & padj.c2vsc1 < 0.05)
write.csv(hyper_c2_vs_c1, file="probes_hyper_cluster_disc.csv")


#visualize differences in methylation in volcano for cluster

cairo_pdf(filename = "Volcano_METHhigh_vs_METHlow_probe_level_clean.pdf", width = 8, height = 6)
ggplot(lm_results_discovery_cluster) + geom_point(aes(x=Est_cluster_1000_newMETHhigh, y=-log10(padj.c2vsc1), color=-log10(padj.c2vsc1), size=abs(Est_cluster_1000_newMETHhigh)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.05, 0)) + 
  scale_radius(range = c(0.05,5))+
  scale_y_continuous(breaks = seq(0,40,10), limits=c(0,40))+
  scale_x_continuous(breaks = seq(-0.5,0.5,0.25), limits=c(-0.5,0.5))+
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = -0.2 , y = 1.3, xend = -0.2, yend = Inf), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = Inf, yend = 1.3), colour = "grey35", linetype = 2) +
  geom_segment(aes(x = 0.2 , y = 1.3, xend = 0.2, yend = Inf), colour = "grey35", linetype = 2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position="none")
dev.off()



#subset summarizedexperiment for individual contrasts not using the entire discovery cohort
se.discovery_for_grade_2 = se.discovery[,!se.discovery$grading2021_new == "3"]
dim(se.discovery_for_grade_2)
colData(se.discovery_for_grade_2)$grading2021_new <- relevel(factor(colData(se.discovery_for_grade_2)$grading2021_new), "1")
str(se.discovery_for_grade_2)

se.discovery_for_grade_3 = se.discovery[,!se.discovery$grading2021_new == "2"]
dim(se.discovery_for_grade_3)
colData(se.discovery_for_grade_3)$grading2021_new <- relevel(factor(colData(se.discovery_for_grade_3)$grading2021_new), "1")
str(se.discovery_for_grade_3)


se.discovery_for_hypermetabolic = se.discovery[, colData(se.discovery)$mol_group != "proliferative"]
colData(se.discovery_for_hypermetabolic)$mol_group <-factor(colData(se.discovery_for_hypermetabolic)$mol_group)
colData(se.discovery_for_hypermetabolic)$mol_group <-
  combineLevels(
    colData(se.discovery_for_hypermetabolic)$mol_group,
    levs = c("Merlinintact", "Immuneenriched"),
    newLabel = "other"
  )
colData(se.discovery_for_hypermetabolic)$mol_group <- relevel(factor(colData(se.discovery_for_hypermetabolic)$mol_group), "other")


se.discovery_for_proliferative = se.discovery[, colData(se.discovery)$mol_group != "hypermetabolic"]
colData(se.discovery_for_proliferative)$mol_group <-factor(colData(se.discovery_for_proliferative)$mol_group)
colData(se.discovery_for_proliferative)$mol_group <-
  combineLevels(
    colData(se.discovery_for_proliferative)$mol_group,
    levs = c("Merlinintact", "Immuneenriched"),
    newLabel = "other"
  )
colData(se.discovery_for_proliferative)$mol_group <- relevel(factor(colData(se.discovery_for_proliferative)$mol_group), "other")




#model differences grade2 vs grade1
smry.discovery.grade2vs1 = DML(se.discovery_for_grade_2, ~grading2021_new)
lm_results_discovery_grade2vs1 = summaryExtractTest(smry.discovery.grade2vs1)
str(lm_results_discovery_grade2vs1)

# adjust p values for grade2 vs grade1
padj.grade2 = p.adjust(lm_results_discovery_grade2vs1$Pval_grading2021_new2, method = "BH", n = length(lm_results_discovery_grade2vs1$Pval_grading2021_new2))
lm_results_discovery_grade2vs1 = cbind(lm_results_discovery_grade2vs1, padj.grade2)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade2 over grade1
# 1 hypomethylated,  45 hypermethylated probes in grade2 vs grade1
hypo_grade2 = lm_results_discovery_grade2vs1 %>% dplyr::filter(Est_grading2021_new2 < -0.2 & padj.grade2 < 0.05)
hyper_grade2 = lm_results_discovery_grade2vs1 %>% dplyr::filter(Est_grading2021_new2 > 0.2 & padj.grade2 < 0.05)
write.csv(hyper_grade2, file="probes_hyper_grade2_disc.csv")




#model differences grade3 vs grade1
smry.discovery.grade3vs1 = DML(se.discovery_for_grade_3, ~grading2021_new)
lm_results_discovery_grade3vs1 = summaryExtractTest(smry.discovery.grade3vs1)
str(lm_results_discovery_grade3vs1)

# adjust p values for grade2 vs grade1
padj.grade3 = p.adjust(lm_results_discovery_grade3vs1$Pval_grading2021_new3, method = "BH", n = length(lm_results_discovery_grade3vs1$Pval_grading2021_new3))
lm_results_discovery_grade3vs1 = cbind(lm_results_discovery_grade3vs1, padj.grade3)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade2 over grade1
# 1 hypomethylated,  45 hypermethylated probes in grade2 vs grade1
hypo_grade3 = lm_results_discovery_grade3vs1 %>% dplyr::filter(Est_grading2021_new3 < -0.2 & padj.grade3 < 0.05)
hyper_grade3 = lm_results_discovery_grade3vs1 %>% dplyr::filter(Est_grading2021_new3 > 0.2 & padj.grade3 < 0.05)
write.csv(hyper_grade3, file="probes_hyper_grade3_disc.csv")




#model differences hypermetabolic vs NF2-WT/Immunogenic
smry.discovery.hypermetabolic = DML(se.discovery_for_hypermetabolic, ~mol_group)
lm_results_discovery_hypermetabolic = summaryExtractTest(smry.discovery.hypermetabolic)
str(lm_results_discovery_hypermetabolic)

# adjust p values for hypermetabolic vs reference
padj.hypermetabolic = p.adjust(lm_results_discovery_hypermetabolic$Pval_mol_grouphypermetabolic, method = "BH", n = length(lm_results_discovery_hypermetabolic$Pval_mol_grouphypermetabolic))
lm_results_discovery_hypermetabolic = cbind(lm_results_discovery_hypermetabolic, padj.hypermetabolic)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in hypermetabolic vs reference
#  hypomethylated,   hypermethylated probes in grade3 vs reference
hypo_hypermetabolic = lm_results_discovery_hypermetabolic %>% dplyr::filter(Est_mol_grouphypermetabolic < -0.2 & padj.hypermetabolic < 0.05)
hyper_hypermetabolic = lm_results_discovery_hypermetabolic %>% dplyr::filter(Est_mol_grouphypermetabolic > 0.2 & padj.hypermetabolic < 0.05)
write.csv(hyper_hypermetabolic, file="probes_hyper_hypermetabolic_disc.csv")



#model differences proliferative vs NF2-WT/Immunogenic
smry.discovery.proliferative = DML(se.discovery_for_proliferative, ~mol_group)
lm_results_discovery_proliferative = summaryExtractTest(smry.discovery.proliferative)
str(lm_results_discovery_proliferative)

# adjust p values for hypermetabolic vs reference
padj.proliferative = p.adjust(lm_results_discovery_proliferative$Pval_mol_groupproliferative, method = "BH", n = length(lm_results_discovery_proliferative$Pval_mol_groupproliferative))
lm_results_discovery_proliferative = cbind(lm_results_discovery_proliferative, padj.proliferative)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in hypermetabolic vs reference
#  hypomethylated,   hypermethylated probes in grade3 vs reference
hypo_proliferative = lm_results_discovery_proliferative %>% dplyr::filter(Est_mol_groupproliferative < -0.2 & padj.proliferative < 0.05)
hyper_proliferative = lm_results_discovery_proliferative %>% dplyr::filter(Est_mol_groupproliferative > 0.2 & padj.proliferative < 0.05)
write.csv(hyper_proliferative, file="probes_hyper_proliferative_disc.csv")






#calculate DMRs
dmContrasts(smry.discovery.cluster)
DMR_cluster = DMR(se.discovery, smry.discovery.cluster, "cluster_1000_newMETHhigh",seg.per.locus = 0.1)
DMR_cluster_hyper = DMR_cluster %>% dplyr::filter(Seg_Est > 0.2 & Seg_Pval_adj < 0.05)
DMR_cluster_hypo = DMR_cluster %>% dplyr::filter(Seg_Est < -0.2 & Seg_Pval_adj < 0.05)

dmr_avg_hyper <- aggregate(
  Seg_Est ~ Seg_Chrm + Seg_Start + Seg_End,
  data = DMR_cluster_hyper,
  FUN = mean
)

dmr_avg_hypo <- aggregate(
  Seg_Est ~ Seg_Chrm + Seg_Start + Seg_End,
  data = DMR_cluster_hypo,
  FUN = mean
)

write.csv(dmr_avg_hyper, file="DMR_hyper.csv")
write.csv(dmr_avg_hypo, file="DMR_hypo.csv")








#####gene enrichment in differentially methylated probes
#test hypo probes in METHhigh
genes.hypo.c2 = testEnrichment(hypo_c2_vs_c1$Probe_ID, KYCG_buildGeneDBs(rownames(beta.disc)))
genes.hypo.c2.sig = genes.hypo.c2 %>% dplyr::filter(FDR < 0.01)
genes.hypo.c2.sig = genes.hypo.c2.sig$gene_name
genes.hypo.c2.sig.df = as.data.frame(genes.hypo.c2.sig)
write.csv(genes.hypo.c2, file="Cluster2_hypomethlyted_genes.csv")
write.csv(genes.hypo.c2.sig, file = "Cluster2_hypomethylated_genes_sig.csv")


#test hyper probes in METHhigh
genes.hyper.c2 = testEnrichment(hyper_c2_vs_c1$Probe_ID, KYCG_buildGeneDBs(rownames(beta.disc)))
genes.hyper.c2.sig = genes.hyper.c2 %>% dplyr::filter(FDR < 0.01)
genes.hyper.c2.sig = genes.hyper.c2.sig$gene_name
genes.hyper.c2.sig.df = as.data.frame(genes.hyper.c2.sig)
write.csv(genes.hyper.c2, file = "Cluster2_hypermethylated_genes.csv")
write.csv(genes.hyper.c2.sig, file = "Cluster2_hypermethylated_genes_sig.csv")


#make Manhatten plot for hypermethylated genes in cluster c2
#get p values ordered by p value
c2_hyper_Man_data = read.csv(file="c2_hypermethylated_Manhatten_mod.csv", header = T)
axisdf = read.csv(file="axisdf.csv")

cairo_pdf(filename = "Manhatten_genes_enriched_hypermethylated_c2.pdf", width = 9, height = 5)
ggplot(c2_hyper_Man_data, aes(x=BPcum, y=P)) +
  
  # Show all points
  geom_point(aes(color=as.factor(CHR)), alpha=0.3, size=3) +
  scale_color_manual(values = rep(c("grey23", "grey60"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 1))+
  ylim(0,69)+
  
  # Add highlighted points
  geom_point(data=subset(c2_hyper_Man_data, PCDHG=="yes"), color="red", size=6,alpha=1) +
  geom_point(data=subset(c2_hyper_Man_data, PCDHB=="yes"), color="darkorange", size=6,alpha=1) +
  geom_point(data=subset(c2_hyper_Man_data, PCDHA=="yes"), color="purple3", size=6,alpha=1) +
  geom_point(data=subset(c2_hyper_Man_data, HOXA=="yes"), color="deepskyblue1", size=6) +
  geom_point(data=subset(c2_hyper_Man_data, HOXB=="yes"), color="deepskyblue2", size=6) +
  geom_point(data=subset(c2_hyper_Man_data, HOXC=="yes"), color="deepskyblue3", size=6) +
  geom_point(data=subset(c2_hyper_Man_data, HOXD=="yes"), color="deepskyblue4", size=6) +
  
  # Custom the theme:
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.position = "none")
dev.off()







#make plot for 23 hypomethylated genes
hypo_genes = read.csv(file = "cluster2_hypo_genes_plot.csv")

cairo_pdf(filename = "Hypo_genes_barplot.pdf", width = 15, height = 4)
ggplot(hypo_genes, aes(x = gene, y = padj)) +
  geom_col(fill = "steelblue", width = 0.7) +
  theme_minimal(base_size = 13) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )
dev.off()









#####perform DAVID functional annotation clustering
###make a dotplot, with coloring showing enrichment and size showing number of genes

hyper_genes_groups = read.csv(file="Hyper_genes_functional_groups.csv", header = T)
str(hyper_genes_groups)
hyper_genes_groups$group = factor(hyper_genes_groups$group, levels = c("Homeobox", "Protocadherin", "Thyroid gland development", "bHLH",
                                                                       "T cell receptor complex", "Sox family"))

cairo_pdf(filename = "Hyper_genes_enrichment_groups_2.pdf", width = 13, height = 4)
ggplot(hyper_genes_groups, aes(x=group, y=enrichment))+
  geom_point(aes(colour = enrichment, size = number_genes))+
  scale_colour_gradient(                          
    low  = "#4575b4",
    high = "#d73027"
  )+
  coord_flip()+
  scale_size(range = c(5, 15),breaks = c(1, 5, 10, 20, 40))+
  theme_classic()
dev.off()








####generate circos plot to visualize differentially methylated CpGs and DMRs
annoEPICv2 = getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
anno = as.data.frame(annoEPICv2$Name)
anno$chr = annoEPICv2$chr
anno$start = annoEPICv2$pos
write.csv(anno, file="EPICv2_anno.csv")

cpgs = read.csv(file="CpGs_for_circos.csv")
dmrs = read.csv(file="DMRs_for_circos.csv")

cpgs <- cpgs %>%
  mutate(
    start = as.numeric(as.character(start)),
    end   = as.numeric(as.character(end))
  )

dmrs <- dmrs %>%
  mutate(
    start = as.numeric(as.character(start)),
    end   = as.numeric(as.character(end))
  )

chr_use <- paste0("chr", 1:22)
cpgs <- cpgs %>%
  filter(chr %in% chr_use)
dmrs <- dmrs %>%
  filter(chr %in% chr_use)

#make plot
cairo_pdf(filename = "Ciros_diff_methylation.pdf", width = 10, height = 10)
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 1,
  track.margin = c(0.005, 0.005),
  cell.padding = c(0, 0, 0, 0),
  gap.after = c(rep(1, 21), 10) 
)

circos.initializeWithIdeogram(
  species = "hg38",
  chromosome.index = chr_use
)

## Shared y-scale for methylation
ylim_meth <- c(-0.5, 0.5)
scale_ticks <- c(-0.5, -0.25, 0, 0.25, 0.5)

#outer track
circos.genomicTrack(
  cpgs,
  ylim = c(-0.5, 0.7),
  track.height = 0.2,
  track.margin = c(0.0, 0.02),
  bg.border = "black",
  panel.fun = function(region, value, ...) {
    
    db <- value[[1]]
    
    col <- ifelse(
      db > 0,
      "red",   # hyper
      "navy"    # hypo
    )
    
    # Draw points
    circos.genomicPoints(
      region,
      db,
      pch = 16,
      cex = 0.5,
      col = col
    )
    
    # Draw zero baseline for all chromosomes
    circos.lines(
      x = CELL_META$xlim, 
      y = c(0, 0),        
      col = "black",
      lwd = 0.5
    )
    
    # Draw y-axis only for chr1
    current_chr <- CELL_META$sector.index
    if (current_chr %in% c("1", "chr1")) {
      circos.yaxis(
        side = "left",
        at = c(-0.5, 0, 0.5),
        labels.cex = 0.4,
        tick.length = 0.02,
        labels.niceFacing = TRUE
      )
    }
  }
)

#inner track DMRs
circos.genomicTrack(
  dmrs,           
  ylim = c(-0.3, 0.5), 
  track.height = 0.17,
  track.margin = c(0.05, 0.02),
  bg.border = "black",
  panel.fun = function(region, value, ...) {
    
    val <- value[[1]]
    
    col <- ifelse(
      val > 0,
      "#EF8A62",  # positive values
      "#67A9CF"   # negative values
    )
    
    # Plot the points
    circos.genomicPoints(
      region,
      val,
      pch = 16,
      cex = 0.7,
      col = col
    )
    
    # Draw zero baseline for all chromosomes
    circos.lines(
      x = CELL_META$xlim, 
      y = c(0, 0),        
      col = "black",
      lwd = 0.5
    )
    
    # Draw y-axis only for chromosome 1
    current_chr <- CELL_META$sector.index
    if (current_chr %in% c("1", "chr1")) {
      circos.yaxis(
        side = "left",
        at = seq(-0.5, 0, 0.5), 
        labels.cex = 0.4,
        tick.length = 0.02,
        labels.niceFacing = TRUE
      )
    }
  }
)
dev.off()

save.image()


####make several track views
#order samples according to cluster
meta_ordered = meta.disc[order(meta.disc$cluster_1000_new),]
meta_ordered = factor(meta_ordered$cluster_1000_new, levels = c("METHlow", "METHhigh"))
meta_ordered <- meta_ordered[nrow(meta_ordered):1, ]
all(rownames(meta_ordered) == colnames(beta.disc))
beta.disc = beta.disc[,match(rownames(meta_ordered), colnames(beta.disc))]
all(rownames(meta_ordered) == colnames(beta.disc))


#all protocadherin clusters
cairo_pdf(filename = "Track_view_PCDH_clusters.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',140786136,141512979, beta.disc,
  show.probeNames = FALSE, show.sampleNames = FALSE, heat.height=0.3)
dev.off()

#HOXA
cairo_pdf(filename = "Track_view_HOXA.pdf", width = 8, height = 6)
visualizeRegion(
  'chr7',27090000,27200000, beta.disc,
  show.probeNames = FALSE,show.sampleNames = FALSE, heat.height=0.8)
dev.off()

#HOXD
cairo_pdf(filename = "Track_view_HOXD.pdf", width = 8, height = 6)
visualizeRegion(
  'chr2',176090000,176200000, beta.disc,
  show.probeNames = FALSE,show.sampleNames = FALSE, heat.height=0.8)
dev.off()


#PITX1 gene
cairo_pdf(filename = "Track_view_PITX1.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',135023734,135039000, beta.disc,
  show.probeNames = FALSE)
dev.off()

#PCDHGC3-5
cairo_pdf(filename = "Track_view_PCDHGC3-5.pdf", width = 8, height = 6)
visualizeRegion(
  'chr5',141465000,141515100, beta.disc,
  show.probeNames = FALSE)
dev.off()

#CCND1
cairo_pdf(filename = "Track_view_CCND1.pdf", width = 8, height = 6)
visualizeGene('CCND1', beta.disc,show.probeNames = FALSE, show.sampleNames = FALSE)
dev.off()

#IGF2
cairo_pdf(filename = "Track_view_IGF2.pdf", width = 8, height = 6)
visualizeGene('IGF2', beta.disc,show.probeNames = FALSE, show.sampleNames = FALSE)
dev.off()

save.image()






#######differential methylation UCSF

#get UCSF data, preprocessed
beta.UCSF = readRDS(file = "Betas_UCSF_preprocessed.rds")
meta.UCSF = read.csv(file="meta_UCSF_current.csv")
rownames(meta.UCSF) = meta.UCSF$ID
all(colnames(beta.UCSF)%in%rownames(meta.UCSF))
all(colnames(beta.UCSF)==rownames(meta.UCSF))
meta.UCSF = meta.UCSF[match(colnames(beta.UCSF), rownames(meta.UCSF)),]

# differential methylation in UCSF cohort
#first, make summarizedExperiment
se.UCSF <- SummarizedExperiment(beta.UCSF, colData = meta.UCSF)
colData(se.UCSF)$cluster <- relevel(factor(colData(se.UCSF)$cluster), "METHlow")


#subset summarizedexperiment for individual contrasts not using the entire UCSF cohort
se.UCSF_for_grade_2 = se.UCSF[,!se.UCSF$grade == "3"]
colData(se.UCSF_for_grade_2)$grade <- relevel(factor(colData(se.UCSF_for_grade_2)$grade), "1")

se.UCSF_for_grade_3 = se.UCSF[,!se.UCSF$grade == "2"]
colData(se.UCSF_for_grade_3)$grade <- relevel(factor(colData(se.UCSF_for_grade_3)$grade), "1")


se.UCSF_for_hypermetabolic = se.UCSF[, colData(se.UCSF)$mol_group != "proliferative"]
colData(se.UCSF_for_hypermetabolic)$mol_group <-factor(colData(se.UCSF_for_hypermetabolic)$mol_group)
colData(se.UCSF_for_hypermetabolic)$mol_group <-
  combineLevels(
    colData(se.UCSF_for_hypermetabolic)$mol_group,
    levs = c("Merlinintact", "Immuneenriched"),
    newLabel = "other"
  )
colData(se.UCSF_for_hypermetabolic)$mol_group <- relevel(factor(colData(se.UCSF_for_hypermetabolic)$mol_group), "other")


se.UCSF_for_proliferative = se.UCSF[, colData(se.UCSF)$mol_group != "hypermetabolic"]
colData(se.UCSF_for_proliferative)$mol_group <-factor(colData(se.UCSF_for_proliferative)$mol_group)
colData(se.UCSF_for_proliferative)$mol_group <-
  combineLevels(
    colData(se.UCSF_for_proliferative)$mol_group,
    levs = c("Merlinintact", "Immuneenriched"),
    newLabel = "other"
  )
colData(se.UCSF_for_proliferative)$mol_group <- relevel(factor(colData(se.UCSF_for_proliferative)$mol_group), "other")



###model methylation differences
####for methylation cluster
smry.UCSF.cluster = DML(se.UCSF, ~cluster)
lm_results_UCSF_cluster = summaryExtractTest(smry.UCSF.cluster)

# adjust p values for status
padj.UCSF.cluster = p.adjust(lm_results_UCSF_cluster$Pval_clusterMETHhigh, method = "BH", n = length(lm_results_UCSF_cluster$Pval_clusterMETHhigh))
lm_results_UCSF_cluster = cbind(lm_results_UCSF_cluster, padj.UCSF.cluster)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in METHhigh vs METHlow
hypo_UCSF_cluster = lm_results_UCSF_cluster %>% dplyr::filter(Est_clusterMETHhigh < -0.2 & padj.UCSF.cluster < 0.05)
hyper_UCSF_cluster = lm_results_UCSF_cluster %>% dplyr::filter(Est_clusterMETHhigh > 0.2 & padj.UCSF.cluster < 0.05)
write.csv(hyper_UCSF_cluster, file="probes_hyper_cluster_UCSF.csv")





#model differences grade2 vs grade1
smry.UCSF.grade2vs1 = DML(se.UCSF_for_grade_2, ~grade)
lm_results_UCSF_grade2vs1 = summaryExtractTest(smry.UCSF.grade2vs1)

# adjust p values for grade2 vs grade1
padj.UCSF.grade2 = p.adjust(lm_results_UCSF_grade2vs1$Pval_grade2, method = "BH", n = length(lm_results_UCSF_grade2vs1$Pval_grade2))
lm_results_UCSF_grade2vs1 = cbind(lm_results_UCSF_grade2vs1, padj.UCSF.grade2)
str(lm_results_UCSF_grade2vs1)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade2 over grade1
hypo_UCSF_grade2 = lm_results_UCSF_grade2vs1 %>% dplyr::filter(Est_grade2 < -0.2 & padj.UCSF.grade2 < 0.05)
hyper_UCSF_grade2 = lm_results_UCSF_grade2vs1 %>% dplyr::filter(Est_grade2 > 0.2 & padj.UCSF.grade2 < 0.05)
write.csv(hyper_UCSF_grade2, file="probes_hyper_grade2_UCSF.csv")




#model differences grade3 vs grade1
smry.UCSF.grade3vs1 = DML(se.UCSF_for_grade_3, ~grade)
lm_results_UCSF_grade3vs1 = summaryExtractTest(smry.UCSF.grade3vs1)

# adjust p values for grade3 vs grade1
padj.UCSF.grade3 = p.adjust(lm_results_UCSF_grade3vs1$Pval_grade3, method = "BH", n = length(lm_results_UCSF_grade3vs1$Pval_grade3))
lm_results_UCSF_grade3vs1 = cbind(lm_results_UCSF_grade3vs1, padj.UCSF.grade3)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in grade3 over grade1
hypo_UCSF_grade3 = lm_results_UCSF_grade3vs1 %>% dplyr::filter(Est_grade3 < -0.2 & padj.UCSF.grade3 < 0.05)
hyper_UCSF_grade3 = lm_results_UCSF_grade3vs1 %>% dplyr::filter(Est_grade3 > 0.2 & padj.UCSF.grade3 < 0.05)
write.csv(hyper_UCSF_grade3, file="probes_hyper_grade3_UCSF.csv")




#model differences hypermetabolic vs NF2-WT/Immunogenic
smry.UCSF.hypermetabolic = DML(se.UCSF_for_hypermetabolic, ~mol_group)
lm_results_UCSF_hypermetabolic = summaryExtractTest(smry.UCSF.hypermetabolic)
str(lm_results_UCSF_hypermetabolic)

# adjust p values for hypermetabolic vs reference
padj.UCSF.hypermetabolic = p.adjust(lm_results_UCSF_hypermetabolic$Pval_mol_grouphypermetabolic, method = "BH", n = length(lm_results_UCSF_hypermetabolic$Pval_mol_grouphypermetabolic))
lm_results_UCSF_hypermetabolic = cbind(lm_results_UCSF_hypermetabolic, padj.UCSF.hypermetabolic)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in hypermetabolic vs reference
hypo_UCSF_hypermetabolic = lm_results_UCSF_hypermetabolic %>% dplyr::filter(Est_mol_grouphypermetabolic < -0.2 & padj.UCSF.hypermetabolic < 0.05)
hyper_UCSF_hypermetabolic = lm_results_UCSF_hypermetabolic %>% dplyr::filter(Est_mol_grouphypermetabolic > 0.2 & padj.UCSF.hypermetabolic < 0.05)
write.csv(hyper_UCSF_hypermetabolic, file="probes_hyper_hypermetabolic_UCSF.csv")



#model differences proliferative vs NF2-WT/Immunogenic
smry.UCSF.proliferative = DML(se.UCSF_for_proliferative, ~mol_group)
lm_results_UCSF_proliferative = summaryExtractTest(smry.UCSF.proliferative)

# adjust p values for proliferative vs reference
padj.UCSF.proliferative = p.adjust(lm_results_UCSF_proliferative$Pval_mol_groupproliferative, method = "BH", n = length(lm_results_UCSF_proliferative$Pval_mol_groupproliferative))
lm_results_UCSF_proliferative = cbind(lm_results_UCSF_proliferative, padj.UCSF.proliferative)

#get probe IDs from modeling that show hypo/hyper methylation of 0.2 and padj < 0.05 in proliferative vs reference
hypo_UCSF_proliferative = lm_results_UCSF_proliferative %>% dplyr::filter(Est_mol_groupproliferative < -0.2 & padj.UCSF.proliferative < 0.05)
hyper_UCSF_proliferative = lm_results_UCSF_proliferative %>% dplyr::filter(Est_mol_groupproliferative > 0.2 & padj.UCSF.proliferative < 0.05)
write.csv(hyper_UCSF_proliferative, file="probes_hyper_proliferative_UCSF.csv")


#plot recall for contrasts in discovery
hyper_grade2 = read.csv(file="probes_hyper_grade2_disc.csv")
hyper_grade3 = read.csv(file="probes_hyper_grade3_disc.csv")
hyper_hypermetabolic = read.csv(file="probes_hyper_hypermetabolic_disc.csv")
hyper_proliferative = read.csv(file="probes_hyper_proliferative_disc.csv")

length(intersect(hyper_c2_vs_c1$Probe_ID, hyper_grade2$Probe_ID))/length(hyper_c2_vs_c1$Probe_ID)
length(intersect(hyper_c2_vs_c1$Probe_ID, hyper_grade3$Probe_ID))/length(hyper_c2_vs_c1$Probe_ID)
length(intersect(hyper_c2_vs_c1$Probe_ID, hyper_hypermetabolic$Probe_ID))/length(hyper_c2_vs_c1$Probe_ID)
length(intersect(hyper_c2_vs_c1$Probe_ID, hyper_proliferative$Probe_ID))/length(hyper_c2_vs_c1$Probe_ID)

recall_disc = data.frame(
  condition = c("grade2", "grade3", "hypermetabolic", "proliferative"),
  recall = c(0.002703132,0.7268246,0.07425664,0.978852)
)

recall_disc$condition = factor(recall_disc$condition, levels = c("grade2", "grade3",
                                                                           "hypermetabolic",
                                                                           "proliferative"))

cairo_pdf(filename = "Recall_discovery_hypermethylated_probes.pdf", width = 5.2, height = 7.5)
ggplot(recall_disc, 
       aes(condition, recall, fill=condition))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#8C6BB1","#810F7C","forestgreen","darkorange2")) +
  scale_y_continuous(trans="sqrt")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()


#plot recall for contrasts in UCSF
hyper_c2_vs_c1_EPIC = mLiftOver(hyper_c2_vs_c1$Probe_ID, "EPIC")

length(intersect(hyper_c2_vs_c1_EPIC, hyper_UCSF_grade2$Probe_ID))/length(hyper_c2_vs_c1_EPIC)
length(intersect(hyper_c2_vs_c1_EPIC, hyper_UCSF_grade3$Probe_ID))/length(hyper_c2_vs_c1_EPIC)
length(intersect(hyper_c2_vs_c1_EPIC, hyper_UCSF_hypermetabolic$Probe_ID))/length(hyper_c2_vs_c1_EPIC)
length(intersect(hyper_c2_vs_c1_EPIC, hyper_UCSF_proliferative$Probe_ID))/length(hyper_c2_vs_c1_EPIC)
length(intersect(hyper_c2_vs_c1_EPIC, hyper_UCSF_cluster$Probe_ID))/length(hyper_c2_vs_c1_EPIC)

recall_UCSF = data.frame(
  condition = c("grade2", "grade3", "hypermetabolic", "proliferative", "cluster"),
  recall = c(0.000736377,0.02209131,0.02117084,0.4305965,0.5149116)
)

recall_UCSF$condition = factor(recall_UCSF$condition, levels = c("grade2", "grade3",
                                                                 "hypermetabolic",
                                                                 "proliferative", "cluster"))

cairo_pdf(filename = "Recall_UCSF_hypermethylated_probes.pdf", width = 6.2, height = 7.5)
ggplot(recall_UCSF, 
       aes(condition, recall, fill=condition))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#8C6BB1","#810F7C","forestgreen","darkorange2","violetred4")) +
  scale_y_continuous(trans="sqrt")+
  theme_classic()+
  theme(legend.position = "none")
dev.off()

save.image()











