library(sesame)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(FactoMineR)
library(factoextra)
library(dplyr)



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







save.image(file = "my_work_space.RData")



