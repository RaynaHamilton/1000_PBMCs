# Normalize scRNA-seq Seurat object and check for uninteresting variables to regress out if necessary
# this file should be run after quality_control.R

#https://hbctraining.github.io/scRNA-seq_online/lessons/06_SC_SCT_normalization.html
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(AnnotationHub)
#filtered_seurat <- load("data/seurat_filtered.RData")
load("data/cycle.rda")

seurat_phase <- NormalizeData(filtered_seurat) #log normalization
seurat_phase <- CellCycleScoring(seurat_phase,g2m.features=g2m_genes,
                                 s.features=s_genes)
#note # of PCs=# of cells
#before PCA we must scale the data since by default the most expressed genes will
#be the most variable otherwise
#find most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase,selection.method="vst",
                                     nfeatures=2000,verbose=FALSE)
#scale counts
seurat_phase <- ScaleData(seurat_phase)

#now PCA...
seurat_phase <- RunPCA(seurat_phase)
DimPlot(seurat_phase,reduction="pca",group.by="Phase",split.by="Phase")+ggtitle("PCA plot by Cell Cycle Phase")
# The PCA biplots look essentially identical regardless of cell cycle phase, so
# we will not regress this out
ggsave("plots/cell_cycle_PCA.png")

#deciding whether to regress out mitochondrial expression...
summary(seurat_phase@meta.data$mitoRatio)
#now cut based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio,
                                     breaks=c(-Inf,as.numeric(summary(seurat_phase@meta.data$mitoRatio))[c(2,3,5)],Inf),
                                     labels=c("Low","Medium","Medium high","High"))
DimPlot(seurat_phase,reduction="pca",group.by="mitoFr",split.by="mitoFr")+ggtitle("PCA plot by Mitochondrial Ratio")

ggsave("plots/mitoRatio_PCA.png")

#Regression does not seem necessary for either variable-you would use the below line to regress
#seurat_phase <- SCTransform(seurat_phase, vars.to.regress = c("mitoRatio"))
seurat_phase@assays#top 10 variable features
saveRDS(seurat_phase, "data/seurat_post_PCA.rds")
#to load use split_seurat <- readRDS("data/split_seurat_regressed.rds")

#https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html
#integration follows normalization if you have multiple samples to compare
