# scRNA-seq clustering by a range of resolutions, visualization to ensure clustering
# is not determined by uninteresting variables like mitochondrial ratio
# this script should be run after scrna_normalization.R

# https://hbctraining.github.io/scRNA-seq_online/lessons/07_SC_clustering_cells_SCT.html

library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
seurat_phase <- RunUMAP(seurat_phase, 
                             dims = 1:40,
                             reduction = "pca")
#seurat_phase <- readRDS("./results/integrated_seurat.rds")
#heatmap to help determine how many PCs to use when differentiating cell types
DimHeatmap(seurat_phase,nfeatures=20,dims=1:9,cells=500,balanced=TRUE,fast=FALSE)+theme(axis.text.y = element_text(size = 1))
ggsave("plots/PCA_heatmap.png")
#cells is nunmber of cells with the most negative or positive scores to use
## Printing out the most variable genes driving PCs
print(x=seurat_phase[["pca"]],dims=1:10,nfeatures=5)

#can also use elbow plot to identify variation threshold
#ElbowPlot(object=seurat_phase,ndims=40)

#find neighbours, then graph-based clusters
# try a few values for resolution
seurat_phase <- FindNeighbors(object=seurat_phase,dims=1:40)
seurat_phase <- FindClusters(object=seurat_phase,resolution=c(0.4,0.6,0.8,1,1.2,1.4))

seurat_phase@meta.data%>%View()
for (val in c(0.3,0.4,0.6,0.8,1,1.2,1.4)){
  Idents(object=seurat_phase) <- paste("RNA_snn_res.",as.character(val),sep="")
  print(DimPlot(seurat_phase,reduction='umap',group.by='ident',label=TRUE,label.size=6)+ggtitle(paste("Resolution=",as.character(val),sep="")))
  ggsave(paste("plots/clustering_resolution=",val,".png"))
  }
#The clusters seem very consistent for a wide range of resolutions beyond 0.4


#https://hbctraining.github.io/scRNA-seq_online/lessons/08_SC_clustering_quality_control.html
#Determine cells per cluster
FetchData(seurat_phase,vars=c("ident","orig.ident"))%>%
  dplyr::count(ident,orig.ident)%>%
  tidyr::spread(ident,n)


#determine whether they cluster by cell cycle stage
DimPlot(seurat_phase,label=TRUE,split.by="Phase")+NoLegend() #all about the same unsurprisingly
ggsave("plots/clustering_by_cell_cycle.png")
#no clustering by cell cycle stage, thankfully

#check if they cluster by other uninteresting sources of variation.
metrics <- c("nUMI","nGene","S.Score","G2M.Score","mitoRatio")
print(FeaturePlot(seurat_phase,reduction="umap",features=metrics,pt.size=0.4,
            order=T,min.cutoff="q10",label=T))
ggsave("plots/clustering_by_other_variables.png")
#minimal clustering by these other variables as well

# https://hbctraining.github.io/scRNA-seq_online/lessons/08_SC_clustering_quality_control.html has list of some cell types and markers
