# http://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html
library(celldex)
library(SingleR)
ref <- BlueprintEncodeData()
ref

pred <- SingleR(test = seurat_phase@assays$RNA@data, ref=ref, labels=ref$label.main)
table(pred$labels)
pred$labels

seurat_phase@meta.data$predicted_celltype <- pred$labels
FeaturePlot(seurat_phase,reduction="umap",features=c("predicted_celltype"),
            order=T,min.cutoff='q10',label=T)

DimPlot(seurat_phase,reduction='umap',group.by='predicted_celltype',label=TRUE,label.size=6)+ggtitle("Cell type Identification by SingleR")
ggsave("plots/automatic_celltype_identification.png")

Idents(object=seurat_phase) <- paste("RNA_snn_res.",as.character(0.4),sep="")
DimPlot(seurat_phase,reduction='umap',group.by='ident',label=TRUE,label.size=6)+ggtitle("Resolution=0.4 clustering")


temp <-plotScoreHeatmap(pred)#export manually from RStudio
ggsave("plots/celltype_identification_heatmap.png",temp)
tab <- table(Assigned=seurat_phase@meta.data$predicted_celltype, Cluster=seurat_phase@meta.data$RNA_snn_res.0.4)

# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
temp <- pheatmap(log2(tab+10), angle_col="0",fontsize_row=15,fontsize_col=15,color=colorRampPalette(c("white", "blue"))(101))
ggsave("plots/celltype_cluster_pheatmap.png",temp)

seurat_phase@meta.data$annotated_celltype <- NA
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.0.4 %in% c(1)]="B cells"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.0.4 %in% c(2)]="Monocytes"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.0.4 ==6]="NK cells"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.0.4 %in% c(7,4,5)]="CD8+ T cells"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.0.4 %in% c(0,3)]="CD4+ T cells"

DimPlot(seurat_phase,reduction='umap',group.by='annotated_celltype',label=TRUE,label.size=6)+ggtitle("Cluster cell type Annotation")
ggsave("plots/annotated_cluster_cell_types.png")

# can also use a  custom reference, see http://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html 1/3 of the way downn