# cluster identification of the small 1000-PBMCs dataset using celldex and SingleR
# this script should be run after clustering.R

# http://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html
library(celldex)
library(SingleR)
ref <- BlueprintEncodeData()

# general, broad cell type labels
pred <- SingleR(test = seurat_phase@assays$RNA@data, ref=ref, labels=ref$label.main)#label.main for broad cell types, label.fine and label.ont are more detailed
table(pred$labels)#B and T cells, NK and monocytes - all expected for PBMCs
seurat_phase@meta.data$predicted_celltype <- pred$labels
DimPlot(seurat_phase,reduction='umap',group.by='predicted_celltype',label=F)+ggtitle("Main cell type Identification by SingleR")
ggsave("plots/main_automatic_celltype_identification.png")

#fine cell type identification
pred <- SingleR(test = seurat_phase@assays$RNA@data, ref=ref, labels=ref$label.fine)
table(pred$labels) # now we can see naive vs memory B cells, Tregs etc
seurat_phase@meta.data$predicted_celltype <- pred$labels
DimPlot(seurat_phase,reduction='umap',group.by='predicted_celltype',label=F)+ggtitle("Fine cell type Identification by SingleR")
ggsave("plots/fine_automatic_celltype_identification.png")

#reminder of what the clustering looked like:
Idents(object=seurat_phase) <- paste("RNA_snn_res.",as.character(1.4),sep="")
DimPlot(seurat_phase,reduction='umap',group.by='ident',label=TRUE,label.size=6)+ggtitle("Resolution=1.4 clustering")

#score heatmap by cell type
temp <-plotScoreHeatmap(pred)
ggsave("plots/celltype_identification_heatmap.png",temp)
tab <- table(Assigned=seurat_phase@meta.data$predicted_celltype, Cluster=seurat_phase@meta.data$RNA_snn_res.1.4)

# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
temp <- pheatmap(log2(tab+10), angle_col="0",fontsize_row=15,fontsize_col=15,color=colorRampPalette(c("white", "blue"))(101))
ggsave("plots/celltype_cluster_pheatmap.png",temp)
#cluster 8 unnfortunately has a mix of Cd4+ and CD8+ cells (perhaps intermediate cell type or not fully differentiated)
# but the remaining clusters have quite clear identities, though the data or resolution does not seem sufficient to fully distinguish
# between subtypes of CD4+ and CD8+

seurat_phase@meta.data$annotated_celltype <- NA
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.1.4 %in% c(0,2)]="CD4+ T cells"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.1.4 %in% c(1)]="Monocytes"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.1.4 %in% c(3)]="Naive B cells"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.1.4 %in% c(4,5,8)]="CD8+ T cells"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.1.4 %in% c(6)]="NK cells"
seurat_phase@meta.data$annotated_celltype[seurat_phase@meta.data$RNA_snn_res.1.4 %in% c(7)]="Memory B cells"

DimPlot(seurat_phase,reduction='umap',group.by='annotated_celltype',label=TRUE,label.size=6)+ggtitle("Cluster cell type Annotation")
ggsave("plots/annotated_cluster_cell_types.png")

# can also use a  custom reference, see http://bioconductor.org/books/3.13/OSCA.basic/cell-type-annotation.html 1/3 of the way downn