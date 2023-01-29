# reading in Seurat object for 1000 PBMC scRNA-seq data, then QC to remove low-
# quality cells, debris and doublets etc.

#https://hbctraining.github.io/scRNA-seq_online/lessons/03_SC_quality_control-setup.html

library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

#read count in as a sparse matrix as there are lots of 0s in this datatype
# Create a Seurat object for each sample

seurat_data <- Read10X(data.dir = paste0("data/raw_feature_bc_matrix"))
pbmcs <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 200)#min.features is min genes detected to not be filtered out as a low quality cell

head(pbmcs)
view(pbmcs@meta.data) #each row is a cell (barcode), with specified number of RNAs and unique RNAs/UMIs

#https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html


#calculate novelty score-number of genes per UMI
pbmcs$log10GenesPerUMI <- log10(pbmcs$nFeature_RNA)/log10(pbmcs$nCount_RNA)
#determine proportion of trascripts that map to mitochondrial genes - should be low for healthy cells
pbmcs$mitoRatio <- PercentageFeatureSet(object=pbmcs,pattern="^MT-") #pattern may be different for nonhuman organisms. Start by trying MT, Mt, mt
pbmcs$mitoRatio <- pbmcs@meta.data$mitoRatio/100
metadata <- pbmcs@meta.data
metadata$cells <- rownames(metadata)

#rename columns
metadata <- metadata%>%dplyr::rename(seq_folder=orig.ident,nUMI=nCount_RNA,
                                     nGene=nFeature_RNA)
pbmcs@meta.data <- metadata
save(pbmcs,file="data/filtered_seurat.RData")


print(paste(nrow(pbmcs@meta.data),"cells of 1000 cells expected.")) #1215 cells after filtering
#now use visualizations to decide which cells to remove

dir.create("plots")
#UMIs/transcripts per cell - should be above 500
metadata%>%ggplot(aes(x=nUMI))+geom_density(alpha=0.2)+
  scale_x_log10()+theme_classic()+ylab("Cell Density")+geom_vline(xintercept=3000)+annotate(geom="text",x=4000,y=1,label="x=3000")+
  theme(plot.title=element_text(hjust=0.5,face="bold"))+ggtitle(paste("UMIs per",nrow(pbmcs@meta.data), "Cells"))
ggsave("plots/initial_UMI_density.png")

#genes per cell.  should be one peak, usually not bimodal with a shoulder (indicates failure or very different types/proportions of cells)
metadata%>%ggplot(aes(x=nGene))+geom_density(alpha=0.2)+
  scale_x_log10()+theme_classic()+ylab("Cell Density")+geom_vline(xintercept=2000)+annotate(geom="text",x=2500,y=1,label="x=2000")+
  theme(plot.title=element_text(hjust=0.5,face="bold"))+ggtitle(paste("Genes per",nrow(pbmcs@meta.data), "Cells"))+geom_vline(xintercept=4300)
ggsave("plots/initial_gene_density.png")

#novelty score(genes per UMI).  .8 is expected for good quality cells
metadata%>%ggplot(aes(x=log10GenesPerUMI))+geom_density(alpha=0.2)+
  theme_classic()+ylab("Cell Density")+geom_vline(xintercept=0.82)+annotate(geom="text",x=.83,y=5,label="x=82")+
  theme(plot.title=element_text(hjust=0.5,face="bold"))+ggtitle("Genes per UMI")
ggsave("plots/initial_genes_per_UMI.png")

#mitochondrial ratio, should be below .2 usually

ggplot(metadata,aes(x=mitoRatio)) + 
geom_density(alpha = 0.2) + 
scale_x_log10() + 
theme_classic() +
geom_vline(xintercept = 0.2)+annotate(geom="text",x=0.25,y=1,label="x=0.2")+
theme(plot.title=element_text(hjust=0.5,face="bold"))+ggtitle(paste("Mitochondrial Ratio per",nrow(pbmcs@meta.data), "Cells"))
ggsave("plots/initial_mitochondrial_ratio.png")
  
ggplot(metadata,aes(x=nUMI,y=nGene,color=mitoRatio))+
  geom_point()+scale_colour_gradient(low="gray90",high="black")+
  stat_smooth(method=lm)+scale_x_log10()+scale_y_log10()+theme_classic()+
  geom_vline(xintercept=500)+geom_hline(yintercept=250)+
  theme(plot.title=element_text(hjust=0.5,face="bold"))+ggtitle("Initial nGenes vs nUMIs")
ggsave("plots/initial_nGenes_vs_nUMIs.png")

#bottom right cells below blue line could be dying or low complexity like RBCs
#bottom left are likely poor quality cells
#note mitochondrial reads are only a significant proportion in low count cells

#filter - adjust these thresholds according to your application
#it is best to consider all metrics, not just one of above in isolation, to avoid excluding viable cells
filtered_seurat <- subset(x=pbmcs,subset=(nUMI>=3000)&(nGene<=4200)&
                            (nGene>=2000)&(log10GenesPerUMI>0.82)&
                            (mitoRatio<0.2))

print(paste(nrow(filtered_seurat@meta.data),"cells of 1000 cells expected."))
#we are down to 1094 cells - we could be more stringent with our thresholds, but 
# the remaining cells do all seem good quality, and based on my experimentation
# this number could only be reduced more by excluding likely viable cells
#remove genes with 0 counts so they don't reduce average expression
counts <- GetAssayData(object=filtered_seurat,slot="counts")
nonzero <- counts>0
#remove genes expressed in 10 or fewer cells
keep_genes <- Matrix::rowSums(nonzero)>=10
filtered_counts <- counts[keep_genes,]
filtered_seurat <- CreateSeuratObject(filtered_counts,meta.data=filtered_seurat@meta.data)
save(filtered_seurat, file="data/seurat_filtered.RData")
