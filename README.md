# 1000_PBMCs
 
Sample analysis of single-cell RNA sequencing data of 1000 peripheral blood mononuclear cells available from [10X Genomics](https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0).  Sequence reads were mapped to a reference transcriptome using cellranger count before quality control and clustering were performed.  This workflow generally follows the Harvard Chan Bioinformatics Core workshop available [here](https://hbctraining.github.io/scRNA-seq_online/schedule/links-to-lessons.html).

The intended script order is quality_control.R --> scrna_normalization.R --> clustering.R --> cluster_identification.R
