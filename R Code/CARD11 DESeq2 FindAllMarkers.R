library(viridis)
library(colorspace)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice) 
library(dplyr)
library(viridis)
library(ggokabeito)
library(scCustomize)
library(colorBlindness)
library(Seurat)
mem.maxVSize(1e10)
############################
## CARD11 : Perform DEA across pseudobulked unsupervised clusters from res 0.15
############################

setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/11725-JM")

B_cells <- readRDS("B_cells.rds")
B_cells$newer.ident <- NA
B_cells$newer.ident[which(B_cells$orig.ident == "0008")] <- "HC"
B_cells$newer.ident[which(B_cells$orig.ident == "0009")] <- "HC"
B_cells$newer.ident[which(B_cells$orig.ident == "0015")] <- "Homozygote"
B_cells$newer.ident[which(B_cells$orig.ident == "0016")] <- "Heterozygote"
B_cells$newer.ident <- factor(x = B_cells$newer.ident, levels = c("HC", "Heterozygote", "Homozygote"))
B_cells$new.ident <- factor(x = B_cells$new.ident, levels = c("HC 1", "HC 2", "WT/G126D", "R331P/R331P"))

# Pseudobulk by cluster (res 0.15) and sample
B_cells$cluster_genotype <- paste( B_cells$newer.ident,B_cells$RNA_snn_res.0.15, sep="-" )

pseudo_B <- AggregateExpression(B_cells, assays = "RNA", group.by = c('new.ident','RNA_snn_res.0.15'), return.seurat = T)
Cells(pseudo_B)
Idents(pseudo_B)<- 'RNA_snn_res.0.15'

Deseq2Output <- FindAllMarkers(pseudo_B, assay = "RNA", min.cells.group = 1, logfc.threshold = 0, test.use = 'DESeq2')

gc()

#Extract genes from CARD11 patient cluster
cluster1 <- DEA[DEA$cluster =='1' & DEA$p_val_adj < 0.05 ,]


# DESeq2 output did not give reliable results possibly because we are averaging each healthy donor, heterozygote and homozygote in cluster 1.


## Perform DEA through the wilcoxon rank-sum test.
WilcoxOutput <- FindAllMarkers(B_cells, group.by = 'RNA_snn_res.0.15', assay = "RNA", logfc.threshold = 0, min.pct = 0.01)

#Extract genes from CARD11 patient cluster.
cluster1 <- WilcoxOutput[WilcoxOutput$cluster =='1' & WilcoxOutput$p_val_adj < 0.05 & WilcoxOutput$avg_log2FC > 0,]
cluster1_GSEA <- WilcoxOutput[WilcoxOutput$cluster =='1',]
FeaturePlot(B_cells, features = c('IGF2BP3','ENSG00000289474','ARHGAP24','MAK','ESYT2'))
# DEGs identified through the Wilcoxon rank sum test seem to be more accurate. 

setwd("/Volumes/pmi/Data/Research/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/SupplementaryData")
write.csv(WilcoxOutput, 'FindAllMarkers output wilcoxon rank-sum test.csv')
write.csv(cluster1, 'Cluster 1 homozygote cluster ORA.csv')
write.csv(cluster1_GSEA, 'Cluster 1 homozygote cluster GSEA.csv')



