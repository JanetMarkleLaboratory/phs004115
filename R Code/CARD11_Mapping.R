
#############################################################
# Samples 11725-JM 0008, 0009, 0015, and 0016 will be analyzed for T cell anergy and B cell expansion related to CARD11 mutation. 
#############################################################

#############################################################
# PART 1: Quality filtering and doublet removal
#############################################################

## Install and load packages
### Seurat v5.1.0 does not recognize the counts layer in the SCT assay. Install v5.0.3 
remotes::install_version("SeuratObject", "5.0.1", repos = c("https://satijalab.r-universe.dev", getOption("repos")), force = TRUE)
remotes::install_version("Seurat", "5.0.3", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scDblFinder")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocParallel")


library(Seurat)
library(SeuratObject)
library(cowplot)
library(ggplot2)
library(BiocManager)
library(tidyverse)
library(dplyr)
library(scDblFinder)
library(BiocParallel)
options(MulticoreParam=MulticoreParam(workers=4))
library(hdf5r)

## Load in the data for CARD11 heterozygote, CARD11 homozygote and two healthy controls. 

setwd("/Volumes/pmi/Data/Research/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/11725-JM")
samples <- c('count_11725-JM-0008', 'count_11725-JM-0009', 'count_11725-JM-0015', 'count_11725-JM-0016')

d10x.data <- sapply(samples, function(i){
  d10x <- Read10X_h5(file.path(i, 'filtered_feature_bc_matrix.h5'))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x), split = "-"), '[[',1L),i,sep = "-")
  d10x})

card.data <- do.call("cbind", d10x.data)

merged <- CreateSeuratObject(
  card.data,
  project = "count_11725-JM",
  min.cells = 3,
  min.features = 10,
  names.field = 2,
  names.delim = "\\-JM-"
)


## Add mito and ribo percentage. Create scatter plots to access quality of cells by sample.

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")

FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  NoLegend() + ggtitle('Number of UMIs per cell vs. Number of Genes per cell') + 
  scale_y_log10() + scale_x_log10() + facet_wrap(~merged$orig.ident)


FeatureScatter(merged, feature1 = "percent.mt", feature2 = "percent.ribo") +
  NoLegend() + ggtitle('Percent mito per cell vs. Percent ribo per cell') + 
  scale_y_log10() + scale_x_log10() + facet_wrap(~merged$orig.ident)



##Split the object for sample-specific quality filtering

#Create a seurat object list
splitobj <- SplitObject(merged, split.by = "orig.ident")
rm(merged)
rm(d10x.data)
   
#Create a list for faster processing with a for loop
ids <- c('control8', 'control9', 'patient15', 'patient16')

#Create new ids
splitobj$control8 <- splitobj$'0008'
splitobj$control9 <- splitobj$'0009'
splitobj$patient15 <- splitobj$'0015'
splitobj$patient16 <- splitobj$'0016'

#Remove old IDs
splitobj$'0008' <- NULL
splitobj$'0009' <- NULL
splitobj$'0015' <- NULL
splitobj$'0016' <- NULL

## Remove low quality cells with a threshold cutoff of the lowest 5% of features 

minCov=1000 #If the minimum nCount is less than 1000, make the cutoff threshold greater than 5% of the quantile distribution. 
#Make a countHIGH variable that has a cutoff lower than 95% of the quantile distribution to remove empty drops, dead cells, or multiplets. Remove the lower quantile of nFeature_RNA which indicate dead or dying cells. 

for (i in ids){
  if(min(splitobj[[i]]$nCount_RNA)>=minCov){
    countLOW=min(splitobj[[i]]$nCount_RNA)
  }else{
    countLOW=quantile(splitobj[[i]]$nCount_RNA, prob=c(0.05))  
  }
  countHIGH=quantile(splitobj[[i]]$nCount_RNA, prob=0.95)
  featureLOW=quantile(splitobj[[i]]$nFeature_RNA, prob=0.05)
  
  ##subset
  splitobj[[i]] <- subset(splitobj[[i]], subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt <= 10)
}
VlnPlot(splitobj$control8, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$control9, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$patient15, features = c('nFeature_RNA', 'nCount_RNA'))
VlnPlot(splitobj$patient16, features = c('nFeature_RNA', 'nCount_RNA'))


## Remove doublets using scDblFinder

for (i in ids){
  sce <- as.SingleCellExperiment(splitobj[[i]])
  
  sce <- scDblFinder(sce, clusters = TRUE, BPPARAM = MulticoreParam(4), samples = "orig.ident")
  
  print(table(sce$orig.ident, sce$scDblFinder.class))
  
  all.equal(Cells(splitobj[[i]]), colnames(sce))
  
  splitobj[[i]] <- AddMetaData(splitobj[[i]], sce$scDblFinder.class, col.name= 'scDblFinder.class')
  
  rm(sce) ; gc()
  
  splitobj[[i]] <- subset(splitobj[[i]], subset = scDblFinder.class == 'singlet')
}

# Make certain that the doublets identified are removed from the object.
table(splitobj$control8$scDblFinder.class)
table(splitobj$control9$scDblFinder.class)
table(splitobj$patient15$scDblFinder.class)
table(splitobj$patient16$scDblFinder.class)


## Save RDS object for later analyses and then clear your global environment to free up RAM.

saveRDS(splitobj,'split_CARD11.rds')

#############################################################
# PART 2: Post-filtering Seurat Workflow: Mapping to a PBMC dataset
#############################################################

## Read in RDS

splitobj <- readRDS(file='11725-JM/split_CARD11.rds')

## Read in cell cycle scoring table for cell cycle regression

exp.mat <- read.table(file = "11725-JM/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

## Perform Seurat workflow

options(future.globals.maxSize = 20000 * 1024^2)
options(future.rng.onMisuse = 'ignore')

splitobj <- lapply(X = splitobj, FUN=NormalizeData)

## Perform cell cycle regression and SCTransform for mapping to a reference
for (i in ids){
  splitobj[[i]] <- JoinLayers(splitobj[[i]])
  splitobj[[i]] <- CellCycleScoring(splitobj[[i]], s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  splitobj[[i]] <- SCTransform(splitobj[[i]], vars.to.regress = c("S.Score", "G2M.Score"))
}

## Read in the reference PBMC data

reference <- readRDS("11725-JM/pbmc_multimodal_2023.rds")
DimPlot(object = reference, reduction = "wnn.umap", group.by = "celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

anchors <- list()
for (i in 1:length(splitobj)) {
  anchors[[i]] <- FindTransferAnchors(
  reference = reference,
  query = splitobj[[i]],
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
  )}

## analyze % of query cells with anchors (# of anchors / # of cells *100)
control1 <- anchors[[1]]@anchors[,"cell2"] %>% unique() %>% length() / 16761 *100
control2 <- anchors[[2]]@anchors[,"cell2"] %>% unique() %>% length() / 22291 *100
patient15 <- anchors[[3]]@anchors[,"cell2"] %>% unique() %>% length() / 4397 *100
patient16 <- anchors[[4]]@anchors[,"cell2"] %>% unique() %>% length() / 14470 *100
control1
#13.7939263767078
control2
#10.8115382889956
patient15
#7.550603
patient16
#10.93988


## lower than 15% mapping suggests that either the reference has a slightly different biological state than the query samples or that there is batch effect. Since, we do not have biological replicates,
## there is no reasonable way to correct for any possible batch effects. Notably, homozygote patient 15 has a lower mapping score.

## Transfer cell type labels from the reference to the query. Project the query data onto the UMAP structure of the reference.
for (i in 1:length(splitobj)){
  splitobj[[i]] <- MapQuery(
    anchorset = anchors[[i]],
    query = splitobj[[i]],
    reference = reference,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca", 
    reduction.model = "wnn.umap"
  )
}

##
## Compute mapping score and merge the list
for (i in 1:length(splitobj)){splitobj[[i]] <- AddMetaData(splitobj[[i]], metadata = MappingScore(anchors[[i]]), col.name = 'mappingscore')}
rm(reference)


object <- merge(splitobj[[1]], splitobj[2:length(splitobj)], merge.dr = "ref.umap")


p1 = DimPlot(object, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE,split.by = "orig.ident") + NoLegend()
p2 = DimPlot(object, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE,split.by = "orig.ident") + NoLegend()
p1 + p2

Idents(object) <- 'predicted.celltype.l1'
DefaultAssay(object) <- 'prediction.score.celltype.l1'
FeaturePlot(object, features = c("B"), ncol = 3, 
            cols = c("lightgrey", "darkred"), split.by = "orig.ident")
DefaultAssay(object) = "SCT"
VlnPlot(object, features = c("CD79A", "RALGPS2", "CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", "IGHM", "MEF2C"), sort = TRUE) + NoLegend()

object$new.ident <- NA
object@meta.data$new.ident[which(object@meta.data$orig.ident == "0008")] <- "Control 1"
object@meta.data$new.ident[which(object@meta.data$orig.ident == "0009")] <- "Control 2"
object@meta.data$new.ident[which(object@meta.data$orig.ident == "0015")] <- "Homozygote"
object@meta.data$new.ident[which(object@meta.data$orig.ident == "0016")] <- "Heterozygote"
object$celltype.genotype <- paste(Idents(object), object$new.ident, sep = "_")

#############################################################
## PART 3: Filtering mapped output and predictive scores
### analyze predicted scores by creating a boxplot and filter accordingly
#############################################################

#### Find the celltype level 1 filtering threshold

df <- data.frame(object$predicted.celltype.l1)
df <- tibble::rownames_to_column(df, "predicted.celltype.l1")
score <- data.frame(object$predicted.celltype.l1.score)
score <- tibble::rownames_to_column(score, "predicted.celltype.l1")
sample <- data.frame(object$new.ident)
sample <- tibble::rownames_to_column(sample, "predicted.celltype.l1")
df <- merge(score, df, by = "predicted.celltype.l1")
df<- merge (df, sample, by = "predicted.celltype.l1")

df %>%
  ggplot( aes(x=object.new.ident, y=object.predicted.celltype.l1.score, fill=object.new.ident)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.25, alpha=0.25) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + ylab("Cell Type L1 Prediction Score Across All Cells")

## Let's look at only B cell prediction scores
B_df <- subset(df, object.predicted.celltype.l1 == "B")
B_df %>%
  ggplot( aes(x=object.new.ident, y=object.predicted.celltype.l1.score, fill=object.new.ident)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.25, alpha=0.25) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + ylab("Cell Type L1 Prediction Score Across All B Cells")

## Remove cells one std dev below the mean prediction score. NOTE: The cutoff threshold is 0.7612352
mean_cells <- mean(df$object.predicted.celltype.l1.score)
stdev <- sd(df$object.predicted.celltype.l1.score)
cutoff_thresh <- mean_cells - stdev

hom<- subset(df, object.new.ident == "Homozygote")
homb <- subset(hom, object.predicted.celltype.l1 =="B")
mean_cells_hom <- mean(hom$object.predicted.celltype.l1.score)
#0.9400908
mean_cells_homb <- mean(homb$object.predicted.celltype.l1.score)
#0.9680142
## The mean B cell predicted score for Homozygote patient is above the overall sample mean at 0.9680142. Therefore, I think it's fine to remove low predictive scores although homozygote
## B cells may be in a different biological state.

df_filtered <- subset(df, object.predicted.celltype.l1.score >= cutoff_thresh)
df_filtered %>%
  ggplot( aes(x=object.new.ident, y=object.predicted.celltype.l1.score, fill=object.new.ident)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.25, alpha=0.25) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + ylab("Cell Type Prediction Score Across All Cells Post-Filtering")
B_df_filtered <- subset(df_filtered, object.predicted.celltype.l1 == "B")
table(B_df$object.predicted.celltype.l1, B_df$object.new.ident)
table(B_df_filtered$object.predicted.celltype.l1, B_df_filtered$object.new.ident)

# > table(B_df$object.predicted.celltype.l1, B_df$object.new.ident)
#Control 1 Control 2 Heterozygote Homozygote
#B      2206      3190         9984       4040

#> table(B_df_filtered$object.predicted.celltype.l1, B_df_filtered$object.new.ident)
#Control 1 Control 2 Heterozygote Homozygote
#B      2188      3172         9914       3841

# calculate percentage of change between filtered and unfiltered B cells. About 5% of predicted B cells were removed from the homozygote dataset (-4.925743%).
control1change <- (2188-2206) /2206 *100
#-0.8159565
control2change <- (3172-3190)/ 3190 *100
#-0.5642633
hetchange <- (9914-9984) /9984 *100
#-0.7011218
homchange <- (3841-4040)/4040 *100
#-4.925743

#### predicted level 2 boxplots

df <- data.frame(object$predicted.celltype.l2)
df <- tibble::rownames_to_column(df, "predicted.celltype.l2")
score <- data.frame(object$predicted.celltype.l2.score)
score <- tibble::rownames_to_column(score, "predicted.celltype.l2")
sample <- data.frame(object$new.ident)
sample <- tibble::rownames_to_column(sample, "predicted.celltype.l2")
df <- merge(score, df, by = "predicted.celltype.l2")
df<- merge (df, sample, by = "predicted.celltype.l2")

df %>%
  ggplot( aes(x=object.new.ident, y=object.predicted.celltype.l2.score, fill=object.new.ident)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.25, alpha=0.25) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + ylab("Cell Type Prediction L2 Score Across All Cells")

## Let's look at only B cell prediction scores
B_df_naive <- subset(df, object.predicted.celltype.l2 == c("B naive"))
B_df_interm <- subset(df, object.predicted.celltype.l2 == c("B intermediate"))
B_df_memory <- subset(df, object.predicted.celltype.l2 == c("B memory"))
B_df <- rbind(B_df_memory, B_df_interm, B_df_naive)

B_df %>%
  ggplot( aes(x=object.new.ident, y=object.predicted.celltype.l2.score, fill=object.predicted.celltype.l2)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.25, alpha=0.25) +
  theme(
    plot.title = element_text(size=11)
  ) +
  xlab("") + ylab("Cell Type L2 Prediction Score Across All B Cells")

## Remove cells one std dev below the mean l2 prediction score. NOTE: The cutoff threshold is 0.5200571
mean_cells <- mean(df$object.predicted.celltype.l2.score)
stdev <- sd(df$object.predicted.celltype.l2.score)
cutoff_thresh <- mean_cells - stdev


rm(splitobj)
p1 = DimPlot(object, reduction = "ref.umap", group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, repel = TRUE,split.by = "orig.ident") + NoLegend()
p2 = DimPlot(object, reduction = "ref.umap", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3 ,repel = TRUE,split.by = "orig.ident") + NoLegend()
p1 + p2

saveRDS(object, 'prefiltered_object.rds')

object.sct.filtered <- subset(object, predicted.celltype.l1.score >= 0.7612352)

table(object.sct.filtered$orig.ident)

#0008  0009  0015  0016 
#14098 18088  3949 13287 
table(object$orig.ident)

#0008  0009  0015  0016 
#16761 22291  4397 14470 
Idents(object.sct.filtered) <- "predicted.celltype.l1"
VlnPlot(object.sct.filtered, features = c("CD79A", "RALGPS2", "CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", "IGHM", "MEF2C"), sort = TRUE) + NoLegend()

#Display a dotplot by cell type

DotPlot(object.sct.filtered, features = c("CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "NCAM1", "CD56", "KLRB1", "NKG7", "CD19", "MS4A1", "CD20", "CD38",
                                          "CD14", "CD68", "FCGR3A", "CD16", "CD1C", "LILRA44", "CD34","HBB", "PPBP"),
        cluster.idents = T, assay = "SCT")

################################
## Differential Expression Analysis
################################
DefaultAssay(object.sct.filtered) <- "SCT"

## Upregulated in ident.1 == higher log2FC and upregulated in ident.2 == lower log2FC
object.sct.filtered <- PrepSCTFindMarkers(object.sct.filtered)
Idents(object.sct.filtered) <- 'celltype.genotype'

## Save the seurat object
saveRDS(object.sct.filtered, 'object.sct.filtered.rds')

object.sct.filtered<- readRDS('object.sct.filtered.rds')
#########################
## Run FindMarkers to create gene lists relevant to BENTA.
#########################

b.cont1hom <- FindMarkers(object.sct.filtered, ident.2 = "B_Control 1", ident.1 = "B_Homozygote",  assay = "SCT")
b.cont2hom <- FindMarkers(object.sct.filtered, ident.2 = "B_Control 2", ident.1 = "B_Homozygote",  assay = "SCT")
b.controls <- FindMarkers(object.sct.filtered, ident.2 = "B_Control 1", ident.1 = "B_Control 2",  assay = "SCT")
b.cont1het <- FindMarkers(object.sct.filtered, ident.2 = "B_Control 1", ident.1 = "B_Heterozygote", assay = "SCT")
b.cont2het <- FindMarkers(object.sct.filtered, ident.2 = "B_Control 2", ident.1 = "B_Heterozygote", assay = "SCT")
b.homhet <- FindMarkers(object.sct.filtered, ident.2 = "B_Homozygote", ident.1 = "B_Heterozygote", assay = "SCT")

t.cont1hom <- FindMarkers(object.sct.filtered, ident.2 = "CD4 T_Control 1", ident.1 = "CD4 T_Homozygote",assay = "SCT")
t.cont2hom <- FindMarkers(object.sct.filtered, ident.2 = "CD4 T_Control 2", ident.1 = "CD4 T_Homozygote",assay = "SCT")
t.controls <- FindMarkers(object.sct.filtered, ident.2 = "CD4 T_Control 1", ident.1 = "CD4 T_Control 2", assay = "SCT")
t.cont1het <- FindMarkers(object.sct.filtered, ident.2 = "CD4 T_Control 1", ident.1 = "CD4 T_Heterozygote", assay = "SCT")
t.cont2het <- FindMarkers(object.sct.filtered, ident.2 = "CD4 T_Control 2", ident.1 = "CD4 T_Heterozygote",assay = "SCT")
t.homhet <- FindMarkers(object.sct.filtered, ident.2 = "CD4 T_Homozygote", ident.1 = "CD4 T_Heterozygote", assay = "SCT")

write.csv(b.homhet, '11725-JM/b.homhet.csv')
write.csv(b.cont1hom, '11725-JM/b.cont1hom.csv')
write.csv(b.cont2hom, '11725-JM/b.cont2hom.csv')
write.csv(b.cont1het, '11725-JM/b.cont1het.csv')
write.csv(b.cont2het, '11725-JM/b.cont2het.csv')
write.csv(b.controls, '11725-JM/b.controls.csv')

library(viridis)
pal <- viridis(n = 10, option = "C", direction = -1)
# Plot B cell markers
FeaturePlot(B_cells, features = c("CD79A", "RALGPS2", "CD79B", "MS4A1", "BANK1", "CD74", "TNFRSF13C", "HLA-DQA1", "IGHM", "MEF2C"), cols =pal)

#####################################
## Do not downsample to 2000 cells per sample and perform unsupervised clustering.
#####################################
Idents(object.sct.filtered) <- 'predicted.celltype.l1'
B_cells <- subset(object.sct.filtered, idents= "B")

DefaultAssay(B_cells) <- "RNA"

B_cells <- NormalizeData(B_cells, assay = "RNA")
B_cells <- FindVariableFeatures(B_cells, assay = "RNA")
B_cells <- ScaleData(B_cells, vars.to.regress = c("S.Score", "G2M.Score"), assay = "RNA", rownames(B_cells))
B_cells <- RunPCA(B_cells, features = rownames(B_cells), assay = "RNA")
B_cells <- FindNeighbors(B_cells, dims = 1:30, assay = "RNA")

B_cells <- FindClusters(
  B_cells, random.seed=12345, verbose=FALSE,
  resolution=seq(0, 1, by=0.05)
)

B_cells <- RunUMAP(B_cells, assay = "RNA", dims = 1:30, reduction = "pca")

head(B_cells)

p <- DimPlot(
  B_cells, reduction="pca", ncol=2,
  group.by=glue::glue("RNA_snn_res.{seq(0, 1, by=0.05)}")
)
p
DimPlot(B_cells, reduction = "umap", group.by = 'RNA_snn_res.0.2', split.by= 'new.ident')
DimPlot(B_cells, reduction = "umap", group.by = 'RNA_snn_res.0.25', split.by= 'new.ident')
DimPlot(B_cells, reduction = "umap", group.by = 'RNA_snn_res.0.3', split.by= 'new.ident')
DimPlot(B_cells, reduction = "umap", group.by = 'RNA_snn_res.0.35', split.by= 'new.ident')
DimPlot(B_cells, reduction = "umap", group.by = 'RNA_snn_res.0.4', split.by= 'new.ident')

table(B_cells$RNA_snn_res.0.3, B_cells$orig.ident)
Idents(B_cells) <- 'RNA_snn_res.0.3'

all.clusters <- FindAllMarkers(B_cells, assay = "RNA")

write.csv(ident.clusters, '11725-JM/ident.clusters.csv')
write.csv(all.clusters, '11725-JM/all.clusters.csv')
#####################################
## Perform DE pairwise analysis.
#####################################

df_contrasts<-t(as.data.frame(combn(levels(as.factor(B_cells@meta.data$new.ident)),2)))
colnames(df_contrasts)<-c("d1","d2")
rownames(df_contrasts)<-c(1:nrow(df_contrasts))
df_contrasts<-as.data.frame(df_contrasts)

B_cells[["RNA"]] <- JoinLayers(B_cells[["RNA"]])
Idents(B_cells) <- 'new.ident'
HC1vsHC2 <- FindMarkers(B_cells, assay = "RNA",ident.1 = 'HC 1', ident.2 = 'HC 2')
HC1vsR331P <- FindMarkers(B_cells, assay = "RNA",ident.1 = 'HC 1', ident.2 = 'R331P/R331P')
HC1vsG126D <- FindMarkers(B_cells, assay = "RNA",ident.1 = 'HC 1', ident.2 = 'WT/G126D')
HC2vsR331P<- FindMarkers(B_cells, assay = "RNA",ident.1 = 'HC 2', ident.2 = 'R331P/R331P')
HC2vsG126D<- FindMarkers(B_cells, assay = "RNA",ident.1 = 'HC 2', ident.2 = 'WT/G126D')
G126DvsR331P<- FindMarkers(B_cells, assay = "RNA",ident.1 = 'WT/G126D', ident.2 = 'R331P/R331P')

write.csv(HC1vsHC2, 'HC1vsHC2.csv')
write.csv(HC1vsR331P, 'HC1vsR331P.csv')
write.csv(HC1vsG126D, 'HC1vsG126D.csv')
write.csv(HC2vsR331P, 'HC2vsR331P.csv')
write.csv(HC2vsG126D, 'HC2vsG126D.csv')
write.csv(G126DvsR331P, 'G126DvsR331P.csv')

saveRDS(B_cells, 'B_cells.rds')
##############
#####select the relevant rows for your analysis here#######
df_contrasts_selec<-df_contrasts[c(1),]
##############
##############
DE_between_condition<-vector(mode = "list", length = nrow(df_contrasts_selec))
DE_between_condition_upregulated<-vector(mode = "list", length = nrow(df_contrasts_selec))
DE_between_condition_downregulated<-vector(mode = "list", length = nrow(df_contrasts_selec))
#############
#############
for (x in c(1:nrow(df_contrasts_selec))) {
  print(x)
  
  cell.count1<-colnames(B_cells)[which(B_cells@meta.data$condition==df_contrasts_selec$d1[x])]
  
  print(length(cell.count1))
  
  cell.count2<-colnames(B_cells)[which(B_cells@meta.data$condition==df_contrasts_selec$d2[x])]
  
  print(length(cell.count2))
  
  min.cell.count<-min(length(cell.count1),length(cell.count2))
  
  if (min.cell.count>=50) {
    
    DE_between_condition[[x]]<-FindMarkers(object = B_cells, ident.1 = df_contrasts_selec$d1[x], ident.2= df_contrasts_selec$d2[x], logfc.threshold = 0.01, min.cells.group=50,test.use = "wilcox")
    
    DE_between_condition[[x]]$Gene<-rownames(DE_between_condition[[x]])
    
    DE_between_condition[[x]]$FC<-2^(DE_between_condition[[x]]$avg_log2FC)
    
    DE_between_condition_upregulated[[x]]<-DE_between_condition[[x]][which(DE_between_condition[[x]]$FC>1),]
    DE_between_condition_downregulated[[x]]<-DE_between_condition[[x]][which(DE_between_condition[[x]]$FC<1),]
  }
  else {
    DE_between_condition_upregulated[[x]]<-data.frame(Result=c("not_enough_cells_to_compare"))
    DE_between_condition_downregulated[[x]]<-data.frame(Result=c("not_enough_cells_to_compare"))
  }
}

#############################
### downsample and subset only B cells for unsupervised clustering
##############################

DefaultAssay(B_cells) <- "RNA"

B_cells <- SplitObject(B_cells, split.by = "orig.ident")
B_cells$control8 <- B_cells$'0008'
B_cells$control9 <- B_cells$'0009'
B_cells$patient15 <- B_cells$'0015'
B_cells$patient16 <- B_cells$'0016'
B_cells$'0008' <- NULL
B_cells$'0009' <- NULL
B_cells$'0015' <- NULL
B_cells$'0016' <- NULL
B_cells$control8 <- subset(B_cells$control8, downsample = 2000)
B_cells$control9 <- subset(B_cells$control9, downsample = 2000)
B_cells$patient15 <- subset(B_cells$patient15, downsample = 2000)
B_cells$patient16 <- subset(B_cells$patient16, downsample = 2000)
B_cells <- merge(B_cells[[1]], B_cells[2:length(B_cells)])

B_cells <- NormalizeData(B_cells)
B_cells <- FindVariableFeatures(B_cells)
B_cells <- ScaleData(B_cells, vars.to.regress = c("S.Score", "G2M.Score"), rownames(B_cells))
B_cells <- RunPCA(B_cells, features = rownames(B_cells), verbose = FALSE)
B_cells<- FindNeighbors(B_cells, dims = 1:30)

for (i in 1:length(seq(0.1,1,0.1))){
  B_cells <- FindClusters(B_cells, cluster.name = "unintegrated_clusters", resolution = seq(0.05,1,0.05)[i])
}

B_cells <- RunUMAP(B_cells, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(B_cells, reduction = "umap.unintegrated", group.by = 'RNA_snn_res.0.35', split.by= 'celltype.genotype')

table(B_cells$RNA_snn_res.0.35, B_cells$orig.ident)
# 0008 0009 0015 0016
#0  124  157    5 1547 # Heterozgyote cluster
#1  825  809    0   35 # Control cluster
#2  687  635    1   11 # Control cluster
#3    3    0  837    6 # Homozygote cluster
#4   15    8  694   16 # Homozygote cluster
#5   98  202    3  236 # Everyone but homozygote
#6  205  153    2  125 # Everyone but homozygote
#7    0    0  410    2 # Homozygote cluster
#8   43   36   48   22 # Everyone


