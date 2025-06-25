library(readxl)
library(Seurat)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(forcats)
library(DOSE)
######################
##Read in the necessary files
######################
setwd("/Volumes/pmi/Data/Research/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/11725-JM")
heatmapList <- read_excel("heatmapList.xlsx")
Universe_FindMarkers <- read_excel("Universe_FindMarkers.xlsx")
geneUniverse <- Universe_FindMarkers$GeneUniverse
length(geneUniverse)
geneUniverse <- AnnotationDbi::select(org.Hs.eg.db, keys = geneUniverse,
                                      columns = c('ENTREZID'), keytype = 'SYMBOL')
universe<-na.omit(geneUniverse)
#universe <- universe[!duplicated(universe$SYMBOL), ]
length(universe$SYMBOL)
View(heatmapList)
nrow(heatmapList)
#610 genes from list and 11201 genes in the universe

B_cells <- readRDS('B_cells.rds')
Idents(B_cells) <- 'new.ident'
######################
# Run two FindMarkers analyses on the gene list
######################
Heatmap_geneList_prefiltered <- FindMarkers(B_cells, assay = 'RNA',features = heatmapList$`Genes Passed Filtering`, 
                                            ident.2 = c("HC 1", "HC 2"), ident.1 = c("WT/G126D","R331P/R331P"))

Heatmap_geneList_R331Pvsall <- FindMarkers(B_cells, assay = 'RNA',
                                      features = heatmapList$`Genes Passed Filtering`, 
                                      ident.2 = c("HC 1", "HC 2", "WT/G126D"), ident.1 = "R331P/R331P")
write.csv(Heatmap_geneList_prefiltered, 'FindMarkers_analysis.1.csv')
write.csv(Heatmap_geneList_R331Pvsall, 'FindMarkers_analysis.2.csv')

Heatmap_geneList_prefiltered<-read.csv('FindMarkers_analysis.1.csv')
Heatmap_geneList_R331Pvsall<- read.csv('FindMarkers_analysis.2.csv')
######################
# positive log2FC indicates Up-regulation among patients, while negative values indicate up-regulation among controls and/or heterozygote.
######################
Heatmap_geneList_patients <- Heatmap_geneList_prefiltered[Heatmap_geneList_prefiltered$avg_log2FC>0,]
Heatmap_geneList_patients$gene <- Heatmap_geneList_patients$X

Heatmap_geneList_HCs <- Heatmap_geneList_prefiltered[Heatmap_geneList_prefiltered$avg_log2FC<0,]
Heatmap_geneList_HCs$gene <- Heatmap_geneList_HCs$X

Heatmap_geneList_R331P <- Heatmap_geneList_R331Pvsall[Heatmap_geneList_R331Pvsall$avg_log2FC>0,]
Heatmap_geneList_R331P$gene <- Heatmap_geneList_R331P$X

Heatmap_geneList_HCs_HET <- Heatmap_geneList_R331Pvsall[Heatmap_geneList_R331Pvsall$avg_log2FC<0,]
Heatmap_geneList_HCs_HET$gene <- Heatmap_geneList_HCs_HET$X
######################
## Perform ORA and compare all gene lists. Remove VAX subcategory from C7 gene set.
######################
C7_gene_sets = msigdbr(species = "Homo sapiens", category = "C7", subcategory = 'IMMUNESIGDB') 
length(unique(C7_gene_sets$gs_name))
C7_gene_sets<- C7_gene_sets %>% 
  dplyr::select(gs_name, entrez_gene)

entrez_patients <- AnnotationDbi::select(org.Hs.eg.db, keys = Heatmap_geneList_patients$gene,
                                       columns = c('ENTREZID'), keytype = 'SYMBOL')

entrez_patients <- na.omit(entrez_patients)
entrez_HCs <- AnnotationDbi::select(org.Hs.eg.db, keys = Heatmap_geneList_HCs$gene,
                                         columns = c('ENTREZID'), keytype = 'SYMBOL')

entrez_HCs <- na.omit(entrez_HCs)

entrez_R331P <- AnnotationDbi::select(org.Hs.eg.db, keys = Heatmap_geneList_R331P$gene,
                                      columns = c('ENTREZID'), keytype = 'SYMBOL')

entrez_R331P <- na.omit(entrez_R331P)
entrez_HCs_HET <- AnnotationDbi::select(org.Hs.eg.db, keys = Heatmap_geneList_HCs_HET$gene,
                                        columns = c('ENTREZID'), keytype = 'SYMBOL')

entrez_HCs_HET <- na.omit(entrez_HCs_HET)
####
comparelist_ALL<-list(entrez_HCs$ENTREZID,entrez_patients$ENTREZID, entrez_HCs_HET$ENTREZID,entrez_R331P$ENTREZID)
names(comparelist_ALL)<-c("Up-regulated in HCs","Up-regulated in CARD11 patients", "Up-regulated in HCs and WT/G126D","Up-regulated in R331P/R331P") 
H_gene_sets = msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
H_clust<-compareCluster(geneCluster = comparelist_ALL, 
                                 fun = "enricher",pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                 universe = universe$ENTREZID, 
                                 minGSSize = 15, maxGSSize = 500, TERM2GENE = H_gene_sets
)


C7_clust<-compareCluster(geneCluster = comparelist_ALL, 
                                  fun = "enricher", pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                                  universe = universe$ENTREZID, 
                                  minGSSize = 15, maxGSSize = 500, TERM2GENE = C7_gene_sets
)

x <- as.data.frame(C7_clust)

top_5_overall <- x %>%
  group_by(Cluster) %>%
  top_n(5, wt = FoldEnrichment)

top_5_overall$Description <-gsub('_', ' ', top_5_overall$Description)
top_5_overall$Description = str_wrap(top_5_overall$Description, width = 10)

top_5_HCs <- top_5_overall[top_5_overall$Cluster == "Up-regulated in HCs",]
top_5_R331P <- top_5_overall[top_5_overall$Cluster == "Up-regulated in R331P/R331P",]
top_5_HCs_Het <- top_5_overall[top_5_overall$Cluster == "Up-regulated in HCs and WT/G126D",]
top_5_Patients <- top_5_overall[top_5_overall$Cluster == "Up-regulated in CARD11 patients",]

scaleFUN <- function(x) sprintf("%.2f", x)

HCs <- ggplot(top_5_HCs,
  aes(RichFactor, fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  facet_grid(~Cluster)+
  geom_point(aes(color=qvalue, size = Count)) +
  scale_color_gradientn(colours=c("darkblue", "lightblue"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +theme_dose(12) +
  ylab("")+ theme(text = element_text(family = 'arial', size = 12),
                  strip.text = element_text(face='bold',family = 'arial', size = 14))+
  scale_size_continuous(range=c(2, 10)) +
  xlab("Rich Factor") + scale_x_continuous(labels = scaleFUN)

HCs_Het <- ggplot(top_5_HCs_Het,
              aes(RichFactor, fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  facet_grid(~Cluster)+
  geom_point(aes(color=qvalue, size = Count)) +
  scale_color_gradientn(colours=c("darkblue", "lightblue"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +theme_dose(12) +
  ylab("")+ theme(text = element_text(family = 'arial', size = 12),
                  strip.text = element_text(face='bold',family = 'arial', size = 14))+
  scale_size_continuous(range=c(2, 10)) + scale_x_continuous(labels = scaleFUN)+
  xlab("Rich Factor") 

R331P <- ggplot(top_5_R331P,
              aes(RichFactor, fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  facet_grid(~Cluster)+
  geom_point(aes(color=qvalue, size = Count)) + scale_x_continuous(labels = scaleFUN)+
  scale_color_gradientn(colours=c("darkblue", "lightblue"),
                        
                        trans = "log10",
                        
                        guide=guide_colorbar(reverse=TRUE, order=1)) +theme_dose(12) +
  ylab("")+ theme(text = element_text(family = 'arial', size = 12),
                  strip.text = element_text(face='bold',family = 'arial', size = 14))+
  scale_size_continuous(range=c(2, 10)) +
  xlab("Rich Factor")

patients <- ggplot(top_5_Patients,
              aes(RichFactor, fct_reorder(Description, RichFactor))) +
  geom_segment(aes(xend=0, yend = Description)) +
  facet_grid(~Cluster)+
  geom_point(aes(color=qvalue, size = Count)) +
  scale_color_gradientn(colours=c("darkblue", "lightblue"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +theme_dose(12) +
  ylab("")+ theme(text = element_text(family = 'arial', size = 12),
                  strip.text = element_text(face='bold',family = 'arial', size = 14))+
  scale_size_continuous(range=c(2, 10)) +
  xlab("Rich Factor") + scale_x_continuous(labels = scaleFUN)

cowplot::plot_grid(HCs, HCs_Het, R331P, patients, label_fontfamily = 'arial', align = "hv")
title <- cowplot::ggdraw() + 
  cowplot::draw_label(
    "C7: Immunologic Signature Gene Set",
    fontface = 'bold', fontfamily = 'arial',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
cowplot::plot_grid(
  title, NULL, HCs, patients, HCs_Het, R331P, label_fontfamily = 'arial', align = "hv",
  ncol = 2,
  # rel_heights values control vertical title margins
  rel_heights = c(0.1, 1,1,1)
)
##### Find genes that overlap the C7 gene sets
genes <- AnnotationDbi::select(org.Hs.eg.db, keys = heatmapList$`Genes Passed Filtering`,
                               columns = c('ENTREZID'), keytype = 'SYMBOL')
genes <- na.omit(genes)
genes_C7<-enricher(genes$ENTREZID, pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                   universe = universe$ENTREZID, 
                   minGSSize = 15, maxGSSize = 500, TERM2GENE = C7_gene_sets
)

genes_C7 <- as.data.frame(genes_C7)
subset_genes_C7 <- genes_C7[genes_C7$ID %in% c('GSE22886_NAIVE_VS_IGG_IGA_MEMORY_BCELL_DN',"GSE22886_NAIVE_VS_IGM_MEMORY_BCELL_DN"),]
df <- strsplit(subset_genes_C7$geneID, split = '/')
df<- unlist(df)
length(df)
df <- unique(df)
length(df)
df
genes <- AnnotationDbi::select(org.Hs.eg.db, keys = df,
                                      columns = c('SYMBOL'), keytype = 'ENTREZID')
View(genes)
write.csv(genes, '53genes.csv')
## 53 genes identified in 20 pathways.

## GSEA rank by log2FC and p-value from selected features

ranking_R331P <- Heatmap_geneList_R331P$avg_log2FC*(-log10(Heatmap_geneList_R331P$p_val_adj+10e-9))

names(ranking_R331P) <- Heatmap_geneList_R331P$gene
ranking_R331P <- sort(ranking_R331P, decreasing = TRUE)

write.csv(ranking_R331P, 'ranking_R331P.csv')


##Run FindMarkers for GSEA with no hyperparameters

GSEA_FindMarkers_Patients <- FindMarkers(B_cells, assay = 'RNA',min.pct = 0.001, logfc.threshold = 0,
                                         ident.2 = c("HC 1", "HC 2"), ident.1 = c("R331P/R331P", "WT/G126D"))
GSEA_FindMarkers_R331P <- FindMarkers(B_cells, assay = 'RNA',min.pct = 0.001, logfc.threshold = 0,
                                         ident.2 = c("HC 1", "HC 2", "WT/G126D"), ident.1 = c("R331P/R331P"))
write.csv(GSEA_FindMarkers_Patients, 'GSEA_FindMarkers_Patients.csv')
write.csv(GSEA_FindMarkers_R331P, 'GSEA_FindMarkers_R331P.csv')
## fgsea: rank by log2FC and p-value




# Create violin plots for 53 genes

B_cells$newer.ident <- NA
B_cells$newer.ident[which(B_cells$orig.ident == "0008")] <- "HC"
B_cells$newer.ident[which(B_cells$orig.ident == "0009")] <- "HC"
B_cells$newer.ident[which(B_cells$orig.ident == "0015")] <- "Homozygote"
B_cells$newer.ident[which(B_cells$orig.ident == "0016")] <- "Heterozygote"
B_cells$newer.ident <- factor(x = B_cells$newer.ident, levels = c("HC", "Heterozygote", "Homozygote"))
B_cells$new.ident <- factor(x = B_cells$new.ident, levels = c("HC 1", "HC 2", "WT/G126D", "R331P/R331P"))
colors <- c(
  "HC" = "blue",
  "HC" = "blue",
  "Homozygote" = "red",
  "Heterozygote" = "orange"
)

log2_expression_data <- GetAssayData(B_cells, assay = "RNA", slot = "data")

log2_expression_data <- as.data.frame(as.table(as.matrix(log2_expression_data)))
colnames(log2_expression_data) <- c("gene", "cell", "Expression")

genes_fit <- c('CD58', 'TOR3A','CFLAR',
               'SP140','ZBTB20','ZBTB38','MFSD10',
               'PXDC1', 'VOPP1','MDFIC', 'KLF10','TUBB4B',
               'ITGB1','IL10RA','GAPDH',
               'TPI1', 'TESC', 'RORA', 'RPLP1',
               'CNP', 'PHB1','ACTG1','EEF2','CD70','PPP1R16B','SAMSN1','SMARCB1')
subset_log2_expression_data <- log2_expression_data[log2_expression_data$gene %in% genes_fit,]
meta.data <- B_cells@meta.data
meta.data$cell <- rownames(meta.data)
length(unique(log2_expression_data$cell))
length(unique(rownames(meta.data)))

merge_log2_expression_data<- merge(subset_log2_expression_data, meta.data, by = "cell")
for (i in genes_fit){
  name <- paste("merge_log2_expression_data", i, sep = '')
  assign(x =name, value = merge_log2_expression_data[merge_log2_expression_data$gene == i,])
}

plot1<- ggplot(merge_log2_expression_dataCD58, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  facet_wrap(~gene, ncol = 5) + coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.x  = element_blank(),
        axis.text.y = element_text(size = 10, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
                #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot2<- ggplot(merge_log2_expression_dataTOR3A, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x  = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot3<- ggplot(merge_log2_expression_dataCFLAR, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 3))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot4<- ggplot(merge_log2_expression_dataSP140, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 2))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot5<- ggplot(merge_log2_expression_dataZBTB20, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 2))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot6<- ggplot(merge_log2_expression_dataZBTB38, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 0.25))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot7<- ggplot(merge_log2_expression_dataMFSD10, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot8<- ggplot(merge_log2_expression_dataPXDC1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 0.25))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot9<- ggplot(merge_log2_expression_dataVOPP1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 2))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot10<- ggplot(merge_log2_expression_dataMDFIC, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 0.5))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot11<- ggplot(merge_log2_expression_dataKLF10, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot12<- ggplot(merge_log2_expression_dataTUBB4B, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 2))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot13<- ggplot(merge_log2_expression_dataITGB1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot14<- ggplot(merge_log2_expression_dataIL10RA, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot15<- ggplot(merge_log2_expression_dataGAPDH, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 2))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot16<- ggplot(merge_log2_expression_dataTPI1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot17<- ggplot(merge_log2_expression_dataTESC, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
  coord_cartesian(ylim = c(0, 0.5))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot18<- ggplot(merge_log2_expression_dataRORA, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 0.5))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot19<- ggplot(merge_log2_expression_dataRPLP1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+ coord_cartesian(ylim = c(0, 6))+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+ facet_grid(~gene)+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot20<- ggplot(merge_log2_expression_dataCNP, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+ facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('') + coord_cartesian(ylim = c(0, 0.25))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot21<- ggplot(merge_log2_expression_dataPHB1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 0.25))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_blank(),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot22<- ggplot(merge_log2_expression_dataACTG1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_text(size = 10, color = "black", family = 'arial',
                                   angle = 45, hjust = 1),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot23<- ggplot(merge_log2_expression_dataEEF2, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 4))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_text(size = 10, color = "black", family = 'arial',
                                   angle = 45, hjust = 1),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot24<- ggplot(merge_log2_expression_dataCD70, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 0.5))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_text(size = 10, color = "black", family = 'arial',
                                   angle = 45, hjust = 1),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot25<- ggplot(merge_log2_expression_dataPPP1R16B, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 2))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_text(size = 10, color = "black", family = 'arial',
                                   angle = 45, hjust = 1),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot26<- ggplot(merge_log2_expression_dataSAMSN1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_text(size = 10, color = "black", family = 'arial',
                                   angle = 45, hjust = 1),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", show.legend = T,
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)
plot27<- ggplot(merge_log2_expression_dataSMARCB1, aes(x= factor(new.ident), y = Expression, fill = factor(newer.ident)))+
  geom_violin(position = "dodge", draw_quantiles = T)+facet_grid(~gene)+
  scale_fill_manual(values = colors)+ labs(x = ' ',y = '')+ xlab('')+
   coord_cartesian(ylim = c(0, 1))+
  theme(strip.background =element_rect(fill="#FFFFD8"),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text.y.right = element_text(size = 10, color = "black", family = 'arial'),
        axis.text.x = element_text(size = 10, color = "black", family = 'arial',
                                   angle = 45, hjust = 1),
        axis.text.y.left = element_text(size = 10, color = "black", family = 'arial'),
        axis.title = element_text(size = 12, color = "black", family = 'arial'),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.5,
               #fatten = ,
               color = "black")+
  NoLegend()+ scale_y_continuous(labels=scaleFUN)

figure<- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8,
             plot9,plot10,plot11,plot12,plot13,plot14,plot15,plot16,plot17,plot18,plot19,plot20,plot21,plot22,
             plot23,plot24,plot25,plot26,plot27, align = "h")
annotate_figure(figure,
                #top = text_grob("Visualizing mpg", color = "red", face = "bold", size = 14),
        
                left = text_grob("Expression (log2)", color = "black", rot = 90, family = 'arial', size = 14)
                #fig.lab = "Figure 1", fig.lab.face = "bold"
)

mean <- 
  merge_log2_expression_data %>% 
  group_by(orig.ident, gene) %>% 
  summarise(
    mean_exp = mean(Expression, na.rm = TRUE),
    sd_exp = sd(Expression, na.rm = TRUE))

write.csv(mean, '58genes.csv')



