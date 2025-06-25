#############################################################
# GSEA: B cells cluster 1 ORA and fGSEA from unsupervised clustering with res 0.15
#############################################################
library(msigdbr)
library(enrichplot)
library(clusterProfiler)
library(Seurat)
library(SeuratObject)
library(org.Hs.eg.db)
library(readxl)
library(fgsea)
mem.maxVSize(1e10)

setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/SupplementaryData")

cluster1_GSEA<- read.csv('Cluster 1 homozygote cluster GSEA.csv')

H_t2g <- msigdbr(species = "Homo sapiens", category = "H")
H_t2g <- split(x = H_t2g$gene_symbol, f = H_t2g$gs_name)

head(H_t2g)


ranked_cluster1_GSEA <- -sign(cluster1_GSEA$avg_log2FC)*(log10(cluster1_GSEA$p_val_adj+1e-100))
max(ranked_cluster1_GSEA)
min(ranked_cluster1_GSEA)

names(ranked_cluster1_GSEA) <- cluster1_GSEA$gene
ranked_cluster1_GSEA <- sort(ranked_cluster1_GSEA, decreasing = TRUE)
head(ranked_cluster1_GSEA)

set.seed(20)
cluster1_GSEA<- fgseaMultilevel(stats = ranked_cluster1_GSEA, pathways = H_t2g,
                                scoreType = 'std',
                                minSize = 10,
                                maxSize = 500,
                                nproc = 1)
#cluster1_GSEA

topPathwaysUp <- cluster1_GSEA[ES > 0& padj <0.05][head(order(pval)), pathway]
topPathwaysDown <- cluster1_GSEA[ES < 0 & padj <0.05][head(order(pval)), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(H_t2g[topPathways], stats = ranked_cluster1_GSEA, fgseaRes = cluster1_GSEA, gseaParam = 0.5)
plotEnrichment(pathway= H_t2g[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']], stats = ranked_cluster1_GSEA)+ggtitle('HALLMARK OXIDATIVE PHOSPHORYLATION')

#Extract the gene ranks involved in significantupregulation/downregulation among patients
library(stringi)
res <- as.data.frame((stri_list2matrix(cluster1_GSEA$leadingEdge)))
colnames(res) <- (cluster1_GSEA$pathway)
E2Fgenes <- na.omit(data.frame('genes' = res$HALLMARK_E2F_TARGETS))
E2Fgenes

ranks_cluster1 <- as.data.frame(ranked_cluster1_GSEA)
ranks_cluster1$genes <- rownames(ranks_cluster1)
GSEA_pathways_ranked_E2F_genes <- merge(ranks_cluster1, E2Fgenes, by = 'genes')


Mycgenes <- na.omit(data.frame('genes' = res$HALLMARK_MYC_TARGETS_V1))
Mycgenes


GSEA_pathways_ranked_Myc_genes <- merge(ranks_cluster1, Mycgenes, by = 'genes')



OxPhosgenes <- na.omit(data.frame('genes' = res$HALLMARK_OXIDATIVE_PHOSPHORYLATION))
OxPhosgenes


GSEA_pathways_ranked_OxPhos_genes <- merge(ranks_cluster1, OxPhosgenes, by = 'genes')


setwd("/Volumes/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/GSEA and ORA for Unsupervised Cluster 1/fGSEA")
write.csv(GSEA_pathways_ranked_OxPhos_genes[order(-GSEA_pathways_ranked_OxPhos_genes$ranked_cluster1_GSEA),], 'Ranked gene list of genes involved in Hallmark Oxidative Phosphorylation.csv')
write.csv(GSEA_pathways_ranked_Myc_genes[order(-GSEA_pathways_ranked_Myc_genes$ranked_cluster1_GSEA),], 'Ranked gene list of genes involved in Hallmark Myc Targets V1.csv')
write.csv(GSEA_pathways_ranked_E2F_genes[order(-GSEA_pathways_ranked_E2F_genes$ranked_cluster1_GSEA),], 'Ranked gene list of genes involved in Hallmark E2F Targets.csv')



#############
### C5 gene set (GO terms)
#############
C5 <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = 'BP' )
C5 <- split(x = C5$gene_symbol, f = C5$gs_name)

head(C5)

set.seed(20)
cluster1_GSEA_GO<- fgseaMultilevel(stats = ranked_cluster1_GSEA, pathways = C5,
                                scoreType = 'std',
                                minSize = 10,
                                maxSize = 500,
                                nproc = 1)
cluster1_GSEA_GO

topPathwaysUp <- cluster1_GSEA_GO[ES > 0& padj <0.05][head(order(pval)), pathway]
topPathwaysDown <- cluster1_GSEA_GO[ES < 0 & padj <0.05][head(order(pval)), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(C5[topPathways], stats = ranked_cluster1_GSEA, fgseaRes = cluster1_GSEA_GO, gseaParam = 0.5)



### C6 gene sets (Genes related to cancer)
C6 <- msigdbr(species = "Homo sapiens", category = "C6")
C6 <- split(x = C6$gene_symbol, f = C6$gs_name)

head(C6)

set.seed(20)
cluster1_GSEA_C6<- fgseaMultilevel(stats = ranked_cluster1_GSEA, pathways = C6,
                                   scoreType = 'std',
                                   minSize = 10,
                                   maxSize = 500,
                                   nproc = 1)
cluster1_GSEA_C6

topPathwaysUp <- cluster1_GSEA_C6[ES > 0& padj <0.05][head(order(pval)), pathway]
topPathwaysDown <- cluster1_GSEA_C6[ES < 0 & padj <0.05][head(order(pval)), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(C6[topPathways], stats = ranked_cluster1_GSEA, fgseaRes = cluster1_GSEA_C6, gseaParam = 0.5)# No terms enriched

### C7 gene sets (Genes related to immunology)
C7 <- msigdbr(species = "Homo sapiens", category = "C7", subcollection = 'IMMUNESIGDB')
C7 <- split(x = C7$gene_symbol, f = C7$gs_name)

head(C7)

set.seed(20)
cluster1_GSEA_C7<- fgseaMultilevel(stats = ranked_cluster1_GSEA, pathways = C7,
                                   scoreType = 'std',
                                   minSize = 10,
                                   maxSize = 500,
                                   nproc = 1)
cluster1_GSEA_C7

topPathwaysUp <- cluster1_GSEA_C7[ES > 0& padj <0.05][head(order(pval)), pathway]
topPathwaysDown <- cluster1_GSEA_C7[ES < 0 & padj <0.05][head(order(pval)), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(C7[topPathways], stats = ranked_cluster1_GSEA, fgseaRes = cluster1_GSEA_C7, gseaParam = 0.5)


### KEGG LEGACY C2 gene sets (Genes sets from canonical pathways)
C2 <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = 'CP', subcategory = 'KEGG_LEGACY')
C2 <- split(x = C2$gene_symbol, f = C2$gs_name)

head(C2)

set.seed(20)
cluster1_GSEA_C2<- fgseaMultilevel(stats = ranked_cluster1_GSEA, pathways = C2,
                                   scoreType = 'std',
                                   minSize = 10,
                                   maxSize = 500,
                                   nproc = 1)
cluster1_GSEA_C2

topPathwaysUp <- cluster1_GSEA_C2[ES > 0& padj <0.05][head(order(pval)), pathway]
topPathwaysDown <- cluster1_GSEA_C2[ES < 0 & padj <0.05][head(order(pval)), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(C2[topPathways], stats = ranked_cluster1_GSEA, fgseaRes = cluster1_GSEA_C2, gseaParam = 0.5)


### REACTOME C2 gene sets (Genes sets from canonical pathways)
reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = 'CP', subcategory = 'REACTOME')
reactome <- split(x = reactome$gene_symbol, f = reactome$gs_name)

head(reactome)

set.seed(20)
cluster1_GSEA_reactome<- fgseaMultilevel(stats = ranked_cluster1_GSEA, pathways = reactome,
                                   scoreType = 'std',
                                   minSize = 10,
                                   maxSize = 500,
                                   nproc = 1)
cluster1_GSEA_reactome

topPathwaysUp <- cluster1_GSEA_reactome[ES > 0& padj <0.05][head(order(pval)), pathway]
topPathwaysDown <- cluster1_GSEA_reactome[ES < 0 & padj <0.05][head(order(pval)), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(reactome[topPathways], stats = ranked_cluster1_GSEA, fgseaRes = cluster1_GSEA_reactome, gseaParam = 0.5)


#######################
### Now, use clusterprofiler to see if the DEGs are overrepresented within genes sets (KEGG, HALLMARK, C7, and C5).
#######################

universe <- read.csv('FindAllMarkers output wilcoxon rank-sum test.csv')
universe <- unique(universe$gene)

H_t2g <- msigdbr(species = "Homo sapiens", category = "H")%>%
  dplyr::select(gs_name, entrez_gene)


All_clusters <- read.csv('FindAllMarkers output wilcoxon rank-sum test.csv')
universe <- select(org.Hs.eg.db, 
                 keys = All_clusters$gene,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL")
universe <- universe[!duplicated(universe), ]

#### Extract only sig DEGs from FindAllMarkers output
DEGs <- All_clusters[All_clusters$avg_log2FC > 0.25 & All_clusters$p_val_adj<0.05,]
length(unique(DEGs$gene))
#7752

#### Map gene symbols to entrezids
entrez <- select(org.Hs.eg.db, 
            keys = DEGs$gene,
            columns = c("ENTREZID", "SYMBOL"),
            keytype = "SYMBOL")

entrez$gene <- entrez$SYMBOL
entrez <-  entrez[!duplicated(entrez), ]
DEGs <- merge(entrez, DEGs, by = 'gene')
DEGs<- na.omit(DEGs)

DEGs_cluster <- split(x = DEGs$ENTREZID, f = DEGs$cluster)

set.seed(20)
ORA_enrichGO <- clusterProfiler::compareCluster(gene = DEGs_cluster, 
                                                OrgDb = "org.Hs.eg.db", 
                                                ont = 'BP', 
                                                fun = 'enrichGO',
                                                universe = na.omit(universe$ENTREZID))

### Find GO terms that are up in clusters 1 and 0 and not up among all other clusters
#Duplicates <- ORA_enrichGO@compareClusterResult[duplicated(ORA_enrichGO@compareClusterResult$Description), ] #We want to identify pathways that are specific to each cluster.
#ORA_enrichGO@compareClusterResult<- ORA_enrichGO@compareClusterResult[ORA_enrichGO@compareClusterResult$Description %!in% Duplicates$Description, ] 
#ORA_enrichGO@compareClusterResult <- ORA_enrichGO@compareClusterResult%>%
#  group_by(Cluster)%>%
#  arrange(GeneRatio) %>%
#  dplyr::slice(1:5)

dotplot(ORA_enrichGO, showCategory = 2)

#### Enriched Hallmark gene sets
set.seed(20)
ORA_Hallmark <- clusterProfiler::compareCluster(gene = DEGs_cluster,  
                                                TERM2GENE = H_t2g,
                                                fun = 'enricher',
                                                universe = na.omit(universe$ENTREZID))
dotplot(ORA_Hallmark, showCategory = Inf)

### C7 immunesigdb
C7 <- msigdbr(species = "Homo sapiens", category = "C7",subcollection = 'IMMUNESIGDB')%>%
  dplyr::select(gs_name, entrez_gene)

set.seed(20)
ORA_C7 <- clusterProfiler::compareCluster(gene = DEGs_cluster,  
                                                TERM2GENE = C7,
                                                fun = 'enricher',
                                                universe = na.omit(universe$ENTREZID))

Duplicates <- ORA_C7@compareClusterResult[duplicated(ORA_C7@compareClusterResult$Description), ] #We want to identify pathways that are specific to each cluster.
ORA_C7@compareClusterResult<- ORA_C7@compareClusterResult[ORA_C7@compareClusterResult$Description %!in% Duplicates$Description, ] 
ORA_C7@compareClusterResult <- ORA_C7@compareClusterResult%>%
  group_by(Cluster)%>%
  arrange(qvalue) %>%
  dplyr::slice(1:5)

dotplot(ORA_C7)


### KEGG
set.seed(20)
ORA_kegg <- clusterProfiler::compareCluster(gene = DEGs_cluster,  
                                          keyType = "kegg",
                                          organism = 'hsa',
                                          fun = 'enrichKEGG',
                                          universe = na.omit(universe$ENTREZID))

Duplicates <- ORA_kegg@compareClusterResult[duplicated(ORA_kegg@compareClusterResult$Description), ] #We want to identify pathways that are specific to each cluster.
ORA_kegg@compareClusterResult<- ORA_kegg@compareClusterResult[ORA_kegg@compareClusterResult$Description %!in% Duplicates$Description, ] 
ORA_kegg@compareClusterResult <- ORA_kegg@compareClusterResult%>%
  group_by(Cluster)%>%
  arrange(qvalue) %>%
  dplyr::slice(1:5)

dotplot(ORA_kegg, showCategory = 5)


### enrichDO
set.seed(20)
ORA_DO <- clusterProfiler::compareCluster(gene = DEGs_cluster,
                                            fun = 'enrichDO',
                                            universe = na.omit(universe$ENTREZID))

Duplicates <- ORA_DO@compareClusterResult[duplicated(ORA_DO@compareClusterResult$Description), ] #We want to identify pathways that are specific to each cluster.
ORA_DO@compareClusterResult<- ORA_DO@compareClusterResult[ORA_DO@compareClusterResult$Description %!in% Duplicates$Description, ]

ORA_DO@compareClusterResult <- ORA_DO@compareClusterResult%>%
  group_by(Cluster)%>%
  arrange(qvalue) %>%
  dplyr::slice(1:5)

dotplot(ORA_DO, showCategory = 5)
