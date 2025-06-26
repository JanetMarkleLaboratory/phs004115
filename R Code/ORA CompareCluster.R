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


#### Map gene symbols to entrezids
entrez <- select(org.Hs.eg.db, 
                 keys = All_clusters$gene,
                 columns = c("ENTREZID", "SYMBOL"),
                 keytype = "SYMBOL")

entrez$gene <- entrez$SYMBOL
entrez <-  entrez[!duplicated(entrez), ]
All_clusters <- merge(entrez, All_clusters, by = 'gene')
All_clusters <- na.omit(All_clusters)

All_clusters<- All_clusters[All_clusters$p_val_adj < 0.05,]
All_clusters <- All_clusters[(All_clusters$avg_log2FC > 0.25)| (All_clusters$avg_log2FC < -0.25),]

All_clusters$group <- "upregulated"
All_clusters$group[All_clusters$avg_log2FC < 0] <- "downregulated"


set.seed(20)
formula_res <- compareCluster(ENTREZID~group+cluster, data=All_clusters, OrgDb = "org.Hs.eg.db", 
                              ont = 'BP', 
                              fun = 'enrichGO',
                              universe = na.omit(universe$ENTREZID))

head(formula_res)

dotplot(formula_res, x="group", showCategory = 2) + 
  facet_grid(~cluster)+
  theme_light()+
  theme(axis.text.x = element_text(angle= 90))+
  xlab('')+
  ggtitle('GO BP terms across Cluster')


# Hallmark gene sets
set.seed(20)
formula_hallmark <- compareCluster(ENTREZID~group+cluster, data=All_clusters, TERM2GENE = H_t2g,
                                   fun = 'enricher',
                                   universe = na.omit(universe$ENTREZID))

head(formula_hallmark)

dotplot(formula_hallmark, x="group") + 
  facet_grid(~cluster)+
  theme_light()+
  theme(axis.text.x = element_text(angle= 90))+
  xlab('')+
  ggtitle('Hallmark Gene Sets across Cluster')



# Disease ontology gene sets
set.seed(20)
ORA_DO <- clusterProfiler::compareCluster(ENTREZID~group+cluster, 
                                          data=All_clusters,
                                          fun = 'enrichDO',
                                          universe = na.omit(universe$ENTREZID))
dotplot(ORA_DO, x="group") + 
  facet_grid(~cluster)+
  theme_light()+
  theme(axis.text.x = element_text(angle= 90))+
  xlab('')+
  ggtitle('Top 10 Disease Ontology Gene Sets across Cluster')



# KEGG 
set.seed(20)
formula_kegg <- clusterProfiler::compareCluster(ENTREZID~group+cluster, 
                                                data=All_clusters,
                                                keyType = "kegg",
                                                organism = 'hsa',
                                                fun = 'enrichKEGG',
                                                universe = na.omit(universe$ENTREZID))
dotplot(formula_kegg, x="group") + 
  facet_grid(~cluster)+
  theme_light()+
  theme(axis.text.x = element_text(angle= 90))+
  xlab('')+
  ggtitle('Top 10 KEGG Gene Sets across Cluster')



#C7 gene set
set.seed(20)
formula_C7 <- clusterProfiler::compareCluster(ENTREZID~group+cluster, 
                                                data=All_clusters,
                                                TERM2GENE = C7,
                                                fun = 'enricher',
                                                universe = na.omit(universe$ENTREZID))
dotplot(formula_C7, x="group") + 
  facet_grid(~cluster)+
  theme_light()+
  theme(axis.text.x = element_text(angle= 90))+
  xlab('')+
  ggtitle('Top 10 C7 ImmuneSigDB Gene Sets across Cluster')

