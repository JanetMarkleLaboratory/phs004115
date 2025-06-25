#############################################################
# PART 7: GSEA: Intermediate and Naïve B cells and CD14 Monocytes
## NOTE: Only Intermediate B cells were deemed interesting to continue further investigation.
#############################################################
library(msigdbr)
library(enrichplot)
library(clusterProfiler)
library(Seurat)
library(SeuratObject)
library(org.Hs.eg.db)
library(readxl)
library(fgsea)
options(enrichment_force_universe=F)

Naive_Markers <- read_excel("GSEA Naive B.xlsx")
H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(H_t2g)

Naive_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = Naive_Markers$gene,
                                      columns = c('ENTREZID'), keytype = 'SYMBOL')%>%
  filter(!duplicated(SYMBOL))
Naive_Markers$ENTREZ <- Naive_entrez$ENTREZID
Naive_Markers <- na.omit(Naive_Markers)

length(Naive_Markers$gene)
#10894
ranked_naive <- Naive_Markers$avg_log2FC*(-log10(Naive_Markers$p_val_adj+10e-9))
names(ranked_naive) <- Naive_Markers$ENTREZ
ranked_naive <- sort(ranked_naive, decreasing = TRUE)
head(ranked_naive)
naive_GSEA<- clusterProfiler::GSEA(ranked_naive, maxGSSize = 500, TERM2GENE = H_t2g)

############################
### Intermediate B cells
############################
Int_B <- read_excel("GSEA Intermediate B.xlsx")
H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(H_t2g)

Int_entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = Int_B$gene,
                                    columns = c('ENTREZID'), keytype = 'SYMBOL')

Int_B$ENTREZ <- Int_entrez$ENTREZID
Int_B <- na.omit(Int_B)
length(Int_B$gene)
#10821
ranked_intermediate <- Int_B$avg_log2FC*(-log10(Int_B$p_val_adj))
names(ranked_intermediate) <- Int_B$ENTREZ
ranked_intermediate <- sort(ranked_intermediate, decreasing = TRUE)
head(ranked_intermediate)
intermediate_GSEA<- clusterProfiler::GSEA(ranked_intermediate, minGSSize = 10, maxGSSize = 500, TERM2GENE = H_t2g)
intermediate_GSEA

m_df<- msigdbr(species = "Homo sapiens", category = "H")

head(m_df)

fgsea_sets<- m_df %>% split(x = .$entrez_gene, f = .$gs_name)
fgsea_int <- fgsea(pathways = fgsea_sets, stats = ranked_intermediate, minSize = 10, maxSize = 500)
gseaplot2(geneSetID = 1,pvalue_table = T,  intermediate_GSEA, color =  'green',title = 'Patient Intermediate B cells')
plotEnrichment(fgsea_sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],ranked_intermediate) + labs(title='Patient Intermediate B cells')
############################
### CD14 Monos
############################
CD14_monos <- read_excel("GSEA CD14 Monos.xlsx")

H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)
head(H_t2g)

CD14_Entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = CD14_monos$gene,
                                     columns = c('ENTREZID'), keytype = 'SYMBOL')%>%
  filter(!duplicated(SYMBOL))

CD14_monos$ENTREZ <- CD14_Entrez$ENTREZID
CD14_monos <- na.omit(CD14_monos)
length(CD14_monos$gene)
#9946
ranked_CD14_monos <- CD14_monos$avg_log2FC*(-log10(CD14_monos$p_val_adj+10e-9))
names(ranked_CD14_monos) <- CD14_monos$ENTREZ
ranked_CD14_monos <- sort(ranked_CD14_monos, decreasing = TRUE)
head(ranked_CD14_monos)
CD14_monos_GSEA<- clusterProfiler::GSEA(ranked_CD14_monos, maxGSSize = 500, 
                                        TERM2GENE = H_t2g)
gseaplot2(geneSetID = 1:2, CD14_monos_GSEA,title = 'Patient CD14 Monocytes')
fgsea_CD14 <- fgsea(pathways = fgsea_sets, stats = ranked_CD14_monos, nperm = 10000)
