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

setwd("/Volumes/pmi/Data/Research/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/11725-JM")

B_cells <- readRDS("B_cells.rds")
B_cells$newer.ident <- NA
B_cells$newer.ident[which(B_cells$orig.ident == "0008")] <- "HC"
B_cells$newer.ident[which(B_cells$orig.ident == "0009")] <- "HC"
B_cells$newer.ident[which(B_cells$orig.ident == "0015")] <- "Homozygote"
B_cells$newer.ident[which(B_cells$orig.ident == "0016")] <- "Heterozygote"
B_cells$newer.ident <- factor(x = B_cells$newer.ident, levels = c("HC", "Heterozygote", "Homozygote"))
B_cells$new.ident <- factor(x = B_cells$new.ident, levels = c("HC 1", "HC 2", "WT/G126D", "R331P/R331P"))


df <- Embeddings(B_cells, reduction = "umap")
df <- as.data.frame(df)
colnames(df) <- c("UMAP1", "UMAP2")
metadata <- B_cells@meta.data
# Include clusters in the metadata
df <- cbind(df, metadata)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


cluster_colors <- c(
  "HC" = "blue",
  "Heterozygote" = "orange",
  "Homozygote" = "red"
)


end_colors <- lighten(cluster_colors, amount = 0.8)


blend_colors <- function(cluster, density, colors, colors_end) {
  base_color <- colors[cluster]
  end_color <- colors_end[cluster]
  blend_color <- colorRampPalette(c(base_color, end_color))(100)[as.integer(density * 100) + 1]
  return(blend_color)
}

dfc1 <- df[df$new.ident  == 'HC 1',]
dfc2 <- df[df$new.ident  == 'HC 2',]
dfp15 <- df[df$new.ident  == 'R331P/R331P',]
dfp16 <- df[df$new.ident  == 'WT/G126D',]

### Control 1 normalized
dfc1$Density <- get_density(dfc1$UMAP1, dfc1$UMAP2, n =500)

dfc1 <- dfc1 %>%
  mutate(normalized_density = (Density - min(Density)) / (max(Density) - min(Density)))


# Apply blending function
dfc1$color <- mapply(blend_colors, dfc1$newer.ident, dfc1$normalized_density, MoreArgs = list(colors = end_colors, colors_end = cluster_colors))

### Control 2 normalized
dfc2$Density <- get_density(dfc2$UMAP1, dfc2$UMAP2, n =500)

dfc2 <- dfc2 %>%
  mutate(normalized_density = (Density - min(Density)) / (max(Density) - min(Density)))


# Apply blending function
dfc2$color <- mapply(blend_colors, dfc2$newer.ident, dfc2$normalized_density, MoreArgs = list(colors = end_colors, colors_end = cluster_colors))

### Patient 15 normalized
dfp15$Density <- get_density(dfp15$UMAP1, dfp15$UMAP2, n =500)

dfp15 <- dfp15 %>%
  mutate(normalized_density = (Density - min(Density)) / (max(Density) - min(Density)))


# Apply blending function
dfp15$color <- mapply(blend_colors, dfp15$newer.ident, dfp15$normalized_density, MoreArgs = list(colors = end_colors, colors_end = cluster_colors))


### Patient 16 normalized
dfp16$Density <- get_density(dfp16$UMAP1, dfp16$UMAP2, n =500)

dfp16 <- dfp16 %>%
  mutate(normalized_density = (Density - min(Density)) / (max(Density) - min(Density)))

# Apply blending function
dfp16$color <- mapply(blend_colors, dfp16$newer.ident, dfp16$normalized_density, MoreArgs = list(colors = end_colors, colors_end = cluster_colors))

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(3, "cm")
)

umap_plot_comb <- ggplot(df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(colour = RNA_snn_res.0.15), size = 0.5) +coord_equal(ylim =c(-6.014122,7.909906), xlim=c(-6.096731,9.633184))+
  scale_colour_manual(breaks = c("0","1","2", "3","4","5","6"), 
                      values = safe_colorblind_palette[1:10])+
  theme_minimal(base_size = 14) +
  annotate(geom = "label", label = "0",
           x = -2, y = -2, size = 8, colour = "black", family = 'arial')+
  annotate(geom = "label", label = "1",
           x = 6, y = 1, size = 8, colour = "black", family = 'arial')+
  annotate(geom = "label", label = "2",
           x = -3.5, y = 2.5, size = 8, colour = "black", family = 'arial')+
  annotate(geom = "label", label = "3",
           x = -2, y = 5.5, size = 8, colour = "black", family = 'arial')+
  annotate(geom = "label", label = "4",
           x = 0.5, y = 0.5, size = 8, colour = "black", family = 'arial')+
  annotate(geom = "label", label = "5",
           x = 3, y = -3, size = 8, colour = "black", family = 'arial')+
  annotate(geom = "label", label = "6",
           x = 2, y = 6.5, size = 8, colour = "black", family = 'arial')+
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")+ggtitle("Combined")+theme(panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank())+theme(aspect.ratio = 1)
b <- ggplot_build(umap_plot_comb)
b$layout$panel_params[[1]]$x.range
b$layout$panel_params[[1]]$y.range

umap_plot_c1 <- ggplot(dfc1, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = color), alpha = 0.4, size = 0.5)+coord_equal(ylim =c(-6.014122,7.909906), xlim=c(-6.096731,9.633184))+
  scale_color_identity()+
  theme_minimal(base_size = 14) +
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none") + ggtitle('HC 1')+theme(panel.grid.major = element_blank(),
                                                          panel.grid.minor = element_blank())+theme(aspect.ratio = 1)

umap_plot_c2 <- ggplot(dfc2, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = color), alpha = 0.4, size = 0.5)+coord_equal(ylim =c(-6.014122,7.909906), xlim=c(-6.096731,9.633184))+
  scale_color_identity()+
  theme_minimal(base_size = 14) +
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")+ ggtitle('HC 2')+theme(panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank())+theme(aspect.ratio = 1)

umap_plot_p15 <- ggplot(dfp15, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = color), alpha = 0.4, size = 0.5) +coord_equal(ylim =c(-6.014122,7.909906), xlim=c(-6.096731,9.633184))+
  scale_color_identity()+
  theme_minimal(base_size = 14) +
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")+ ggtitle('R331P/R331P')+theme(panel.grid.major = element_blank(),
                                                                panel.grid.minor = element_blank())+theme(aspect.ratio = 1)

umap_plot_p16 <- ggplot(dfp16, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = color), alpha = 0.4, size = 0.5)  +coord_equal(ylim =c(-6.014122,7.909906), xlim=c(-6.096731, 9.633184))+
  scale_color_identity()+
  theme_minimal(base_size = 14) +
  guides(x = axis, y = axis)+ xlab('UMAP 1')+ylab('UMAP 2')+
  #theme_set(theme_bw(18))+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(arrow = arrow(type='closed', length = unit(3,'pt'))),
        axis.title = element_text(hjust = 0, size = 12),
        strip.text = element_text(size = 12, color = "black", family = 'arial'),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")+ggtitle("WT/G126D") +theme(panel.grid.major = element_blank(),
                                                             panel.grid.minor = element_blank())+theme(aspect.ratio = 1)

umap_plot_comb$theme$text$family = 'arial'
umap_plot_c1$theme$text$family = 'arial'
umap_plot_c2$theme$text$family = 'arial'
umap_plot_p15$theme$text$family = 'arial'
umap_plot_p16$theme$text$family = 'arial'

umap_plot_comb$theme$text$size = 16
umap_plot_c1$theme$text$size = 16
umap_plot_c2$theme$text$size = 16
umap_plot_p15$theme$text$size = 16
umap_plot_p16$theme$text$size = 16

## Define function for legend
gradient <- function(cell_type, colors, colors_end){
  base_color = colors[cell_type]
  end_color = colors_end[cell_type]
  colorRampPalette(c(base_color, end_color))(10)
}
names <- c(
  "HCs" = "blue",
  "WT/G126D" = "orange",
  "R331P/R331P" = "red"
)
## Create Legend dataframe
legend <- data.frame(
  Patient = rep(names(names), each = 10),
  density = rep(seq(0,1,length.out = 10), length(cluster_colors)),
  color = unlist(lapply(names(cluster_colors), function(ct) gradient(ct, end_colors, cluster_colors)))
)
legendHCs <- legend[legend$Patient  == 'HCs',]
legendHet <- legend[legend$Patient  == "WT/G126D",]
legendHom <- legend[legend$Patient  == "R331P/R331P",]
legend_plot_HCs <- ggplot(legendHCs, aes(x = Patient, y = density, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  labs(title = "Density", x = "", y = "") +
  theme_minimal(base_size = 12, base_family = 'arial') +
  theme(axis.text.x = element_blank(),plot.title = element_text(family = 'arial', size = 12, hjust = .5),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none"
  )
legend_plot_Het <- ggplot(legendHet, aes(x = Patient, y = density, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  labs(title = "Density", x = "", y = "") +
  theme_minimal(base_size = 12, base_family = 'arial') +
  theme(axis.text.x = element_blank(), plot.title = element_text(family = 'arial', size = 12, hjust = .5),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none"
  )
legend_plot_Hom <- ggplot(legendHom, aes(x = Patient, y = density, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  labs(title = "Density", x = "", y = "") +
  theme_minimal(base_size = 12, base_family = 'arial') +
  theme(axis.text.x = element_blank(),plot.title = element_text(family = 'arial', size = 12,
                                                                hjust = .5),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none"
  )


legend_plot_HCs <- legend_plot_HCs + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

legend_plot_Het <- legend_plot_Het + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

legend_plot_Hom <- legend_plot_Hom + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
# Convert the custom legend plot to a grob
legend_grob_HCs <- ggplotGrob(legend_plot_HCs)
legend_grob_Het <- ggplotGrob(legend_plot_Het)
legend_grob_Hom <- ggplotGrob(legend_plot_Hom)

legend_xmin <- 10  # Adjust as needed
legend_xmax <- 12 # Adjust to control width of the legend
legend_ymin <- -9  # Position at the bottom of the plot
legend_ymax <- 9 # Adjust to control height of the legend

# Add the custom legend to the main plot
#combined_plot_c1 <- umap_plot_c1 +
#  annotation_custom(
#    grob = legend_grob_HCs,
#    xmin = legend_xmin,
#    xmax = legend_xmax,
#    ymin = legend_ymin,
#    ymax = legend_ymax
#  )
combined_plot_c2 <- umap_plot_c2 +
  annotation_custom(
    grob = legend_grob_HCs,
    xmin = legend_xmin,
    xmax = legend_xmax,
    ymin = legend_ymin,
    ymax = legend_ymax
  )
combined_plot_p15 <- umap_plot_p15 +
  annotation_custom(
    grob = legend_grob_Hom,
    xmin = legend_xmin,
    xmax = legend_xmax,
    ymin = legend_ymin,
    ymax = legend_ymax)
  
combined_plot_p16 <- umap_plot_p16 +
  annotation_custom(
    grob = legend_grob_Het,
    xmin = legend_xmin,
    xmax = legend_xmax,
    ymin = legend_ymin,
    ymax = legend_ymax
  )

cowplot::plot_grid(umap_plot_comb, NULL, NULL,NULL, NULL, NULL, umap_plot_c1, NULL, NULL, umap_plot_c2, legend_grob_HCs, NULL, umap_plot_p16,legend_grob_Het, NULL, umap_plot_p15, 
                   legend_grob_Hom, labels=NULL, ncol = 6, rel_widths=c(6,.75,.5,6,.75,.5),
                   rel_heights = c(6,6,6))

##########################
### stacked bar plots
##########################

cell_type_proportions <- df %>%
  group_by(new.ident, RNA_snn_res.0.15) %>%
  summarise(count = n()) %>%
  group_by(new.ident) %>%
  mutate(total_count = sum(count),
         proportion = count / total_count)
# check if cluster_proportions is correct
x <- sum(cell_type_proportions$proportion[1:7])
x

bar_plot <- ggplot(cell_type_proportions, aes(x = new.ident, y = proportion, fill = RNA_snn_res.0.15)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(#title = "Distribution of Patient Cells per Cluster",
    x = "",
    y = "Frequency",
    fill = "Cluster") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.location = "panel", axis.ticks.x = element_blank(),
        legend.title = element_text(size = 14, family = 'arial'),
        legend.text = element_text(size = 12, family = 'arial'),
        axis.title.y = element_text(size = 12, family = 'arial'),
        axis.text = element_text(size = 12, family = 'arial'))
bar_plot$theme$text$family = 'arial'

# Display the plot
hcl.colors(n = 25)
safe_colorblind_palette <- c("#117733", "#4B0055", 
                             "#000000", "#888888", "#F0E442",
                             "#CC79A7", "#28E2E5")
colorspace::swatchplot(safe_colorblind_palette, cvd = TRUE)
bar <- bar_plot + scale_fill_manual(breaks = c("0","1","2", "3","4","5","6"), 
                                    values = safe_colorblind_palette[1:7])
bar <- bar + theme(legend.key.size = unit(0.4, "cm"))
bar <- bar+theme(panel.grid.major.x = element_blank())

bar

#######################
## FindAllMarkers
#######################
Idents(B_cells) <- 'RNA_snn_res.0.15'
FindAll_cluster <- FindAllMarkers(B_cells, assay = 'RNA')
Heatmap_cluster <- FindAll_cluster[FindAll_cluster$p_val_adj < 0.05,]

write.csv(FindAll_cluster, 'FindMarkers_RNA_snn_res.0.15.csv')
Heatmap_cluster %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10

genes = unique(top10$gene)
Cluster = unique(top10$cluster)

mat = matrix(0,nrow = length(genes), ncol = length(Cluster))
rownames(mat) = genes
colnames(mat) = Cluster

##Cluster samples and heatmap
for(i in 1:nrow(FindAll_cluster)){mat[rownames(mat)==FindAll_cluster[i,7],colnames(mat)==FindAll_cluster[i,6]] <- FindAll_cluster[i,2]}

library(dendsort)
library(pheatmap)
library(ComplexHeatmap)
# perform hierarchical clustering and 
row_dend = dendsort(hclust(dist(mat), method = "ward.D2"))
col_dend = dendsort(hclust(dist(t(mat)), method = "ward.D2"))

annot <- data.frame(Cluster = c("0", "1", "2", "3", "4", "5","6"))
breaks = seq(-2,2,0.1)
my_colour = list(
  Expression = c('Avg log2FC' = colorRampPalette(c("navy", "white", "red"))(length(breaks))),
  Cluster = c("6"= "#28E2E5", "3" ="#888888","5"= "#CC79A7",'0' = "#117733","4" = "#F0E442", "2" = "#000",'1' = "#4B0055")
)

# Top 10 genes by cluster
ComplexHeatmap::pheatmap(mat,
                                    color=my_colour$Expression, breaks = breaks,
                                    cluster_cols=TRUE, border_color='NA', fontsize_row=12,
                                    treeheight_row=25, treeheight_col=25,name = 'Avg Log2FC',
                                    clustering_method='ward.D2', show_colnames=TRUE, show_rownames = TRUE,
                                    annotation_col=annot, annotation_colors = my_colour, 
                                    fontfamily = 'arial', fontface_row = 'italic')


ComplexHeatmap::pheatmap(mat, breaks = breaks, fontsize_row=2,
                         column_title = NULL,color = my_colour$Expression, labels_col = my_colour$Cluster,
                         border_color='NA',
                         treeheight_row=10, treeheight_col=10,name = 'Avg Log2FC',
                         show_colnames=FALSE,
                         heatmap_legend_param = list(direction = 'horizontal',
                                                     title = "Avg Log2FC", at = c(-2,-1,0,1,2), 
                                                     labels_gp = gpar(col = "black", fontsize = 10, fontfamily='arial'),
                                                     labels = c("-2", "-1", "0","1","2")
                         ),fontfamily = 'arial', fontface_row = 'italic')


heatmap <- ComplexHeatmap::pheatmap(mat, color=my_colour$Expression, breaks = breaks,
                         cluster_cols=TRUE, border_color='black', fontsize_row=9, cluster_rows = TRUE,
                         treeheight_row=50, treeheight_col=50,name = 'Avg Log2FC', 
                         clustering_method='ward.D2', show_colnames=FALSE, show_rownames = TRUE,
                         
                         top_annotation = HeatmapAnnotation(Cluster = anno_block(gp = gpar(fill = my_colour$Cluster),
                                                                             labels = c("6", "3", "5","0","4","2","1"), 
                                                                             labels_gp = gpar(col = "white", fontsize = 12,
                                                                                              fontfamily = 'arial'))),
                        column_split = annot$Cluster,
                                    heatmap_legend_param = list(direction = "horizontal", at =c(-2,-1,0,1,2), legend_height = unit(5, "cm"),title = 'Avg Log2FC', 
                                                                col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                                                                labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 12, fontface = 'bold'),
                                                                labels = c(-2, -1, 0,1,2), border = TRUE
                                    ),fontfamily = 'arial', fontface_row = 'italic') 


draw(heatmap, heatmap_legend_side= "bottom", padding = unit(c(2, 2, 2, 15), "mm"))



