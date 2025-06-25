##########################
### load in the necessary libraries and set the working directory
##########################
install.packages("extrafont")
library(extrafont)
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

getwd()
setwd("/Volumes/pmi/Data/Research/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/11725-JM")


##########################
### Predicted cell type level 1 scores boxplots
##########################
object <- readRDS('/Volumes/Markle Lab/Data/Christine Mariskanish/11725-JM experiment/CARD11/11725-JM/prefiltered_object.rds')

object$newer.ident <- NA
object$newer.ident[which(object$orig.ident == "0008")] <- "HC"
object$newer.ident[which(object$orig.ident == "0009")] <- "HC"
object$newer.ident[which(object$orig.ident == "0015")] <- "Homozygote"
object$newer.ident[which(object$orig.ident == "0016")] <- "Heterozygote"
object$newer.ident <- factor(x = object$newer.ident, levels = c("HC", "Heterozygote", "Homozygote"))
object$new.ident <- factor(x = object$new.ident, levels = c("HC 1", "HC 2", "WT/G126D", "R331P/R331P"))

df <- data.frame(object$newer.ident)
df <- tibble::rownames_to_column(df, "newer.ident")
score1 <- data.frame(object$predicted.celltype.l1.score)
score1 <- tibble::rownames_to_column(score1, "newer.ident")
score2 <- data.frame(object$predicted.celltype.l2.score)
score2 <- tibble::rownames_to_column(score2, "newer.ident")
sample <- data.frame(object$new.ident)
sample <- tibble::rownames_to_column(sample, "newer.ident")
df <- merge(score1, df, by = "newer.ident")
df <- merge(score2, df, by = "newer.ident")
df<- merge(df, sample, by = "newer.ident")

names <- c(
  "HC" = "blue",
  "Heterozygote" = "orange",
  "Homozygote" = "red"
)

boxplot.l1 <- df %>%
  ggplot( aes(x=factor(object.new.ident), y=object.predicted.celltype.l1.score, fill=factor(object.newer.ident))) +
  geom_violin()+ geom_point(size = 0)+
  geom_boxplot(width=0.25, color="black", alpha=0.01,outlier.colour="blue", outlier.size=0.1, outlier.shape = NA)+
  scale_fill_manual(values =names) +
  theme_bw()+coord_cartesian(ylim = c(0,1)) +
  theme(text=element_text(family="arial", size=12),
        legend.position="none"
  ) +
  xlab("") + ylab("Cell Type Level 1 Prediction Score")+ 
  theme(panel.grid.major.x = element_blank())+
  geom_hline(yintercept = (mean(df$object.predicted.celltype.l1.score) - sd(df$object.predicted.celltype.l1.score)), linetype = 2)
boxplot.l1


##########################
### Predicted cell type level 2 scores boxplot
##########################

boxplot.l2 <- df %>%
  ggplot( aes(x=object.new.ident, y=object.predicted.celltype.l2.score, fill=object.newer.ident)) +
  geom_violin() + geom_point(size = 0)+geom_boxplot(width=0.25, color="black", alpha=0.01, 
                                                    outlier.colour="blue", outlier.size=0.1, outlier.shape = NA)+
  scale_fill_manual(values=names) +
  theme_bw()+coord_cartesian(ylim = c(0,1)) +
  theme(text=element_text(family="arial", size=12),
        legend.position="none"
  ) +
  xlab("") + ylab("Cell Type Prediction Level 2 Score")+ theme(panel.grid.major.x = element_blank())#+

boxplot.l1|boxplot.l2


boxplot.l1 + geom_text(x = 1.5, y= 0, 
                       label = paste("Kruskal-Wallis p = < 0.05"),
                       size = 5, family = 'arial')
cowplot::plot_grid(boxplot.l1 , boxplot.l2, ncol = 1, nrow = 2, rel_widths=c(7,7),
                   rel_heights = c(7,7))

## PangalaoDb
Idents(object.sct.filtered) <- 'predicted.celltype.l1'

breakpoints <- c(-4, seq(-4, 4, 0.2))
colors <- c("red", RColorBrewer::brewer.pal(9, "Blues"))

dot <- scCustomize::DotPlot_scCustom(object.sct.filtered, assay = 'RNA',
                                     features = c("CD4", "CD27", "CD28", "CD3D", "CD3E","CD8A", "CD8B", # T cells
                                                  "CXCR6", "KLRB1", "NKG7", # NK
                                                  "CD14",  #Monocyte and cd16, cd56 = NK
                                                  "IL3RA", "CD86", #DC
                                                  "MS4A1", "CD79A", "CD79B", "CD19", "CD38" # B cells
                                     ), group.by = "predicted.celltype.l1", 
                                     x_lab_rotate = TRUE, remove_axis_titles = FALSE)
dot$theme$text$family = 'arial'
dot <- dot + theme(legend.key.size = unit(0.4, "cm"))
dot$labels$y = 'Cell Type Level 1'
dot$labels$x = 'Marker Genes'
dot$theme$axis.title.x$face = 'bold'
dot$theme$axis.title.y$face = 'bold'
dot$theme$text$size = 12
dot$theme$axis.text$size = 10


