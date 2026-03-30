
getwd()
setwd("/local/workdir/cqr3/human_FT_samples/All_Healthy_Human_FT/")
getwd()


#### Load Libraries and Functions ####
library(dplyr) # Click Manually - 1.1.0
library(ggplot2) # Click Manually - 3.4.0
library(patchwork)
library(Seurat) # Click Manually? - 4.3.0
library(harmony)

library(cowplot)
library(SoupX)
library(DoubletFinder)
library(data.table)
library(parallel)
library(tidyverse) # Check Later: Error w/ libiconv
library(SoupX)
library(DropletUtils)

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(viridisLite)
library(Polychrome)
library(circlize)
library(NatParksPalettes)

library(phateR)
library(monocle3)

#### Load Objects ####
Human_HQ_Cells <- readRDS(file = "20241212_FT_Healthy_Cells_Final.rds", refhook = NULL)

FT_Named <- RenameIdents(Human_HQ_Cells, 
                         '0' = "T/NK Cell 1", 
                         '1' = "Fibroblast 1", 
                         '2' = "Secretory Epithelial 1",
                         '3' = "Fibroblast 2", 
                         '4' = "T/NK Cell 2", 
                         '5' = "T/NK Cell 3", 
                         '6' = "Macrophage",
                         '7' = "Ciliated Epithelial",
                         '8' = "Blood Endothelial", 
                         '9' = "Pericyte", 
                         '10' = "Smooth Muscle 1", 
                         '11' = "Lymph Endothelial", 
                         '12' = "Mast", 
                         '13' = "Secretory Epithelial 2", 
                         '14' = "B Cell", 
                         '15' = "Cycling Immune", 
                         '16' = "T/NK Cell 4",
                         '17' = "T/NK Cell 5",
                         '18' = "Smooth Muscle 2",
                         '19' = "Secretory Epithelial 3")


FT_Named@active.ident <- factor(x = FT_Named@active.ident, 
                                levels = 
                                  c("Fibroblast 1",
                                    "Fibroblast 2", 
                                    "Smooth Muscle 1", 
                                    "Smooth Muscle 2", 
                                    "Pericyte",
                                    "Blood Endothelial",
                                    "Lymph Endothelial",
                                    "Secretory Epithelial 1", 
                                    "Secretory Epithelial 2", 
                                    "Secretory Epithelial 3",
                                    "Ciliated Epithelial",
                                    "T/NK Cell 1", 
                                    "T/NK Cell 2",
                                    "T/NK Cell 3",
                                    "T/NK Cell 4",
                                    "T/NK Cell 5",
                                    "Cycling Immune", 
                                    "B Cell", 
                                    "Macrophage", 
                                    "Mast"))

Mouse_HQ_Cells <- readRDS(file = "20220810_Filtered_Distal_MO_Cells.rds" , refhook = NULL)

UT_Named <- RenameIdents(Mouse_HQ_Cells, 
                             '0' = "Fibroblast 1", 
                             '1' = "Fibroblast 2", 
                             '2' = "Secretory Epithelial",
                             '3' = "Smooth Muscle", 
                             '4' = "Ciliated Epithelial 1", 
                             '5' = "Fibroblast 3", 
                             '6' = "Stem-like Epithelial 1",
                             '7' = "Stem-like Epithelial 2",
                             '8' = "Ciliated Epithelial 2", 
                             '9' = "Blood Endothelial", 
                             '10' = "Pericyte", 
                             '11' = "Intermediate Epithelial", 
                             '12' = "T/NK Cell", 
                             '13' = "Epithelial/Fibroblast", 
                             '14' = "Macrophage", 
                             '15' = "Erythrocyte", 
                             '16' = "Luteal",
                             '17' = "Mesothelial",
                             '18' = "Lymph Endothelial/Epithelial") # Remove cluster due few data points and suspected doublet


UT_Named@active.ident <- factor(x = UT_Named@active.ident, 
                                    levels = c('Fibroblast 1',
                                               'Fibroblast 2',
                                               'Fibroblast 3',
                                               'Smooth Muscle',
                                               'Pericyte',
                                               'Blood Endothelial',
                                               'Lymph Endothelial/Epithelial',
                                               'Epithelial/Fibroblast',
                                               'Stem-like Epithelial 1',
                                               'Stem-like Epithelial 2',
                                               'Intermediate Epithelial',
                                               'Secretory Epithelial',
                                               'Ciliated Epithelial 1',
                                               'Ciliated Epithelial 2',
                                               'T/NK Cell',
                                               'Macrophage',
                                               'Erythrocyte',
                                               'Mesothelial',
                                               'Luteal'))

UT_Named <- SetIdent(UT_Named, value = UT_Named@active.ident)



UT_Named <- subset(UT_Named, 
                       idents = c('Fibroblast 1',
                                  'Fibroblast 2',
                                  'Fibroblast 3',
                                  'Smooth Muscle',
                                  'Pericyte',
                                  'Blood Endothelial',
                                  'Epithelial/Fibroblast',
                                  'Stem-like Epithelial 1',
                                  'Stem-like Epithelial 2',
                                  'Intermediate Epithelial',
                                  'Secretory Epithelial',
                                  'Ciliated Epithelial 1',
                                  'Ciliated Epithelial 2',
                                  'T/NK Cell',
                                  'Macrophage',
                                  'Erythrocyte',
                                  'Mesothelial',
                                  'Luteal'))


#### Fig 1C ####

Fibroblasts <- c('#FF9D00' , '#FFB653' , '#FFCB9A')   # Oranges
Muscle <- c('#E55451' , '#FFB7B2') # Reds
Endothelial <- c('#8D021F')  # Reds
FiboEpi <- "#ffad77" # Reddish Brown
Epi <-c('#6E3E6E','#8A2BE2','#604791','#CCCCFF','#DA70D6','#DF73FF') # Blues/Purples
Immune <- c( "#238B45", "#7FFFD4", '#FC86AA') # Yellowish Brown
Meso <- 'darkgrey' # Pink
Lut <- "#9DCC00" # Green

colors <- c(Fibroblasts, Muscle, Endothelial, FiboEpi, Epi, Immune, Meso, Lut)

p1 <- DimPlot(
  UT_Named,
  reduction='umap',
  cols=colors,
  pt.size = 0.5,
  label.size = 4,
  label.color = "black",
  repel = TRUE,
  label=F) +
  NoLegend() +
  labs(x="UMAP_1",y="UMAP_2")

LabelClusters(p1, id="ident", color = "black", repel = T , size = 4, box.padding = .75)

ggsave(filename = "20260316_FIG1b_all_distal_umap.pdf", plot = p1, width = 12, height = 8, dpi = 600)


#### Fig 1D ####


oranges <- c("#FEB24C", "#FD8D3C")
deeper_reds <- c("#A50F15", "#67000D", "#FC4E2A", "#E31A1C")
green <- "#06402b"
blues <- c("#3182BD", "#08519C", "#6BAED6",  "#6A5ACD")
greens <- c("#005824", "#238B45", "#41AE76", "#66C2A4", 
            "#228B22","#98fb98", "#2CFF05" , "#7FFFD4",
            "#CCFF00")


# Combine the colors into one vector
colors <- c(oranges, deeper_reds, green, blues, greens)

all_human <- DimPlot(object = FT_Named, 
                     reduction = 'umap', 
                     group.by = "ident",
                     repel = TRUE, 
                     label = F, 
                     pt.size = 0.01, 
                     label.size = 5,
                     raster = F) + 
  labs(title = 'Colored by Cluster')


ggsave(filename = "20250122_all_human_UMAP.pdf" , plot = all_human , width = 15, height = 10, dpi =600)


#### Fig 1F ####



mo_features <- c("Spp1","Fkbp5",
                 "Msln","Lrrn4",   
                 "Kit", 
                 "Hbb-bs","Jchain" ,"Aif1",
                 "Mki67","Klrd1","Cd4", "Cd8a",
                 "Ptprc", 
                 "Foxj1",  "Ovgp1", "Pax8",
                 "Krt8", "Epcam",
                 "Lyve1", "Vwf",
                 "Mcam",
                 "Acta2","Myh11", 
                 "Dcn","Col1a1")               


hu_features <- c( "SPP1","FKBP5",
               "MSLN","LRRN4",   
               "KIT", 
               "HBB","JCHAIN" ,"AIF1",
               "MKI67","KLRD1","CD4", "CD8A",
               "PTPRC", 
               "FOXJ1",  "OVGP1", "PAX8",
               "KRT8", "EPCAM",
               "LYVE1", "VWF",
               "MCAM",
               "ACTA2","MYH11", 
               "DCN","COL1A1")              

mo_dp <- DotPlot(object = UT_Named,                    # Seurat object
                  assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                  features = mo_features,                  # List of features (select one from above or create a new one)
                  cols =  c('grey','#FFA500'),
                  col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
                  col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
                  dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
                  dot.scale = 9,                        # Scale the size of the points
                  group.by = NULL,              # How the cells are going to be grouped
                  split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
                  scale = TRUE,                         # Whether the data is scaled
                  scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
                  scale.min = NA,                       # Set lower limit for scaling
                  scale.max = NA                        # Set upper limit for scaling
)+    
  labs(x = NULL, y = NULL)+
  #scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)
  #scale_fill_continuous(values=c("grey", '#785AA3'))+
  #scale_fill_gradientn(colors=c("#1984c5", "#e2e2e2", "#c23728"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  #theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 14 , face = "italic"))+
  theme(axis.text.y = element_text(size = 14))+
  scale_y_discrete(limits = levels(UT_Named))+
  theme(legend.title = element_text(size = 14))+
  coord_flip()


ggsave(filename = "20260311_FIG1f_all_distal_mouse_dotplot.pdf", plot = mo_dp, width = 10, height = 10, dpi = 600)




hu_dp <- DotPlot(object = FT_Named,                    # Seurat object
                  assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                  features = features,                  # List of features (select one from above or create a new one)
                  cols =  c('grey','#785AA3'),
                  col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
                  col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
                  dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
                  dot.scale = 9,                        # Scale the size of the points
                  group.by = NULL,              # How the cells are going to be grouped
                  split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
                  scale = TRUE,                         # Whether the data is scaled
                  scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
                  scale.min = NA,                       # Set lower limit for scaling
                  scale.max = NA                        # Set upper limit for scaling
)+    
  labs(x = NULL, y = NULL)+
  #scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)
  #scale_fill_continuous(values=c("grey", '#785AA3'))+
  #scale_fill_gradientn(colors=c("#1984c5", "#e2e2e2", "#c23728"))+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  #theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 14 , face = "italic"))+
  theme(axis.text.y = element_text(size = 14))+
  scale_y_discrete(limits = levels(FT_Named))+
  theme(legend.title = element_text(size = 14))+
  coord_flip()





ggsave(filename = "FIG1f_all_distal_human_dotplot.pdf", plot = hu_dp, width = 10, height = 10, dpi = 600)




