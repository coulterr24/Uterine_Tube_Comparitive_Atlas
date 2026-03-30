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

library(UCell)
library(tidyr)



##### Load All Objects ####

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



Mouse_HQ_Cells <- readRDS(file = "20220805_Mouse_HQ_Cells.rds" , refhook = NULL)

UT_Named <- RenameIdents(Mouse_HQ_Cells, 
                         '0' = "Fibroblast 1", 
                         '1' = "Fibroblast 2", 
                         '2' = "Secretory Epithelial",
                         '3' = "Stem-like Epithelial", 
                         '4' = "Fibroblast 3", 
                         '5' = "Smooth Muscle", 
                         '6' = "Ciliated Epithelial 1",
                         '7' = "Ciliated Epithelial 2",
                         '8' = "Blood Endothelial", 
                         '9' = "T-Cell", 
                         '10' = "Pericyte", 
                         '11' = "Epithelial/Fibroblast", 
                         '12' = "Macrophage", 
                         '13' = "Mesothelial", 
                         '14' = "Erythrocyte", 
                         '15' = "Lymphatic Endothelial", 
                         '16' = "Luteal")

UT_Named@active.ident <- factor(x = UT_Named@active.ident, 
                                levels = 
                                  rev(c("Luteal",
                                        "Erythrocyte",
                                        "Macrophage",
                                        "T-Cell",
                                        "Mesothelial",
                                        "Epithelial/Fibroblast",
                                        "Stem-like Epithelial",
                                        "Secretory Epithelial",
                                        "Ciliated Epithelial 2",
                                        "Ciliated Epithelial 1",
                                        "Lymphatic Endothelial",
                                        "Blood Endothelial",
                                        "Pericyte",
                                        "Smooth Muscle", 
                                        "Fibroblast 3", 
                                        "Fibroblast 2", 
                                        "Fibroblast 1")))

#### Supp Figure 4A ####


Fibroblasts <- c('#FF9D00' , '#FFB653' , '#FFCB9A')   # Oranges
Muscle <- c('#E55451' , '#FFB7B2') # Reds
Endothelial <- c('#8D021F')  # Reds
FiboEpi <- "#ffad77" # Reddish Brown
Epi <-c('#6E3E6E','#8A2BE2','#604791','#CCCCFF','#DA70D6','#DF73FF') # Blues/Purples
Immune <- c( "#238B45", "#7FFFD4", '#FC86AA') # Yellowish Brown
Meso <- 'darkgrey' # Pink
Lut <- "#9DCC00" # Green

colors <- c(Fibroblasts, Muscle, Endothelial, FiboEpi, Epi, Immune, Meso, Lut)


table <- table(UT_Named@meta.data$Location ,
               UT_Named@active.ident)    # Create a table of counts

IDs = c("Fibroblast 1",
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
        "T Cell 1", 
        "T Cell 2",
        "T Cell 3",
        "T Cell 4",
        "T Cell 5",
        "Cycling Immune", 
        "B Cell", 
        "Macrophage", 
        "Mast")

df <- data.frame(table) 

table2 <- table(loc@active.ident)


age <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
              aes(x = Var1,              # Variable to plot on the x-axis
                  y = Freq,              # Variable to plot on the y-axis
                  fill = factor(Var2,    # Variable to fill the bars
                                levels = IDs))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", IDs,
                    values = hu_colors)+
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis


ggsave(filename = "FigS4A_location_mouse_cell_dist.pdf", plot = age, width = 10, height = 14, dpi = 600)



#### Supp Figure 4B ####


oranges <- c("#FEB24C", "#FD8D3C")
deeper_reds <- c("#A50F15", "#67000D", "#FC4E2A", "#E31A1C")
green <- "#06402b"
blues <- c("#3182BD", "#08519C", "#6BAED6",  "#6A5ACD")
greens <- c("#005824", "#238B45", "#41AE76", "#66C2A4", 
            "#228B22","#98fb98", "#2CFF05" , "#7FFFD4",
            "#CCFF00")


# Combine the colors into one vector
colors <- c(oranges, deeper_reds, green, blues, greens)


table <- table(FT_Named@meta.data$Location ,
               FT_Named@active.ident)    # Create a table of counts

IDs = c("Luteal",
        "Erythrocyte",
        "Macrophage",
        "T-Cell",
        "Mesothelial",
        "Epithelial/Fibroblast",
        "Stem-like Epithelial",
        "Secretory Epithelial",
        "Ciliated Epithelial 2",
        "Ciliated Epithelial 1",
        "Lymphatic Endothelial",
        "Blood Endothelial",
        "Pericyte",
        "Smooth Muscle", 
        "Fibroblast 3", 
        "Fibroblast 2", 
        "Fibroblast 1")

df <- data.frame(table) 


#df <- na.omit(df)


#table2 <- table(FT_Named@active.ident)


source <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
                 aes(x = Var1,              # Variable to plot on the x-axis
                     y = Freq,              # Variable to plot on the y-axis
                     fill = factor(Var2,    # Variable to fill the bars
                                   levels = IDs))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", IDs,
                    values = colors)+
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis

ggsave(filename = "202501125_cell_human_Location_dist.pdf", plot = source, width = 10, height = 14, dpi = 600)













## Stacked Bar Location #

table <- table(Epi_Filter2@meta.data$Location ,
               Epi_Filter2@active.ident)    # Create a table of counts

IDs = unique(Epi_Filter2@active.ident)

df <- data.frame(table) 

table2 <- table(
  
  y <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
              aes(x = Var1,              # Variable to plot on the x-axis
                  y = Freq,              # Variable to plot on the y-axis
                  fill = factor(Var2,    # Variable to fill the bars
                                levels = c("KRT7+ Secretory", 
                                           "PAX8+/LGR5+ Progenitor",
                                           "TEX14+ Progenitor",
                                           "HPSE2+ Progenitor",
                                           "APOA1+ Progenitor",
                                           "KRT5+ Progenitor",
                                           "AGBL4+/MAPK8+ Pre-Ciliated", 
                                           "CFAP299+ Ciliated",
                                           "TPPP3+ Ciliated")),
              )) + # Order of the stacked bars
    theme_classic() +               # ggplot2 theme
    # Bar plot
    geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
             stat = "identity",     # Height of bars represent values in the data
             size = 1) +            # Size of bars
    # Color scheme
    scale_fill_manual("Location", c("KRT7+ Secretory",
                                    "PAX8+/LGR5+ Progenitor",
                                    "TEX14+ Progenitor",
                                    "HPSE2+ Progenitor",
                                    "APOA1+ Progenitor",
                                    "KRT5+ Progenitor",
                                    "AGBL4+/MAPK8+ Pre-Ciliated", 
                                    "CFAP299+ Ciliated",
                                    "TPPP3+ Ciliated"),
                      values = c("#B20224", #1,
                                 '#35EFEF',
                                 "#F28D86", #2
                                 "#2188F7", #3
                                 "#59D1AF", #5
                                 '#00A1C6',
                                 "#EA68E1", #4
                                 "#A374B5", #8
                                 "#9000C6"))+  
    labs(x = NULL,                     # x-axis label
         y = "Fraction of Cells") +    # y-axis label
    theme(text = element_text(size = 15),                                      # Text size throughout the plot
          axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
          axis.text.y = element_text(color = 'black', hjust = 1))    )         # Text color and horizontal adjustment on y-axis
  
  ggsave(filename = "20251020_epi_pre-cilia_location_stacked_bar.pdf" , plot = y , width = 10, height = 12, dpi =600)
  
#### Supp Figure 4C ####

  
  FT_sub <- subset(FT_sub, idents = c("Macrophage"))
  
  
  
  FT_sub$Location <- factor(FT_sub$Location, levels = c('Isthmus', 'Ampulla', 'Fimbria'))
  
  
  
  Mphage_CCL20 <- VlnPlot(FT_sub, features = "CCL20", 
                          group.by = 'Location', 
                          #slot = "counts",
                          pt.size = 0)
  
  ggsave(filename = "202501027_human_CCL20_MPhage_cells_violin.pdf" , plot = Mphage_CCL20 , width = 15, height = 10, dpi =600)
  
  
 
  
  
  
  
  
#### Supp Figure 4D ####
  
  
  
  
  UT_sub <- subset(UT_Named, idents = c("Macrophage"))
  
  
  
  UT_sub$Location <- factor(UT_sub$Location, levels = c('Proximal', 'Distal'))
  
  
  
  Mphage_Ccl20 <- VlnPlot(UT_sub, features = "Ccl20", 
                          group.by = 'Location', 
                          #slot = "counts",
                          pt.size = 0)
  
  ggsave(filename = "202501027_mouse_Ccl20_MPhage_cells_violin.pdf" , plot = Mphage_Ccl20 , width = 15, height = 10, dpi =600)
  
  
#### Supp Figure 4E ####
  
  
  
  
  Fibro_sub <- subset(FT_Named, idents = c("Fibroblast 1",
                                           "Fibroblast 2"))
  
  
  
  Idents(Fibro_sub) <-  "Location"
  
  Fibro_sub <- subset(Fibro_sub, idents = c("Isthmus",
                                            "Ampulla",
                                            "Fimbria"))
  
  Fibro_sub@active.ident <- factor(Fibro_sub$Location, levels = c("Isthmus",
                                                                  "Ampulla",
                                                                  "Fimbria"))
  
  
  genes <- c('KCND2',
             'POSTN',
             'KCNIP4',
             'NEGR1',
             'SCARA5',
             'CFD',
             'NOVA1',
             'GREB1L',
             'MYH11'
  )
  
  
  fibro_reg_dp <- DotPlot(object = Fibro_sub,                    # Seurat object
                          assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                          features = genes,                 # List of features (select one from above or create a new one)
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
    scale_y_discrete(limits = levels(Fibro_sub))+
    theme(legend.title = element_text(size = 14))+
    coord_flip()
  
  
  ggsave(filename = "20260120_Reg_Fibrobalst_Dot_Plot_HU.pdf", plot = fibro_reg_dp, width = 6, height = 8, dpi = 600)
  
  
  
  
  
#### Supp Figure 4F ####
  
  
  
  
  Macro_sub <- subset(FT_Named, idents = c("Macrophage"))
  
  
  
  Idents(Macro_sub) <-  "Location"
  
  Macro_sub <- subset(Macro_sub, idents = c("Isthmus",
                                            "Ampulla",
                                            "Fimbria"))
  
  Macro_sub@active.ident <- factor(Macro_sub$Location, levels = c("Isthmus",
                                                                  "Ampulla",
                                                                  "Fimbria"))
  
  
  genes <- c('CCL20',
             'TIMP1',
             'LGALS1',
             'LIMK2',
             'ADGRG3',
             'FCGR3B',
             'CD163',
             'C1QB',
             'MRC1'
  )
  
  
  macro_reg_dp <- DotPlot(object = Macro_sub,                    # Seurat object
                          assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                          features = genes,                 # List of features (select one from above or create a new one)
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
    scale_y_discrete(limits = levels(Macro_sub))+
    theme(legend.title = element_text(size = 14))+
    coord_flip()
  
  
  ggsave(filename = "20260120_Reg_Macrophage_Dot_Plot_HU.pdf", plot = macro_reg_dp, width = 6, height = 8, dpi = 600)
  
  
#### Supp Figure 4G ####

  ## Epi Subset
  
  
  Epi_Filter <- subset(FT_Named, 
                       idents = c("Secretory Epithelial 1", 
                                  "Secretory Epithelial 2",
                                  "Secretory Epithelial 3",
                                  "Ciliated Epithelial"))
  
  
  
  
  
  # Need to Normalize and Run PCA before Running Harmony ##
  
  Epi_Filter <- NormalizeData(object = Epi_Filter, assay = 'RNA')
  Epi_Filter <- FindVariableFeatures(object = Epi_Filter, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
  Epi_Filter <- ScaleData(object = Epi_Filter, assay = 'RNA')
  Epi_Filter <- RunPCA(object = Epi_Filter, assay = 'RNA', features = VariableFeatures(object = Epi_Filter),
                       reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')
  
  
  # Calculate the number of PCs that contain some proportion (95%) of the variance
  npcs <- function(seu, var.toal=0.95, reduction="pca"){
    if(is.null(seu@reductions[[reduction]])){
      cat("Reduction", reduction, "not found!")
      return(NULL)
    }
    tmp.var <- (seu@reductions[[reduction]]@stdev)^2
    var.cut <- var.toal*sum(tmp.var)
    n.pcs=0
    var.sum = 0
    while(var.sum < var.cut){
      n.pcs = n.pcs + 1
      var.sum <- var.sum + tmp.var[n.pcs]
    }
    return(n.pcs)
  }
  
  pcs <- npcs(seu = Epi_Filter, var.toal = 0.95, reduction = 'pca_RNA')
  
  ## Run Harmony ##
  
  gc()
  
  Epi_Filter <- RunHarmony(Epi_Filter, 
                           group.by.vars = c("Sample.ID"), 
                           reduction = "pca_RNA", 
                           assay = "RNA", 
                           plot_convergence = TRUE)
  
  ## Normal Seurat clustering ##
  
  Epi_Filter <- FindNeighbors(object = Epi_Filter, 
                              reduction = "harmony",
                              k.param = 70, 
                              dims = 1:pcs, 
                              graph.name = 'harmony_snn')
  
  Epi_Filter <- FindClusters(object = Epi_Filter, 
                             resolution = 0.4, 
                             graph.name = 'harmony_snn')
  
  
  ## Set Python UMAP via reticulate ##
  
  umap.method = 'umap-learn'
  metric = 'correlation'
  
  Epi_Filter <- RunUMAP(object = Epi_Filter, 
                        reduction = "harmony", dims = 1:pcs)
  
  
  DimPlot(object = Epi_Filter, 
          reduction = 'umap', 
          group.by = "seurat_clusters", 
          repel = TRUE, 
          label = TRUE, 
          pt.size = 0.5, 
          label.size = 5,
          raster = F) + 
    labs(title = 'Colored by Cluster')
  
  
  DimPlot(object = Epi_Filter, 
          reduction = 'umap', 
          group.by = "Location", 
          repel = TRUE, 
          label = TRUE, 
          pt.size = 0.5, 
          label.size = 5,
          raster = F) + 
    labs(title = 'Colored by Cluster')
  

  
  EpiMarkers <- FindAllMarkers(Epi_Filter, 
                               max.cells.per.ident = 5000)
  
  
  
  
  # Remove LQ Cells - 4(Fibroblast) / 6(Basal) / 9 (Endothelial) / 10 (MPhage-like) / 11 (Blood-Like)
  
  Epi_Filter2 <- subset(Epi_Filter, 
                        idents = c("0", 
                                   "1",
                                   "2",
                                   "3",
                                   "5",
                                   "7",
                                   '8'))
  
  ## Reprocess 2nd Object ##
  # Need to Normalize and Run PCA before Running Harmony ##
  
  Epi_Filter2 <- NormalizeData(object = Epi_Filter2, assay = 'RNA')
  Epi_Filter2 <- FindVariableFeatures(object = Epi_Filter2, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
  Epi_Filter2 <- ScaleData(object = Epi_Filter2, assay = 'RNA')
  Epi_Filter2 <- RunPCA(object = Epi_Filter2, assay = 'RNA', features = VariableFeatures(object = Epi_Filter2),
                        reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')
  
  
  pcs <- npcs(seu = Epi_Filter2, var.toal = 0.95, reduction = 'pca_RNA')
  
  ## Run Harmony ##
  
  gc()
  
  Epi_Filter2 <- RunHarmony(Epi_Filter2, 
                            group.by.vars = c("Sample.ID"), 
                            reduction = "pca_RNA", 
                            assay = "RNA", 
                            plot_convergence = TRUE)
  
  ## Normal Seurat clustering ##
  
  Epi_Filter2 <- FindNeighbors(object = Epi_Filter2, 
                               reduction = "harmony",
                               k.param = 60, 
                               dims = 1:pcs, 
                               graph.name = 'harmony_snn')
  
  Epi_Filter2 <- FindClusters(object = Epi_Filter2, 
                              resolution = 0.8, 
                              graph.name = 'harmony_snn')
  
  
  ## Set Python UMAP via reticulate ##
  
  umap.method = 'umap-learn'
  metric = 'correlation'
  
  Epi_Filter2 <- RunUMAP(object = Epi_Filter2, 
                         reduction = "harmony", dims = 1:pcs)
  
  
  DimPlot(object = Epi_Filter2, 
          reduction = 'umap', 
          group.by = "seurat_clusters", 
          repel = TRUE, 
          label = TRUE, 
          pt.size = 0.5, 
          label.size = 5,
          raster = F) + 
    labs(title = 'Colored by Cluster')
  
  DimPlot(object = Epi_Filter2, 
          reduction = 'umap', 
          group.by = "Location", 
          repel = TRUE, 
          label = TRUE, 
          pt.size = 0.5, 
          label.size = 5,
          raster = F) + 
    labs(title = 'Colored by Cluster')

  
  saveRDS(Epi_Filter2, file = "20250702_AvF_Epi_Cells.rds")
  
  
  ### Annotation of Epi Cells ###
  
  
  library(hash)
  
  Epi_Filter2 <- readRDS(file = "20250702_AvF_Epi_Cells.rds", refhook = NULL)
  
  Epi_Filter2$cluster <- Idents(Epi_Filter2)
  
  Fimb_Epi_Filter <- readRDS(file = "20250123_Distal_Epi_Cells.rds", refhook = NULL)
  
  
  Fimb_Epi_Filter <- RenameIdents(Fimb_Epi_Filter, 
                                  '0' = "PAX8+/LGR5+ Progenitor", 
                                  '1' = "KRT7+ Secretory", 
                                  '2' = "CFAP299+ Ciliated",
                                  '3' = "TPPP3+ Ciliated", 
                                  '4' = "KRT5+ Progenitor", 
                                  '5' = "APOA1+ Progenitor", 
                                  '6' = "TEX14+ Progenitor",
                                  '7' = "AGBL4+/MAPK8+ Pre-Ciliated",
                                  '8' = "HPSE2+ Progenitor")
  
  
  Fimb_Epi_Filter@active.ident <- factor(x = Fimb_Epi_Filter@active.ident, 
                                         levels = 
                                           c("KRT7+ Secretory", 
                                             "TEX14+ Progenitor",
                                             "HPSE2+ Progenitor",
                                             "KRT5+ Progenitor",
                                             "APOA1+ Progenitor",
                                             "PAX8+/LGR5+ Progenitor",
                                             "AGBL4+/MAPK8+ Pre-Ciliated", 
                                             "CFAP299+ Ciliated",
                                             "TPPP3+ Ciliated"))
  
  ncol(Fimb_Epi_Filter) #15941 cells
  #View(epi)
  Idents(Fimb_Epi_Filter) <- Fimb_Epi_Filter@active.ident
  DimPlot(Fimb_Epi_Filter, reduction = "umap")
  
  
  
  
  #Duplicate net object seurat clusters metadata column###
  Celltype_IDs <- as.character(Epi_Filter2@active.ident)
  Epi_Filter2$combined_clusters <- Celltype_IDs
  
  
  ##epithelial cells ##
  epidf <- data.frame(values=Fimb_Epi_Filter@active.ident)
  epidf[] <- lapply(epidf, as.character)
  
  #Match barcodes and replace cell type label based on barcode
  epidf_hash <- hash(keys=rownames(epidf), values=epidf$values) #create 'dictionary' of barcodes and cell type ID
  
  for(i in 1:nrow(Epi_Filter2@meta.data)){     #for loop to match barcodes, and if they do, replace with new cell type ID
    if((length(epidf_hash[[colnames(Epi_Filter2)[i]]]) == 1)==TRUE){
      Epi_Filter2@meta.data[i,25] <- epidf_hash[[colnames(Epi_Filter2)[i]]]   #[i,n of metadata]
    }
  }
  
  View(Epi_Filter2$combined_clusters)
  
  DimPlot(Epi_Filter2, group.by = 'combined_clusters',
          cols = c('grey','grey','grey','grey','grey','grey',
                   'grey','grey','grey','grey', 'grey',
                   "#B20224", #1
                   "#F28D86", #2
                   "#2188F7", #3
                   '#00A1C6',
                   "#59D1AF", #5
                   '#35EFEF',
                   "#EA68E1", #4
                   "#A374B5", #8
                   "#9000C6"),
          pt.size = .9,
          shuffle = T) #9)
  
  FeaturePlot(Epi_Filter2, features = c("KRT5"), pt.size = 0.75 , raster = F)+
    scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
    theme(plot.title = element_text(size = 32,face = "bold.italic"))
  
  
  Idents(Epi_Filter2) <- Epi_Filter2$combined_clusters
  
  Epi_Filter2 <- RenameIdents(Epi_Filter2, 
                              '0' = "KRT7+ Secretory", 
                              '1' = "PAX8+/LGR5+ Progenitor", 
                              '2' = "CFAP299+ Ciliated",
                              '3' = "TPPP3+ Ciliated", 
                              '4' = "APOA1+ Progenitor", 
                              '5' = "TEX14+ Progenitor", 
                              '6' = "KRT5+ Progenitor",
                              '7' = "KRT5+ Progenitor",                               
                              '8' = "HPSE2+ Progenitor",
                              '9' = "AGBL4+/MAPK8+ Pre-Ciliated",
                              '10' = "PAX8+/LGR5+ Progenitor", 
                              '11' = "TPPP3+ Ciliated",
                              'PAX8+/LGR5+ Progenitor' = "PAX8+/LGR5+ Progenitor", 
                              'KRT7+ Secretory' = "KRT7+ Secretory", 
                              'CFAP299+ Ciliated' = "CFAP299+ Ciliated",
                              'TPPP3+ Ciliated' = "TPPP3+ Ciliated", 
                              'KRT5+ Progenitor' = "KRT5+ Progenitor", 
                              'APOA1+ Progenitor' = "APOA1+ Progenitor", 
                              'TEX14+ Progenitor' = "TEX14+ Progenitor",
                              'AGBL4+/MAPK8+ Pre-Ciliated' = "AGBL4+/MAPK8+ Pre-Ciliated",
                              'HPSE2+ Progenitor' = "HPSE2+ Progenitor")
  
  
  
  Epi_Filter2@active.ident <- factor(x = Epi_Filter2@active.ident, 
                                     levels = 
                                       c("KRT7+ Secretory", 
                                         "PAX8+/LGR5+ Progenitor",
                                         "TEX14+ Progenitor",
                                         "HPSE2+ Progenitor",
                                         "APOA1+ Progenitor",
                                         "KRT5+ Progenitor",
                                         "AGBL4+/MAPK8+ Pre-Ciliated", 
                                         
                                         "CFAP299+ Ciliated",
                                         "TPPP3+ Ciliated"))
  
  
  colors <- c('#35EFEF',
              "#F28D86", #2
              "#2188F7", #3
              "#59D1AF", #5
              '#00A1C6',
              "#EA68E1", #4
              "#B20224", #1
              "#A374B5", #8
              "#9000C6") #9
  
  
  DimPlot(Epi_Filter2, 
          cols = colors,
          pt.size = .9,
          shuffle = T) #9)
  
  DimPlot(Epi_Filter2, 
          group.by = 'Location',
          pt.size = .9,
          shuffle = T) #9)
  
  
  FeaturePlot(Epi_Filter2, features = c("FAM183A"), pt.size = 0.75 , raster = F)+
    scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
    theme(plot.title = element_text(size = 32,face = "bold.italic"))
  
  
  
  
  Epi_Filter2$annotations <- Epi_Filter2@active.ident
  

  ## Stacked Bar Location #
  
  table <- table(Epi_Filter2@meta.data$Location ,
                 Epi_Filter2@active.ident)    # Create a table of counts
  
  IDs = unique(Epi_Filter2@active.ident)
  
  df <- data.frame(table) 
  
  table2 <- table(Epi_Filter2@active.ident)
  
  library(forcats)
  
  y <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
              aes(x = Var1,              # Variable to plot on the x-axis
                  y = Freq,              # Variable to plot on the y-axis
                  fill = factor(Var2,    # Variable to fill the bars
                                levels = c("KRT7+ Secretory", 
                                           "PAX8+/LGR5+ Progenitor",
                                           "TEX14+ Progenitor",
                                           "HPSE2+ Progenitor",
                                           "APOA1+ Progenitor",
                                           "KRT5+ Progenitor",
                                           "AGBL4+/MAPK8+ Pre-Ciliated", 
                                           "CFAP299+ Ciliated",
                                           "TPPP3+ Ciliated")),
              )) + # Order of the stacked bars
    theme_classic() +               # ggplot2 theme
    # Bar plot
    geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
             stat = "identity",     # Height of bars represent values in the data
             size = 1) +            # Size of bars
    # Color scheme
    scale_fill_manual("Location", c("KRT7+ Secretory",
                                    "PAX8+/LGR5+ Progenitor",
                                    "TEX14+ Progenitor",
                                    "HPSE2+ Progenitor",
                                    "APOA1+ Progenitor",
                                    "KRT5+ Progenitor",
                                    "AGBL4+/MAPK8+ Pre-Ciliated", 
                                    "CFAP299+ Ciliated",
                                    "TPPP3+ Ciliated"),
                      values = c("#B20224", #1,
                                 '#35EFEF',
                                 "#F28D86", #2
                                 "#2188F7", #3
                                 "#59D1AF", #5
                                 '#00A1C6',
                                 "#EA68E1", #4
                                 "#A374B5", #8
                                 "#9000C6"))+  
    labs(x = NULL,                     # x-axis label
         y = "Fraction of Cells") +    # y-axis label
    theme(text = element_text(size = 15),                                      # Text size throughout the plot
          axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
          axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis
  
  ggsave(filename = "20251020_epi_pre-cilia_location_stacked_bar.pdf" , plot = y , width = 10, height = 12, dpi =600)
  
  
