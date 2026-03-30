
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

Mo_Epi <- readRDS(file = "20220817_Mouse_Distal_Epi_Cells.rds" , refhook =  NULL)

Mo_Epi_Named <- RenameIdents(Mo_Epi, 
                             '0' = "Spdef+ Secretory", 
                             '1' = "Slc1a3+ Stem/Progenitor", 
                             '2' = "Cebpdhigh/Foxj1- Progenitor",
                             '3' = "Ciliated 1", 
                             '4' = "Ciliated 2", 
                             '5' = "Pax8low/Prom1+ Cilia-forming", 
                             '6' = "Fibroblast-like",
                             '7' = "Slc1a3med/Sox9+ Cilia-forming",
                             '8' = "Selenop+/Gstm2high Secretory")

Mo_Epi_Named@active.ident <- factor(x = Mo_Epi_Named@active.ident, levels = c( c("Spdef+ Secretory",
                                                                                 "Selenop+/Gstm2high Secretory",
                                                                                 "Slc1a3+ Stem/Progenitor",
                                                                                 "Cebpdhigh/Foxj1- Progenitor",
                                                                                 "Slc1a3med/Sox9+ Cilia-forming",
                                                                                 "Pax8low/Prom1+ Cilia-forming", 
                                                                                 "Fibroblast-like",
                                                                                 "Ciliated 1",
                                                                                 "Ciliated 2")))




Mo_Epi_Named <- subset(Mo_Epi_Named, 
                       idents = c("Spdef+ Secretory",
                                  "Selenop+/Gstm2high Secretory",
                                  "Slc1a3+ Stem/Progenitor",
                                  "Cebpdhigh/Foxj1- Progenitor",
                                  "Slc1a3med/Sox9+ Cilia-forming",
                                  "Pax8low/Prom1+ Cilia-forming", 
                                  "Ciliated 1",
                                  "Ciliated 2"))









Hu_Epi <- readRDS(file = "20250123_Distal_Epi_Cells.rds" , refhook =  NULL)


Hu_Epi_Named <- RenameIdents(Hu_Epi, 
                             '0' = "PAX8+/LGR5+ Progenitor", 
                             '1' = "KRT7+ Secretory", 
                             '2' = "CFAP299+ Ciliated",
                             '3' = "TPPP3+ Ciliated", 
                             '4' = "KRT5+ Progenitor", 
                             '5' = "APOA1+ Progenitor", 
                             '6' = "TEX14+ Progenitor",
                             '7' = "AGBL4+/MAPK8+ Pre-ciliated",
                             '8' = "HPSE2+ Progenitor")



Hu_Epi_Named@active.ident <- factor(x = Hu_Epi_Named@active.ident, 
                                    levels = 
                                      c("KRT7+ Secretory", 
                                        "PAX8+/LGR5+ Progenitor",
                                        "TEX14+ Progenitor",
                                        "HPSE2+ Progenitor",
                                        "APOA1+ Progenitor",
                                        "KRT5+ Progenitor",
                                        "AGBL4+/MAPK8+ Pre-ciliated", 
                                        "CFAP299+ Ciliated",
                                        "TPPP3+ Ciliated"))



#### Figure 2A ####

epi_umap <- DimPlot(object = Mo_Epi_Named,                # Seurat object  
                    reduction = 'umap',                 # Axes for the plot (UMAP, PCA, etc.) 
                    repel = TRUE,                       # Whether to repel the cluster labels
                    label = FALSE,                       # Whether to have cluster labels 
                    cols = c("#35EFEF", #1
                             "#00A1C6", #2
                             "#2188F7", #3
                             "#EA68E1", #4
                             "#59D1AF", #5
                             "#B20224", #6
                             "#F28D86", #7
                             "#A374B5", #8
                             "#9000C6"), #9
                    
                    pt.size = 0.6,                      # Size of each dot is (0.1 is the smallest)
                    label.size = 0.5) +                   # Font size for labels    
  # You can add any ggplot2 1customizations here
  labs(title = 'Colored by Cluster')+        # Plot title
  NoLegend() +
  labs(x="UMAP_1",y="UMAP_2")

ggsave(filename = "Fig2a_epi_umap.pdf", plot = epi_umap, width = 15, height = 12, dpi = 600)


#### Figure 2B ####



colors <- c('#35EFEF',
            "#F28D86", #2
            "#2188F7", #3
            "#59D1AF", #5
            '#00A1C6',
            "#EA68E1", #4
            "#B20224", #1
            "#A374B5", #8
            "#9000C6") #9




epi_human <- DimPlot(object = Hu_Epi_Named, 
                     reduction = 'umap', 
                     #group.by = "seurat_clusters", 
                     repel = TRUE, 
                     cols = colors,
                     label = F, 
                     pt.size = 0.6, 
                     label.size = 5,
                     raster = F) + 
  labs(title = 'Colored by Cluster')


ggsave(filename = "20250205_epi_human_UMAP.pdf" , plot = epi_human , width = 15, height = 12, dpi =600)

