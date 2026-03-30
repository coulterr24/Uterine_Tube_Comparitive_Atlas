
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



#### Load Objects #### 

LQ_Cells <- readRDS(file = "20241211_Healthy_LQ_Cells.rds", refhook = NULL)




LQ_Cells$Sample.ID <- factor(x = LQ_Cells$Sample.ID, 
                             levels = 
                               c('D33',
                                 'P33',
                                 'D34',
                                 'P34',
                                 'D35',
                                 'P35',
                                 'D36',
                                 'P36',
                                 'D37',
                                 'P37',
                                 'D38',
                                 'P38',
                                 'p1',
                                 'p2',
                                 'p3',
                                 'p4',
                                 'p5l',
                                 'p5r',
                                 'p6f',
                                 'p6a',
                                 'p6i',
                                 'p7',
                                 'p7f',
                                 'p8',
                                 'FT1f',
                                 'FT1a',
                                 'FT1i',
                                 'FT2f',
                                 'FT2a',
                                 'FT2i',
                                 'FT3f',
                                 'FT3a',
                                 'FT3i',
                                 'FT4',
                                 'D3I',
                                 'D3A',
                                 'D3F',
                                 'D4I',
                                 'D4F',
                                 'D5I',
                                 'D5A',
                                 'D5F',
                                 'D6I',
                                 'D6A',
                                 'D6F',
                                 'D7A',
                                 'D7F'
                               ))






HQ_Cells <- readRDS(file = "20241212_FT_Healthy_Cells_Final.rds", refhook = NULL)


FT_Named <- RenameIdents(HQ_Cells, 
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








#### Supp Figure 2A ####

purples <- c(
  "#FCFBFD",
  "#EFEDF5",
  "#DADAEB",
  "#BCBDDC",
  "#9E9AC8",
  "#807DBA",
  "#6A51A3",
  "#5B3F99",
  "#4B2C8A",
  "#3B1B7A",
  "#2B0A6A",
  "#1F004B"
)

oranges <- c(
  "#FFF5EB",
  "#FEE6CE",
  "#FDD0A2",
  "#FDAE6B",
  "#FD8D3C",
  "#F16913",
  "#D94801",
  "#C44100",
  "#AD3A00",
  "#963300",
  "#7F2C00",
  "#662400"
)

greens <- c(
  "#F7FCF5",
  "#E5F5E0",
  "#C7E9C0",
  "#A1D99B",
  "#74C476",
  "#41AB5D",
  "#238B45",
  "#006D2C",
  "#004F1F",
  "#003814"
)

yellows <- c(
  "#FFF6CC",  # very light pastel yellow
  "#FFEDB3",
  "#FFE39A",
  "#FFD980",
  "#FCD345",  # reference mid-yellow
  "#F2C83E",
  "#E8BD37",
  "#DEB230",
  "#D4A729",
  "#CA9C22",
  "#C0911B",
  "#B68614",
  "#AF8902"   # deep golden endpoint
)

colors <- c(purples,oranges,greens,yellows)

unfiltered_MT <- VlnPlot(LQ_Cells, features = c("percent.mt"), group.by = 'Sample.ID', pt.size = 0,
                         cols = colors)+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+     # Change X-Axis Text Size
  theme(axis.text.y = element_text(size = 16))+     # Change Y-Axis Text Size
  theme(axis.title.y = element_text(size = 18))+    # Change Y-Axis Title Text Size
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank())

ggsave(filename = "FIGs2a_unfiltered_MT.pdf", plot = unfiltered_MT, width = 18, height = 9, dpi = 600)


#### Supp Figure 2B ####

unfiltered_nFeature <- VlnPlot(LQ_Cells, features = c("nFeature_RNA"), group.by = 'Sample.ID', pt.size = 0,
                               cols = colors)+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 18))+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank()) # Change object to visualize other samples

ggsave(filename = "FIGs2b_unfiltered_nFeature.pdf", plot = unfiltered_nFeature, width = 18, height = 9, dpi = 600)


#### Supp Figure  2C #### 

unfiltered_nCount <- VlnPlot(LQ_Cells, features = c("nCount_RNA"), group.by = 'Sample.ID', pt.size = 0,
                             cols = colors)+
  theme(legend.position = 'none')+
  theme(axis.text.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 16))+
  theme(axis.title.y = element_text(size = 18))+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))+
  theme(axis.title.x = element_blank()) # Change object to visualize other samples

ggsave(filename = "FIGs2c_unfiltered_nCount.pdf", plot = unfiltered_nCount, width = 12, height = 9, dpi = 600)


#### Supp Figure 2D #### 



table <- table(FT_Named@active.ident ,
               FT_Named@meta.data$Classify_Doublets)    # Create a table of counts

df <- data.frame(table) 


doublet <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
                  aes(x = Var1,              # Variable to plot on the x-axis
                      y = Freq,              # Variable to plot on the y-axis
                      fill = factor(Var2,    # Variable to fill the bars
                                    levels = c("Doublet","Singlet")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Doublet", limits = c("Doublet","Singlet"),
                    values = c('#8B0000','#808080')) +
  
  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))+
  coord_flip()



ggsave(filename = "FIGs2d2_doublet_quant.pdf", plot = doublet, width = 10, height = 16, dpi = 600)


#### Supp Figure 2E ####


table <- table(Hu_Epi_Named@active.ident ,
               Hu_Epi_Named@meta.data$Classify_Doublets)    # Create a table of counts

df <- data.frame(table) 




epi_doublet <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
                      aes(x = Var1,              # Variable to plot on the x-axis
                          y = Freq,              # Variable to plot on the y-axis
                          fill = factor(Var2,    # Variable to fill the bars
                                        levels = c("Doublet","Singlet")))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Doublet", limits = c("Doublet","Singlet"),
                    values = c('#8B0000','#808080')) +
  
  # Add plot labels
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))+       # Text color and horizontal adjustment on y-axis
  
  coord_flip()



ggsave(filename = "FIGs2e2_epi_doublet_quant.pdf", plot = epi_doublet, width = 10, height = 16, dpi = 600)


