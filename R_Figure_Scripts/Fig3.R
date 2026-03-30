
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







hu_genes <- c('GAS5','EEF1B2','EEF2', 'MIF', 'KRT7', 'S100A16',  #Secretory_Corr
           'UPK1A','ECRG4','GSTM1', 'GSTM2',  #Secretory_Adj
           'PAX8','ITGA6', #Literature Stem
           'ARHGAP26','ESR1','ARIH1','LAMC2', #Stem Corr
           'CRISP3','NRG3','GPC6', 'ROBO1', 'PLA2R1', #Stem Ad Hu
           'SLC1A3', 'TRPM3','MIR100HG','FN1','NCAM1','FAT3', #Stem Ad Mo
           'TEX14','NUP210L','MAFB', #hu Prog adj1
           'HPSE2','UNC5D','LSAMP', #hu Prog adj2
           'APOA1','PRSS33','PKHD1', #hu Prog adj3
           'KRT5','S100A2','SFN', #hu Prog adj4
           'KLF6','FOSB','CDH1', #Prog Corr Mo
           'RHOJ','IFRD1','HBEGF','KALRN', #Prog Adj Mo
           'AGBL4','NEK10','TTC29','KIF6','KNDC1', #PreCil Corr
           'NELL2','CCBE1','GDF15', #PreCil Adj Hu
           'KCTD8','MMP8','SNHG11','AREG', #PreCil Adj Mo
           'FOXJ1','TP73', #Cilia 2 lit
           'CFAP299','CDHR3','DNAH12', #Cilia 1 Corr
           'FAM183A','TPPP3','PIFO','CAPSL' #Cilia 2 Corr
           
           
           
)


mo_genes <- c('Gas5','Eef1b2','Eef2', 'Mif', 'Krt7', 'S100a16',  #Secretory_Corr
           'Upk1a','Ecrg4','Gstm1', 'Gstm2',  #Secretory_Adj
           'Pax8','Itga6', #Literature Stem
           'Arhgap26','Esr1','Arih1','Lamc2', #Stem Corr'
           'Crisp3','Nrg3','Gpc6', 'Robo1', 'Pla2r1',
           'Slc1a3','Trpm3','Mir100hg','Fn1','Ncam1','Fat3', #Stem Adj Mo
           'Tex14','Nup210l','Mafb', #hu Prog adj1
           'Hpse2','Unc5d','Lsamp', #hu Prog adj2
           'Apoa1','Prss33','Pkhd1', #hu Prog adj3
           'Krt5','S100a2','Sfn', #hu Prog adj4
           'Klf6','Fosb','Cdh1', #Prog Corr Mo
           'Rhoj','Ifrd1','Hbegf','Kalrn', #Prog Adj Mo
           'Agbl4','Nek10','Ttc29','Kif6','Kndc1', #PreCil Corr
           'Nell2','Ccbe1','Gdf15',#PreCil Adj Hu
           'Kctd8','Mmp8','Snhg11','Areg', #PreCil Adj Mo
           'Foxj1','Trp73', #Lit Cilia
           'Cfap299','Cdhr3','Dnah2', #Cilia 1 Corr
           'Fam183b','Tppp3','Pifo','Capsl' #Cilia 2 Corr
           
)           
           

## Mouse DP ##
genes <- c( 'Fshr', 'Rhox8', 'Nr5a2', 'Cyp19a1','Hoxb9', 'Peg3', 'Ccn1', 'Wnt6')

all_dp <- DotPlot(object = Mo_Epi_Named,                    # Seurat object
                  assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                  features = mo_genes,                  # List of features (select one from above or create a new one)
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
  scale_y_discrete(limits = levels(Mo_Epi_Named))+
  theme(legend.title = element_text(size = 14))+
  coord_flip()

ggsave(filename = "20260311_Correlation_Dot_Plot_mouse.pdf", plot = all_dp, width = 8, height = 22, dpi = 600)








all_dp <- DotPlot(object = Hu_Epi_Named,                    # Seurat object
                  assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                  features = hu_genes,                 # List of features (select one from above or create a new one)
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
  scale_y_discrete(limits = levels(Hu_Epi_Named))+
  theme(legend.title = element_text(size = 14))+
  coord_flip()


ggsave(filename = "20251223_Correlation_Dot_Plot_human.pdf", plot = all_dp, width = 8, height = 22, dpi = 600)






