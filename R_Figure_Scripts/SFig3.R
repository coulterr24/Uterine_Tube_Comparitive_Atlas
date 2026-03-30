
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


#### Meta-Data Clusters all Human ####



oranges <- c("#FEB24C", "#FD8D3C")
deeper_reds <- c("#A50F15", "#67000D", "#FC4E2A", "#E31A1C")
green <- "#06402b"
blues <- c("#3182BD", "#08519C", "#6BAED6",  "#6A5ACD")
greens <- c("#005824", "#238B45", "#41AE76", "#66C2A4", 
            "#228B22","#98fb98", "#2CFF05" , "#7FFFD4",
            "#CCFF00")


# Combine the colors into one vector
colors <- c(oranges, deeper_reds, green, blues, greens)

#### Supp Figure 3A ####


table <- table(FT_Named@meta.data$Source ,
               FT_Named@active.ident)    # Create a table of counts

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
        "T/NK Cell 1", 
        "T/NK Cell 2",
        "T/NK Cell 3",
        "T/NK Cell 4",
        "T/NK Cell 5",
        "Cycling Immune", 
        "B Cell", 
        "Macrophage", 
        "Mast")

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

ggsave(filename = "202501125_cell_human_source_dist.pdf", plot = source, width = 10, height = 14, dpi = 600)



#### Supp Figure 3B ####


table <- table(FT_Named@meta.data$Source ,
               FT_Named@meta.data$Preservation)    # Create a table of counts

IDs <- c("Surgery", 
         "Cryo", 
         "Cadaver")

df <- data.frame(table) 
df$Var1 <- factor(df$Var1, 
                  levels = c('Nikitin' , 
                             'Dinh' , 
                             'Ulrich',
                             'Lengyel'))

colors <- c('orange' , '#50ad95' , 'grey')

#colors <- c('#0000a2' , '#e9c716' , '#bc272d' , '#50ad9f')

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

ggsave(filename = "202501125_Preservation_human_source_dist.pdf", plot = source, width = 10, height = 14, dpi = 600)



#### Supp Figure 3C ####


table <- table(FT_Named@meta.data$Source ,
               FT_Named@meta.data$age_group)    # Create a table of counts

IDs <- c("Pre-Menopause", 
         "Peri-Menopause", 
         "Menopause")

df <- data.frame(table) 
df$Var1 <- factor(df$Var1, 
                  levels = c('Nikitin' , 
                             'Dinh' , 
                             'Ulrich',
                             'Lengyel'))

colors <- c('#0000CD' , 'purple' , 'red')

#colors <- c('#0000a2' , '#e9c716' , '#bc272d' , '#50ad9f')

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

ggsave(filename = "202501125_source_human_source_dist.pdf", plot = source, width = 10, height = 14, dpi = 600)




#### Supp Figure 3D ####


table <- table(FT_Named@meta.data$Source ,
               FT_Named@meta.data$Location)    # Create a table of counts


IDs <- c('Isthmus' ,
         'Mid-Portion' , 
         'Ampulla' , 
         'Infundibulum', 
         'Fimbria' , 
         'All')



df <- data.frame(table) 
df$Var1 <- factor(df$Var1, 
                  levels = c('Nikitin' , 
                             'Dinh' , 
                             'Ulrich',
                             'Lengyel'))


df <- na.omit(df)


colors <- c("#f4e6ff", "#ddb8ff", "#c78aff", "#b05cff", "#9a2eff", "#8400ff")
#colors <- c('#0000a2' , '#e9c716' , '#bc272d' , '#50ad9f')

#df <- na.omit(df)


#table2 <- table(FT_Named@active.ident)


location <- ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
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

ggsave(filename = "202501125_location_human_source_dist.pdf", plot = location, width = 10, height = 14, dpi = 600)
























