
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

library(tidyr)

#### Mouse Pseudotime Objects and Analysis ####


HQ_Distal_Epi_Filter <- readRDS(file = "20250123_Distal_Epi_Cells_pseudotime.rds", refhook = NULL)



## Psuedotime and Lineage Assignment ##

cellID <- rownames(HQ_Distal_Epi_Filter@reductions$phate@cell.embeddings)
phate_embeddings <- HQ_Distal_Epi_Filter@reductions$phate@cell.embeddings
pseudotime_vals <- HQ_Distal_Epi_Filter@meta.data$Pseudotime

combined_data <- data.frame(cellID, phate_embeddings, pseudotime_vals)

# Calculate the Average PHATE_1 Value for Pseudotime Points = 0 #
avg_phate_1 <- mean(phate_embeddings[pseudotime_vals == 0, 1])

# Pseudotime Values lower than avge PHATE_1 Embedding will be Negative to split lineages
combined_data$Split_Pseudo <- ifelse(phate_embeddings[, 1] < avg_phate_1, -pseudotime_vals, pseudotime_vals)

# Define Lineage #
combined_data$lineage <- ifelse(combined_data$phate_1 < avg_phate_1, "Secretory",
                                ifelse(combined_data$phate_1 > avg_phate_1, "Ciliogenic", "Progenitor"))


HQ_Distal_Epi_Filter$Pseudotime_Adj <- combined_data$Split_Pseudo
HQ_Distal_Epi_Filter$Lineage <- combined_data$lineage

## Set Bins ##

Pseudotime_Lineage <- HQ_Distal_Epi_Filter

bins <- cut_number(Pseudotime_Lineage@meta.data$Pseudotime_Adj , 40) # Evenly distribute bins 

Pseudotime_Lineage@meta.data$Bin <- bins # Metadata for Bins



table <- data.frame(Pseudotime_Lineage@meta.data$Sample,Pseudotime_Lineage$Bin) 

#saveRDS(Pseudotime_Lineage, file = "20250528_Human_PseudotimeLin.rds")



## Set Idents to PSeudoime Bin ##

time_ident <- SetIdent(Pseudotime_Lineage, value = Pseudotime_Lineage@meta.data$Bin)

av.exp <- AverageExpression(time_ident, return.seurat = T)$RNA # Calculate Avg log normalized expression

#### Figure 4A ####

phate_dif <- DimPlot(Pseudotime_Lineage , 
                     reduction = "phate", 
                     cols = c("#B20224", 
                              "#35EFEF", 
                              "#00A1C6", 
                              "#A374B5", 
                              "#9000C6", 
                              "#EA68E1", 
                              "#59D1AF", 
                              "#2188F7", 
                              "#F28D86"),
                     pt.size = 0.7,
                     shuffle = TRUE,
                     seed = 0,
                     label = FALSE)+  
  labs(title = 'Colored by Cluster')+        # Plot title
  NoLegend() +
  labs(x="UMAP_1",y="UMAP_2")

ggsave(filename = "Fig4a_epi_phate.pdf", plot = phate_dif, width = 15, height = 12, dpi = 600)


#### Figure 4B ####

cds <- readRDS(file = "20221101_Distal_Epi_PHATE_Monocle3.rds" , refhook = NULL)

pseudtotime <- plot_cells(cds, 
                          color_cells_by = "pseudotime",
                          label_cell_groups=FALSE,
                          label_leaves=FALSE,
                          label_branch_points=FALSE,
                          graph_label_size=0,
                          cell_size = .01,
                          cell_stroke = 1)+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.line.x = element_blank())+
  theme(axis.line.y = element_blank())+
  theme(axis.ticks.x = element_blank())+
  theme(axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_blank())+
  theme(legend.text = element_text(size = 12))

ggsave(filename = "Fig4b_epi_pseudtotime.pdf", plot = pseudtotime, width = 18, height = 12, dpi = 600)




#### Figure 4C ####

colors <- c(
  "#FF0000", "#FF1A00", "#FF3300", "#FF4D00", "#FF6600", "#FF7F00",
  "#FF9900", "#FFB300", "#FFCC00", "#FFE600", "#FFFF00", "#E6FF00",
  "#CCFF00", "#B3FF00", "#99FF00", "#80FF00", "#66FF00", "#4DFF00",
  "#33FF00", "#1AFF00", "#00FF00", "#00FF1A", "#00FF33", "#00FF4D",
  "#00FF66", "#00FF80", "#00FF99", "#00FFB3", "#00FFCC", "#00FFE6",
  "#00FFFF", "#1A00FF", "#3300FF", "#4D00FF", "#6600FF", "#7F00FF",
  "#9900FF", "#B300FF", "#CC00FF", "#E600FF"
)



rainbow_pseudo <- DimPlot(Pseudotime_Lineage , reduction = "phate", 
                          cols = colors2,
                          pt.size = 0.7,
                          shuffle = TRUE,
                          seed = 0,
                          label = F,
                          group.by = "Bin")+    
  NoLegend()


ggsave(filename = "Fig4d_epi_pseudtotime_bins.pdf", plot = rainbow_pseudo, width = 18, height = 12, dpi = 600)











####


Distal_Epi_Filter <- readRDS(file = "20241216_Distal_Epi_Cells2.rds", refhook = NULL)


SaveRDS(cds_from_seurat, file = "20250123_Distal_Epi_PHATE_monocle.rds")


pseudo <- pseudotime(cds_from_seurat)

HQ_Distal_Epi_Filter@meta.data$Pseudotime <- pseudo # Add to Seurat Metadata

saveRDS(HQ_Distal_Epi_Filter, file = "20250123_Distal_Epi_Cells_pseudotime.rds")

#### Figure 4G-1 ####

features <- c('Spdef' , 'Mif' , 
              'Ifitm2' , 'Selenop' ,
              'Krt7', 'Ovgp11', 'Pax8' , "Slc1a3" , "Itga6",
              'Sox5' , 'Klf6' , 'Krt5' , 'Trp53','Trp73',
              "Prom1", "Foxj1" , "Dnah12" , "Cfap126",
              "Tppp3", "Fam183b")


time_ident <- SetIdent(Pseudotime_Lineage, value = Pseudotime_Lineage@meta.data$Bin)

av.exp <- AverageExpression(time_ident, return.seurat = T)$RNA # Calculate Avg log normalized expression

bin_list <- unique(Pseudotime_Lineage@meta.data$Bin) 

plot_info <- as.data.frame(av.exp[features,]) # Call Avg Expression for features


z_score <- transform(plot_info, SD=apply(plot_info,1, mean, na.rm = TRUE))
z_score <- transform(z_score, MEAN=apply(plot_info,1, sd, na.rm = TRUE))

z_score1 <- (plot_info-z_score$MEAN)/z_score$SD



plot_info$y <- rownames(plot_info) # y values as features
z_score1$y <- rownames(plot_info)


plot_info <- gather(data = plot_info, x, expression, bin_list) #set plot
z_score1 <- gather(data = z_score1, x, z_score, bin_list) #set plot




# calculate max and min z-scores
max_z <- max(z_score1$z_score, na.rm = TRUE)
min_z <- min(z_score1$z_score, na.rm = TRUE)

# set color for outliers
outlier_color <- ifelse(z_score1$z_score > max_z | z_score1$z_score < min_z, ifelse(z_score1$z_score > 0, "#AD1F24", "#51A6DC"), "#e2e2e2")




# Set different na.value options for positive and negative values
na_color_pos <- "#AD1F24" # color for positive NA values
na_color_neg <- "#51A6DC" # color for negative NA values

#custom_bin_names <- c(paste0("S", 20:1), paste0("C", 1:20))


custom_bin_names <- c(paste0("B", 1:40))


figure <- ggplot(z_score1, aes(x, y, fill = z_score)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradientn(colors=c("#1984c5", "#e2e2e2", "#c23728"), 
                       name = "Average Expression \nZ-Score", limits = c(-3,3), 
                       na.value = ifelse(is.na(z_score1) & z_score1 > 0, na_color_pos, 
                                         ifelse(is.na(z_score1) & z_score1 < 0, na_color_neg, "grey50")),
                       oob = scales::squish)+
  scale_x_discrete(limits= sort(bin_list) , labels= custom_bin_names)+
  scale_y_discrete(limits= rev(features))+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 0, hjust = 0.5, size = 10, face = "bold"),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold.italic"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-0.5,1,1,1), "cm"))


ggsave(filename = "Fig4g_mo_pseudotime_binned.pdf" , plot = figure , width = 18, height = 9, dpi =600)



#### Human Pseudotime Objects and Analysis ####


HQ_Distal_Epi_Filter <- readRDS(file = "20250123_Distal_Epi_Cells_pseudotime.rds", refhook = NULL)

#### Figure 4D ####

phate <-  DimPlot(object = Distal_Epi_Named, 
                  reduction = 'phate', 
                  group.by = "Sample.ID", 
                  repel = TRUE, 
                  cols = colors,
                  label = F, 
                  pt.size = 0.9, 
                  label.size = 5,
                  raster = F) + 
  labs(title = 'Colored by Cluster')


ggsave(filename = "Fig4d_epi_human_PHATE.pdf" , plot = phate , width = 15, height = 10, dpi =600)

#### Figure 4E ####

cds_from_seurat <- readRDS(filename = "Fig4e_epi_human_PHATE_Monocle3.pdf" , refhook = NULL)


plot_cells<- plot_cells(cds_from_seurat, 
                        color_cells_by = 'cluster',
                        label_groups_by_cluster=TRUE,
                        label_leaves=FALSE,
                        label_branch_points=TRUE,
                        graph_label_size=4,
                        group_label_size = 0,
                        cell_size = 0.9,)+ 
  scale_color_manual(values = colors)


ggsave(filename = "20250123_epi_human_PHATE_Monocle3.pdf" , plot = plot_cells , width = 15, height = 12, dpi =600)



#### Figure 4F ####

HQ_Distal_Epi_Filter <- readRDS(file = "20250123_Distal_Epi_Cells_pseudotime.rds", refhook = NULL)



## Psuedotime and Lineage Assignment ##

cellID <- rownames(HQ_Distal_Epi_Filter@reductions$phate@cell.embeddings)
phate_embeddings <- HQ_Distal_Epi_Filter@reductions$phate@cell.embeddings
pseudotime_vals <- HQ_Distal_Epi_Filter@meta.data$Pseudotime

combined_data <- data.frame(cellID, phate_embeddings, pseudotime_vals)

# Calculate the Average PHATE_1 Value for Pseudotime Points = 0 #
avg_phate_1 <- mean(phate_embeddings[pseudotime_vals == 0, 1])

# Pseudotime Values lower than avge PHATE_1 Embedding will be Negative to split lineages
combined_data$Split_Pseudo <- ifelse(phate_embeddings[, 1] < avg_phate_1, -pseudotime_vals, pseudotime_vals)

# Define Lineage #
combined_data$lineage <- ifelse(combined_data$phate_1 < avg_phate_1, "Secretory",
                                ifelse(combined_data$phate_1 > avg_phate_1, "Ciliogenic", "Progenitor"))


HQ_Distal_Epi_Filter$Pseudotime_Adj <- combined_data$Split_Pseudo
HQ_Distal_Epi_Filter$Lineage <- combined_data$lineage

## Set Bins ##

Pseudotime_Lineage <- HQ_Distal_Epi_Filter

bins <- cut_number(Pseudotime_Lineage@meta.data$Pseudotime_Adj , 40) # Evenly distribute bins 

Pseudotime_Lineage@meta.data$Bin <- bins # Metadata for Bins



table <- data.frame(Pseudotime_Lineage@meta.data$Sample,Pseudotime_Lineage$Bin) 

write.csv(table, file = "20250604_epi_human_pseudotime_labels_for_SAMap.csv")

#saveRDS(Pseudotime_Lineage, file = "20250528_Human_PseudotimeLin.rds")



## Set Idents to PSeudoime Bin ##

time_ident <- SetIdent(Pseudotime_Lineage, value = Pseudotime_Lineage@meta.data$Bin)

av.exp <- AverageExpression(time_ident, return.seurat = T)$RNA # Calculate Avg log normalized expression

# Calculates Average Expression for Each Bin #
# if you set return.seurat=T, NormalizeData is called which by default performs log-normalization #
# Reported as avg log normalized expression #



colors2 <- c( "#C000FF" , "#A000FF" , "#0020FF" , "#00A0FF", '#00F0FF',"#00FFA0" , 
              "#A0FF00" , "#FFFF00" , "#FFC000" , "#FF8000" , "#FF0000",
              "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000",
              "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000", "#FF0000",
              "#FF0000", '#FF0000',
              '#FF8000','#FFA000','#FFE000','#FFFF00','#C0FF00',
              '#A0FF00','#00FF00','#00FFA0','#00F0FF',
              '#00A0FF','#0020FF','#4000FF','#8000FF',
              '#A000FF','#C000FF')



length(colors2)



rainbow_pseudo <- DimPlot(Pseudotime_Lineage , reduction = "phate", 
                          cols = colors2,
                          pt.size = 0.7,
                          shuffle = TRUE,
                          seed = 0,
                          label = F,
                          group.by = "Bin")+    
  NoLegend()




ggsave(filename = "Fig4f_epi_pseudotime_bin_plot_colors2.pdf" , plot = rainbow_pseudo , width = 15, height = 12, dpi =600)


#### Figure 4G-2 ####


## Genes ##


features <- c('SPDEF' , 'MIF' , 
              'IFITM2' , 'SELENOP' ,
              'KRT7', 'OVGP1', 'PAX8' , "SLC1A3" , "ITGA6",
              'SOX5' , 'KLF6' , 'KRT5' , 'TP53','TP73',
              "PROM1", "FOXJ1" , "DNAH12" , "CFAP126",
              "TPPP3", "FAM183A")

time_ident <- SetIdent(Pseudotime_Lineage, value = Pseudotime_Lineage@meta.data$Bin)

av.exp <- AverageExpression(time_ident, return.seurat = T)$RNA # Calculate Avg log normalized expression

bin_list <- unique(Pseudotime_Lineage@meta.data$Bin) 

plot_info <- as.data.frame(av.exp[features,]) # Call Avg Expression for features


z_score <- transform(plot_info, SD=apply(plot_info,1, mean, na.rm = TRUE))
z_score <- transform(z_score, MEAN=apply(plot_info,1, sd, na.rm = TRUE))

z_score1 <- (plot_info-z_score$MEAN)/z_score$SD



plot_info$y <- rownames(plot_info) # y values as features
z_score1$y <- rownames(plot_info)


plot_info <- gather(data = plot_info, x, expression, bin_list) #set plot
z_score1 <- gather(data = z_score1, x, z_score, bin_list) #set plot




# calculate max and min z-scores
max_z <- max(z_score1$z_score, na.rm = TRUE)
min_z <- min(z_score1$z_score, na.rm = TRUE)

# set color for outliers
outlier_color <- ifelse(z_score1$z_score > max_z | z_score1$z_score < min_z, ifelse(z_score1$z_score > 0, "#AD1F24", "#51A6DC"), "#e2e2e2")




# Set different na.value options for positive and negative values
na_color_pos <- "#AD1F24" # color for positive NA values
na_color_neg <- "#51A6DC" # color for negative NA values

#custom_bin_names <- c(paste0("S", 20:1), paste0("C", 1:20))


custom_bin_names <- c(paste0("B", 1:40))


figure <- ggplot(z_score1, aes(x, y, fill = z_score)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradientn(colors=c("#1984c5", "#e2e2e2", "#c23728"), 
                       name = "Average Expression \nZ-Score", limits = c(-3,3), 
                       na.value = ifelse(is.na(z_score1) & z_score1 > 0, na_color_pos, 
                                         ifelse(is.na(z_score1) & z_score1 < 0, na_color_neg, "grey50")),
                       oob = scales::squish)+
  scale_x_discrete(limits= sort(bin_list) , labels= custom_bin_names)+
  scale_y_discrete(limits= rev(features))+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 0, hjust = 0.5, size = 10, face = "bold"),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold.italic"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-0.5,1,1,1), "cm"))


ggsave(filename = "Fig4g_hu_pseudotime_binned.pdf" , plot = figure , width = 18, height = 9, dpi =600)




list <- 1:40
colors = c(colors2)
df <- data.frame(data = list, color = colors)

binning <- ggplot(df, aes(x = 1:40, y = 1, fill = color)) + 
  geom_bar(stat = "identity",position = "fill", color = "black", size = 1, width = 1) +
  scale_fill_identity() +
  theme_void()+ 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Pseudotime Bin ")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust =1, vjust = .75, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))

library(egg)

psuedotime_lineage <- ggarrange(`binning`,
                                figure,
                                ncol=1,
                                heights = c(20 , (20*length(features)),
                                            widths = c(1)),
                                padding = unit(0.01))



ggsave(filename = "Fig4g_hu_pseudotime_binned_and_annotated.pdf", 
       plot = psuedotime_lineage, width = 16, height = 20, dpi = 600)


