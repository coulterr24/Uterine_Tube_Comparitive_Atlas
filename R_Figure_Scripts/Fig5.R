
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


#### Figure 5B ####



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



Idents(Pseudotime_Lineage) <- Pseudotime_Lineage@meta.data$Bin

DimPlot(Pseudotime_Lineage)

Pseudotime_Lineage2 <- RenameIdents(Pseudotime_Lineage, 
                                    '[-40.2979,-35.3956]' = "Secretory" , 
                                    '(-35.3956,-32.525]' =  "Secretory" , 
                                    '(-32.525,-26.4568]' =  "Secretory" , 
                                    '(-26.4568,-22.9389]'  =  "Secretory", 
                                    '(-22.9389,-21.0667]' = 'Secretory',
                                    '(-21.0667,-19.4567]' =    "Secretory" , 
                                    '(-19.4567,-17.9258]' =   "Secretory" , 
                                    '(-17.9258,-16.5403]' =  "Secretory" , 
                                    '(-16.5403,-15.645]' =    "Secretory" , 
                                    '(-15.645,-14.6551]' =   "Secretory" , 
                                    '(-14.6551,-13.5987]' =    "Progenitor",
                                    '(-13.5987,-12.4014]' =  "Progenitor", 
                                    '(-12.4014,-11.1275]' =  "Progenitor", 
                                    '(-11.1275,-10.0766]' =   "Progenitor", 
                                    '(-10.0766,-9.26489]' =    "Progenitor", 
                                    '(-9.26489,-8.49165]' =     "Progenitor", 
                                    '(-8.49165,-7.57875]' =  "Progenitor",
                                    '(-7.57875,-6.67719]' =     "Progenitor", 
                                    '(-6.67719,-5.32179]' =   "Progenitor", 
                                    '(-5.32179,-3.4411]' =   "Progenitor", 
                                    '(-3.4411,-1.87455]' =   "Progenitor", 
                                    '(-1.87455,0.0818802]' =   "Progenitor", 
                                    '(0.0818802,4.49992]' =  "Progenitor",
                                    '(4.49992,6.1506]' =    "Progenitor", 
                                    '(6.1506,13.224]' =   'Progenitor',
                                    '(13.224,21.8303]' =       'Ciliated',
                                    '(21.8303,40.9479]'  =    'Ciliated',
                                    '(40.9479,46.2216]'   = 'Ciliated',
                                    '(46.2216,48.7064]'  =  'Ciliated',
                                    '(48.7064,50.1256]' =     'Ciliated',
                                    '(50.1256,50.1257]'  =   'Ciliated',
                                    '(50.1257,51.6946]'  =  'Ciliated',
                                    '(51.6946,52.5965]'  =   'Ciliated',
                                    '(52.5965,53.891]'  =    'Ciliated',
                                    '(53.891,54.7246]'  =  'Ciliated',
                                    '(54.7246,55.1099]'   =    'Ciliated',
                                    '(55.1099,58.4321]'  =    'Ciliated',
                                    '(58.4321,65.0151]'  =   'Ciliated',
                                    '(65.0151,74.0136]'  =  'Ciliated',
                                    '(74.0136,76.0726]' =    'Ciliated')

sort(unique(Pseudotime_Lineage@active.ident))

DimPlot(Pseudotime_Lineage2)


colors2 <- c( "Secretory" , 
              "Secretory" , 
              "Secretory" , 
              "Secretory", 
              'Secretory',
              "Secretory" , 
              "Secretory" , 
              "Secretory" , 
              "Secretory" , 
              "Secretory" , 
              "Progenitor",
              "Progenitor", 
              "Progenitor", 
              "Progenitor", 
              "Progenitor", 
              "Progenitor", 
              "Progenitor",
              "Progenitor", 
              "Progenitor", 
              "Progenitor", 
              "Progenitor", 
              "Progenitor", 
              "Progenitor",
              "Progenitor", 
              'Progenitor',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated',
              'Ciliated')















gene_list <- c('ROBO1', 'ITGAV', 'NRCAM', 'BMPR2', 'ABCC4',
               'CLIC4', 'ITGB8', 'INSR', 'PTPRJ', 'ST14', 'MYH9',
               'MET', 'NEO1', 'BMPR1A', 'MYH10', 'EGFR', 'ANO6', 
               'ITGA6', 'MSN', 'NLGN1', 'PLA2R1', 'CD55', 'RALA',
               'ITGA2', 'FGFR2', 'LDLR', 'NOTCH2', 'SORL1',
               'PLAT', 'ALCAM', 'EPHA2', 'TMC1', 'WNT5B', 'PDIA3',
               'PDGFC', 'C3', 'PLAUR', 'CPM')



avg_exp <- AverageExpression(Pseudotime_Lineage2, features = gene_list, return.seurat = FALSE)
x <- avg_exp$RNA[gene_list, ]

write.csv(x, file = "20250923_Surface_Marker_Specificity.csv")


#### Figure 5C ####


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


time_ident <- SetIdent(Pseudotime_Lineage, value = Pseudotime_Lineage@meta.data$Bin)

av.exp <- AverageExpression(time_ident, return.seurat = T)$RNA # Calculate Avg log normalized expression

# Calculates Average Expression for Each Bin #
# if you set return.seurat=T, NormalizeData is called which by default performs log-normalization #
# Reported as avg log normalized expression #


features <- c('KRT7', 
              'ROBO1', 'PLA2R1' , 
              "FOXJ1")

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


ggsave(filename = "Fig5c_hu_pseudotime_binned_analgous_stem_markers.pdf" , plot = figure , width = 18, height = 9, dpi =600)











