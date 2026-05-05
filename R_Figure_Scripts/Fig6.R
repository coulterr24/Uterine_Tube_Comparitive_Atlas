
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
library (UCell)

#### Figure 6A ####


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


FT_sub <- FT_Named
Idents(FT_sub) <- "Location"
FT_sub <- subset(FT_sub, idents = c("Isthmus","Ampulla","Fimbria"))

Idents(FT_sub) <- "seurat_clusters"
FT_sub <- RenameIdents(FT_sub, 
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
FT_sub@active.ident <- factor(x = FT_sub@active.ident, 
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

FT_sub <- subset(FT_sub, idents = c("Fibroblast 1",
                                    "Fibroblast 2"))



FT_sub$Location <- factor(FT_sub$Location, levels = c('Isthmus', 'Ampulla', 'Fimbria'))

colors <- c('#C71585' , '#8A2BE2' , '#4B00FF')

Fibro_POSTN <- VlnPlot(FT_sub, features = "POSTN", 
                       group.by = 'Location', 
                       pt.size = 0,
                       cols = colors)


ggsave(filename = "Fig6A_human_fibro_POSTN_vln.pdf" , plot = Fibro_POSTN , width = 15, height = 10, dpi =600)




#### Figure 6C ####

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



UT_sub <- subset(UT_Named, idents = c("Fibroblast 1",
                                      "Fibroblast 2",
                                      "Fibroblast 3"))

Idents(UT_sub) <- "Location"

UT_sub <- subset(UT_sub, idents = c('Proximal', 'Distal'))


UT_sub$Location <- factor(UT_sub$Location, levels = c('Proximal', 'Distal'))


UT_sub_POSTN <- VlnPlot(UT_sub, features = "Postn", 
                        group.by = 'Location', 
                        pt.size = 0)



ggsave(filename = "Fig6C_mouse_fibro_POSTN_vln.pdf" , plot = UT_sub_POSTN , width = 15, height = 10, dpi =600)


#### Figure 6F ####


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


features <- c('YAP1', 'WWTR1',
              'TRPV4', 
              'PIEZO1' , 'PIEZO2')

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



### UCell Features ###



Gene_List <- read.csv("UCell_Gene_Lists.csv")

all_gene_list <- as.list(Gene_List)

all_gene_list <- lapply(all_gene_list, function(x) x[x != ""])

Pseudotime_Lineage <- AddModuleScore_UCell(Pseudotime_Lineage, features = all_gene_list , maxRank=6000)

gene_signature.names <- paste0(names(all_gene_list), "_UCell")

VlnPlot(Pseudotime_Lineage, features = gene_signature.names[17], 
        group.by = , 
        pt.size = 0)


## Genes ##
for (i in 1:length(gene_signature.names)){
  Pseudotime_Lineage@meta.data[gene_signature.names[i]] <- Pseudotime_Lineage@meta.data[gene_signature.names[i]]
}

gene_signature.names[1]
# Create Bin List and expression of features #

bin_list <- unique(Pseudotime_Lineage@meta.data$Bin) 

plot_info



plot_info <- Pseudotime_Lineage@meta.data %>%
  select(Bin, all_of(gene_signature.names)) %>%
  pivot_longer(cols = all_of(gene_signature.names), names_to = "Signature", values_to = "Score")

plot_info$Signature <- factor(plot_info$Signature, levels = gene_signature.names)





# Step 1: Compute Bin means
bin_means <- plot_info %>%
  group_by(Bin, Signature) %>%
  summarise(mean_Score = mean(Score, na.rm = TRUE), .groups = "drop")

# Step 2: Compute z-scores of the Bin means
bin_means_z <- bin_means %>%
  group_by(Signature) %>%
  mutate(z_Score = (mean_Score - mean(mean_Score, na.rm = TRUE)) / sd(mean_Score, na.rm = TRUE)) %>%
  ungroup()



# calculate max and min z-scores
max_z <- max(bin_means_z$z_Score, na.rm = TRUE)
min_z <- min(bin_means_z$z_Score, na.rm = TRUE)

# set color for outliers
outlier_color <- ifelse(bin_means_z$z_Score > max_z | bin_means_z$z_Score< min_z, ifelse(bin_means_z$z_Score > 0, "#AD1F24", "#51A6DC"), "#e2e2e2")




# Set different na.value options for positive and negative values
na_color_pos <- "cyan" # color for positive NA values
na_color_neg <- "magenta" # color for negative NA values

#custom_bin_names <- c(paste0("S", 20:1), paste0("C", 1:20))


custom_bin_names <- c(paste0("B", 1:40))


figure2 <- ggplot(bin_means_z, aes(Bin, Signature, fill = z_Score)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradientn(colors=c("cyan", "#e2e2e2", "magenta"), 
                       name = "Average Expression \nZ-Score", limits = c(-3,3), 
  )+
  scale_x_discrete(limits= sort(bin_list) , labels= custom_bin_names)+
  scale_y_discrete(limits= gene_signature.names)+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 0, hjust = 0.5, size = 10, face = "bold"),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold.italic"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-0.5,1,1,1), "cm"))



library(egg)

psuedotime_lineage <- ggarrange(figure,
                                figure2 , 
                                ncol=1,
                                heights = c((20*length(features)),(20*length(gene_signature.names)),
                                            widths = c(1)),
                                padding = unit(0.01))



ggsave(filename = "Fig6F_pseudotime_maechanotransduction_and_wounding.pdf", 
       plot = psuedotime_lineage, width = 16, height = 20, dpi = 600)




