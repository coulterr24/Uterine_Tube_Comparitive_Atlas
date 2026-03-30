
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




Distal_Epi_Filter <- readRDS(file = "20241216_Distal_Epi_Cells.rds", refhook = NULL)


Distal_Epi_Named <- RenameIdents(Distal_Epi_Filter, 
                         '0' = "LGR5+ Secretory", 
                         '1' = "KRT17+ Secretory", 
                         '2' = "MAPK10+ Ciliated",
                         '3' = "TPPP3+ Ciliated", 
                         '4' = "Fibroblast-like", 
                         '5' = "KRT5+ Secretory", 
                         '6' = "RUNX3+ Basal",
                         '7' = "LMNTD1+ Ciliated",
                         '8' = "Endothelial-like")


Distal_Epi_Named@active.ident <- factor(x = Distal_Epi_Named@active.ident, 
                                levels = 
                                  c("LGR5+ Secretory", 
                                    "KRT17+ Secretory", 
                                    "KRT5+ Secretory",
                                    "LMNTD1+ Ciliated",
                                    "MAPK10+ Ciliated",
                                    "TPPP3+ Ciliated", 
                                    "RUNX3+ Basal",
                                    "Fibroblast-like", 
                                    "Endothelial-like"))




Distal_Epi_Filter <- readRDS(file = "20250123_Distal_Epi_Cells.rds", 
                             refhook = NULL)


Distal_Epi_Named <- RenameIdents(Distal_Epi_Filter, 
                                 '0' = "PAX8+/LGR5+ Progenitor", 
                                 '1' = "KRT7+ Secretory", 
                                 '2' = "CFAP299+ Ciliated",
                                 '3' = "TPPP3+ Ciliated", 
                                 '4' = "KRT5+ Progenitor", 
                                 '5' = "APOA1+ Progenitor", 
                                 '6' = "TEX14+ Progenitor",
                                 '7' = "AGBL4+/MAPK8+ Pre-ciliated",
                                 '8' = "HPSE2+ Progenitor")


Distal_Epi_Named@active.ident <- factor(x = Distal_Epi_Named@active.ident, 
                                        levels = 
                                          c("PAX8+/LGR5+ Progenitor",
                                            "TEX14+ Progenitor",
                                            "HPSE2+ Progenitor",
                                            
                                            "APOA1+ Progenitor",
                                            "KRT5+ Progenitor",
                                            "AGBL4+/MAPK8+ Pre-ciliated", 
                                            "KRT7+ Secretory", 
                                            "CFAP299+ Ciliated",
                                            "TPPP3+ Ciliated"))



## PHATE Embedding #####



assay = 'RNA'
data <- t(GetAssayData(object = Distal_Epi_Filter, slot = 'data', assay = assay))
phate_output <- phateR::phate(data)

# save PHATE result to csv
write.csv(as.data.frame(phate_output), "20241216_distal_epi_phate_embeddings")

# save PHATE result to Seurat
phate_output <- as.matrix(phate_output)
colnames(phate_output) <- paste0("PHATE_", 1:ncol(phate_output))
rownames(phate_output) <- colnames(Distal_Epi_Filter)
phate.reduction <- CreateDimReducObject(
  embeddings = phate_output,
  key = "PHATE_",
  assay = assay
)
Distal_Epi_Filter[["phate"]] <- phate.reduction




### Only this below ###

## Retry and address batch corrections ##

# https://github.com/KrishnaswamyLab/phateR/issues/26

library(batchelor)

sce <- as.SingleCellExperiment(Distal_Epi_Named, assay = NULL)


set.seed(246)

phate.out <- fastMNN(sce, 
                     batch = c(sce$Sample.ID)
                     )

tmp.phate <- phate(
  t(assay(phate.out , "reconstructed")),
  npca = 41,
  knn = 5,
  decay = 100,
  t = 150,
  gamma = 1,
  n.jobs = 6,
  #sample_idx = "Sample.ID",
  #kernel_symm='mnn'
) 

colnames(tmp.phate$embedding) <- paste0("phate_", 1:2)
tmp.phate$params$data <- NULL

Distal_Epi_Filter[["phate"]] <- CreateDimReducObject(
  embeddings = tmp.phate$embedding,
  key = "phate_", 
  assay = 'RNA',
  misc = tmp.phate$params
)

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


ggsave(filename = "20250123_epi_human_PHATE.pdf" , plot = phate , width = 15, height = 10, dpi =600)



DimPlot(object = Distal_Epi_Filter, 
        reduction = 'phate', 
        group.by = "Sample.ID", 
        repel = TRUE, 
        #cols = magma(n = 5, direction = 1, begin=.7,end=0.05),
        label = T, 
        pt.size = 0.3, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cancer Stage')


DimPlot(object = Distal_Epi_Filter3, 
        reduction = 'phate', 
        group.by = "Age", 
        repel = TRUE, 
        #cols = magma(n = 5, direction = 1, begin=.7,end=0.05),
        label = T, 
        pt.size = 0.1, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cancer Stage')


DimPlot(object = Distal_Epi_Filter5, 
        reduction = 'phate', 
        group.by = "Classify_Doublets", 
        repel = TRUE, 
        #cols = magma(n = 5, direction = 1, begin=.7,end=0.05),
        label = T, 
        pt.size = 0.1, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cancer Stage')


Beeg_PHATE <- Distal_Epi_Named

Beeg_PHATE@reductions[["phate"]]@cell.embeddings <- Beeg_PHATE@reductions[["phate"]]@cell.embeddings*100






Feature_PHATE <- FeaturePlot(Beeg_PHATE, reduction = "phate", features = c("NR4A3"), 
            pt.size = 0.9, order = T)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

#ggsave(filename = "FigX_NR4A3_PHATE.pdf", plot = Feature_PHATE, width = 12, height = 9, dpi = 600)


saveRDS(Beeg_PHATE, file = "20250123_HQ_Distal_Epi_Cells.rds")



HQ_Distal_Epi_Filter <- readRDS(file = "20250123_HQ_Distal_Epi_Cells.rds", refhook = NULL)




FeaturePlot(HQ_Distal_Epi_Filter, reduction = "phate", features = c("MCIDAS"), pt.size = 0.4)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))



DimPlot(object = HQ_Distal_Epi_Filter, 
        reduction = 'phate', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        #cols = magma(n = 5, direction = 1, begin=.7,end=0.05),
        label = T, 
        pt.size = 0.1, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cluster')



#### Monocle3 ####

library(tidyr)


gene_annotation <- as.data.frame(rownames(HQ_Distal_Epi_Filter@reductions[["harmony"]]@feature.loadings),
                                 row.names = rownames(HQ_Distal_Epi_Filter@reductions[["harmony"]]@feature.loadings))

colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(HQ_Distal_Epi_Filter@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = HQ_Distal_Epi_Filter@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "sample"


cell_metadata <- separate(cell_metadata, col = 'sample', 
                          into = c('Patient', 'Barcode'), sep = '_')

New_matrix <- HQ_Distal_Epi_Filter@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(HQ_Distal_Epi_Filter@reductions[["harmony"]]@feature.loadings), ]



expression_matrix <- New_matrix


cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


list_cluster <- HQ_Distal_Epi_Filter@active.ident
names(list_cluster) <- HQ_Distal_Epi_Filter@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster



cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- (HQ_Distal_Epi_Filter@reductions[["phate"]]@cell.embeddings*10)

cds_from_seurat

# Process cds file #

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)


colors <- c("#B20224", #1
            "#F28D86", #2
            "#2188F7", #3
            '#00A1C6',
            "#59D1AF", #5
            '#35EFEF',
            "#EA68E1", #4
            "#A374B5", #8
            "#9000C6") #9


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



cds_from_seurat <- order_cells(cds_from_seurat)


pseudotime<- plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.9)


ggsave(filename = "20250123_epi_human_PHATE_Monocle3_pseudotime.pdf" , plot = pseudotime , width = 15, height = 12, dpi =600)


saveRDS(cds_from_seurat, file = "20250123_Distal_Epi_PHATE_monocle.rds")


pseudo <- pseudotime(cds_from_seurat)

HQ_Distal_Epi_Filter@meta.data$Pseudotime <- pseudo # Add to Seurat Metadata

saveRDS(HQ_Distal_Epi_Filter, file = "20250123_Distal_Epi_Cells_pseudotime.rds")




