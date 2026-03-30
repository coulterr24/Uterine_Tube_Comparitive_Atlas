
getwd()
setwd("/local/workdir/cqr3/human_FT_samples/All_Healthy_Human_FT/")
getwd()

#### Load Libraries and Functions ####
library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(ggplot2)
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

#### Metadata ####
meta <- fread("Human_FT_Healthy_meta_2.csv")



# Generate/preprocess objects with Seurat pipeline ####
# SoupX ####
soup.list <- list()
soup.mat.list <- list()

num.cores <- 6
cl <- makeCluster(num.cores)

# Load souplist objects from cell ranger outputs
soup.list <- mclapply( 
  as.list(paste0(meta$data.dir,"/outs")),
  FUN = SoupX::load10X,
  keepDroplets=TRUE,
  mc.cores=num.cores
)

soup.list.est <- mclapply(
  soup.list,
  FUN = function(sc){
    return(tryCatch(autoEstCont(sc), error=function(e) NULL))
  },
  mc.cores=num.cores
)

adj.mat.list <- mclapply(
  soup.list.est,
  FUN = function(sc){
    return(tryCatch(adjustCounts(sc), error=function(e) NULL))
  },
  mc.cores=num.cores
)
stopCluster(cl)


#Save adjust matrices to disk!
for(i in 1:length(adj.mat.list)){
  if(!is.null(adj.mat.list[[i]]) & !file.exists(paste0(meta$data.dir[i],"/outs/soupx/matrix.mtx"))){
    DropletUtils:::write10xCounts(
      path=paste0(meta$data.dir[i],"/outs/soupx"), #path to each sample's cellranger count output
      adj.mat.list[[i]]
    )
  }
}

# save Rho values
rhos <- list()
for(i in 1:length(soup.list.est)){
  rhos[[i]] <- mean(soup.list.est[[i]]$metaData$rho)
}
rhos <- do.call(rbind,rhos)
meta$soupx.rho <- rhos

# Save metadata file with the Rho values
write.table(x = meta, 
            file = "Human_FT_Healthy_meta_2_rho.csv",
            sep = ",",
            row.names = FALSE)

# Read in the count matrices ####
mat.list <- list()
soupx.used <- list()

for(i in 1:length(meta$data.dir)){ 
  if(file.exists(paste0(meta$data.dir[i], '/outs/soupx'))){ #first try soupx
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i],"/outs/soupx"))
    soupx.used[[i]] <- T
  }else if(file.exists(paste0(meta$data.dir[i], '/outs/filtered_feature_bc_matrix'))){ #then try cellranger outputs
    cat("Reading #",i, ": ", meta$data.dir[i], ' \n')
    mat.list[[i]] <- Read10X(data.dir = paste0(meta$data.dir[i], '/outs/filtered_feature_bc_matrix'))
    soupx.used[[i]] <- F
  }else{
    cat("Data not found for # ", i, " (", meta$data.dir[i], ")", "\n")
    soupx.used[[i]] <- NULL
  }
}

meta$soupx.used <- unlist(soupx.used)
cat(sum(unlist(lapply(mat.list, ncol))),"cells (total) loaded...\n")

# Seurat ####
seur.list <- list()

# Initialize Seurat objects
num.cores <- 5 # Can go higher (maybe 15-20 cores) - 5 cores takes ~40min
cl <- makeCluster(num.cores)

seur.list <- mclapply( 
  mat.list,
  FUN = function(mat){
    return(CreateSeuratObject(
      counts = mat, 
      project = 'human_FT'
    ))
  },
  mc.cores = num.cores
)  
stopCluster(cl)

for(i in 1:length(seur.list)){
  cat(' #####################################\n',
      '### Processing dataset number ', i, '###\n',
      '#####################################\n')
  # Add meta data
  for(md in colnames(meta)){
    seur.list[[i]][[md]] <- meta[[md]][i]
  }
  # add %MT
  seur.list[[i]][["percent.mt"]]  <- PercentageFeatureSet(seur.list[[i]], pattern = "MT-") 
  
}

cat((sum(unlist(lapply(mat.list, ncol)))-sum(unlist(lapply(seur.list, ncol)))),"cells (total) removed...\n")

# Preprocess seurat objects
seuPreProcess <- function(seu, assay='RNA', n.pcs=50, res=0.8){
  # NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  pca.name = paste0('pca_', assay)
  pca.key = paste0(pca.name,'_')
  umap.name = paste0('umap_', assay)
  
  seu = NormalizeData(
    seu
  ) %>% FindVariableFeatures(
    assay = assay,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = F
  ) %>% ScaleData(
    assay = assay
  ) %>% RunPCA(
    assay = assay,
    reduction.name = pca.name,
    reduction.key = pca.key,
    verbose = F,
    npcs = n.pcs
  )
  
  #find pcs to use
  tmp.var <- (seu@reductions[[pca.name]]@stdev)^2
  var.cut <- 0.95*sum(tmp.var)
  j=0
  var.sum = 0
  while(var.sum < 0.95*var.cut){
    j = j + 1
    var.sum <- var.sum + tmp.var[j]
  }
  n.pcs.use = j
  
  # FindNeighbors %>% RunUMAP, FindClusters
  seu <- FindNeighbors(
    seu,
    reduction = pca.name,
    dims = 1:n.pcs.use,
    force.recalc = TRUE,
    verbose = FALSE
  ) %>% RunUMAP(
    reduction = pca.name,
    dims = 1:n.pcs.use,
    reduction.name=umap.name
  )
  
  seu@reductions[[umap.name]]@misc$n.pcs.used <- n.pcs.use
  
  seu <- FindClusters(object = seu,resolution = res)
  seu[[paste0('RNA_res.',res)]] <- as.numeric(seu@active.ident)
  
  return(seu)
}

seur.list <- lapply(seur.list, seuPreProcess)

#### Identify Multiplet Rate w/ Linear Model ####

estimateDoubletRate <- function(seur.list, doublet.dist=NULL){
  if(is.null(doublet.dist)){
    doublet.dist <- data.frame(cells.recovered=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
                               multi.rate=     c(0.4,  0.8,  1.6,  2.3,  3.1,  3.9,  4.6,  5.4,  6.1,  6.9,   7.6)
    )
  }
  
  fit <- lm(multi.rate~cells.recovered, data = doublet.dist)
  fit.func <- function(x){
    return(as.numeric(fit$coefficients['cells.recovered']*x + fit$coefficients['(Intercept)']))
  }
  
  ncells.list <- lapply(seur.list, ncol)
  out.list <- lapply(ncells.list, fit.func)
  return(unlist(out.list))
}

meta$Multiplet_Rate <- estimateDoubletRate(seur.list)

#### DoubletFinder on RNA for Eack Individual Dataset ####

bcmvn <- list()
pK <- list()
homotypic.prop <- list()
nExp_poi <- list()
nExp_poi.adj <- list()

## Prepare Doublet Estimation ##

for(i in 1:length(seur.list)){
  cat(' --------------------------------------------\n',
      '--- DoubletFinder for dataset number ', i, '---\n',
      '--------------------------------------------\n')
  
  ## pK Identification (no ground-truth)
  bcmvn[[i]]<- paramSweep_v3(             
    # pN-pK parameter sweeps on a 10,000-cell subset
    seu=seur.list[[i]], 
    PCs = 1:seur.list[[i]]@reductions$umap_RNA@misc$n.pcs.used, 
  ) %>% summarizeSweep(
    # computing the bimodality coefficient across pN and pK parameter space
    GT = FALSE
  ) %>% find.pK() 
  # Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value tested during doubletFinder_ParamSweep
  
  ## Pull out max of bcmvn
  pK[[i]] <- as.numeric(as.character(bcmvn[[i]]$pK[bcmvn[[i]]$BCmetric==max(bcmvn[[i]]$BCmetric)])) # ugly, but functional...
  
  ## Homotypic Doublet Proportion Estimate
  homotypic.prop[[i]] <- modelHomotypic(seur.list[[i]]$seurat_clusters) 
  
  nExp_poi[[i]] <- round(meta$Multiplet_Rate[i]/100*length(colnames(seur.list[[i]])))  # Rate is 0.0242 for 5k "Target for Cell Recovery"
  nExp_poi.adj[[i]] <- round(nExp_poi[[i]]*(1-homotypic.prop[[i]]))
}

## Run DoubletFinder ##

for(i in 1:length(seur.list)){
  seur.list[[i]] <- 
    doubletFinder_v3( # just changed it so the output metadata column name is customizable
      seu=seur.list[[i]], 
      PCs = 1:seur.list[[i]]@reductions$umap_RNA@misc$n.pcs.used,
      pN = 0.25, #default value
      pK= pK[[i]], 
      nExp = nExp_poi.adj[[i]],  
      reuse.pANN = F
    )
}



# Merging and Integration ####
# Merge datasets ####

# Add a prefix to the cell IDs so that cell IDs are NOT merged when the datasets are merged
Sample.ID <- meta$Sample.ID

# Merge the datasets together
All.merged <- merge(seur.list[[1]], y = seur.list[2:47], add.cell.ids = Sample.ID[1:47])

# Save as RDS file
saveRDS(All.merged, file = "20241211_FT_Healthy_All_Cells.rds")

# Load RDS
All.merged <- readRDS(file = "20241211_FT_Healthy_All_Cells.rds", refhook = NULL)

# Seurat workflow + Harmony integration ####

# Need to normalize and run PCA before running Harmony
All.merged <- NormalizeData(object = All.merged, assay = 'RNA')
All.merged <- FindVariableFeatures(object = All.merged, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
All.merged <- ScaleData(object = All.merged, assay = 'RNA')
All.merged <- RunPCA(object = All.merged, assay = 'RNA', features = VariableFeatures(object = All.merged),
                     reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')

# Determine the 'dimensionality' of the dataset
ElbowPlot(All.merged,
          reduction = 'pca_RNA',
          ndims = 50) 

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

max_pcs <- npcs(seu = All.merged, var.toal = 0.95, reduction = 'pca_RNA')

# Run Harmony
gc()

All.merged <- RunHarmony(All.merged, 
                         group.by.vars = "Sample.ID", 
                         reduction = "pca_RNA", 
                         assay = "RNA", 
                         plot_convergence = TRUE)

# Normal Seurat clustering
All.merged <- FindNeighbors(object = All.merged, 
                            reduction = "harmony",
                            k.param = 70, 
                            dims = 1:max_pcs, 
                            graph.name = 'harmony_snn')

All.merged <- FindClusters(object = All.merged, 
                           resolution = 0.7, 
                           graph.name = 'harmony_snn')

# Set Python UMAP via reticulate
umap.method = 'umap-learn'
metric = 'correlation'

All.merged <- RunUMAP(object = All.merged, reduction = "harmony", dims = 1:max_pcs)

# UMAP 
DimPlot(object = All.merged, 
        reduction = 'umap', 
        group.by = "Classify_Doublets", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.0001, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')

FeaturePlot(object = All.merged,features = "PI16")

## Function to Reduce Meta Data Columns ##

df <- as.data.frame(All.merged@meta.data)

Classify_Doublets <- select(df, contains("DF.c"))
Classify_Doublets2 <- data.frame(na.omit(unlist(Classify_Doublets)))

DF_values <- select(df, contains("pANN_0"))
pANN_Values <- data.frame(na.omit(unlist(DF_values)))

df$Classify_Doublets <- Classify_Doublets2[,1]
df$DF_values <- pANN_Values[,1]

df$Classify_Doublets 
df$DF_values


df <- df[,!names(df) %in% colnames(Classify_Doublets)]
df <- df[,!names(df) %in% colnames(DF_values)]

df


All.merged@meta.data <- df

## Save File ##

saveRDS(All.merged, file = "20241211_Healthy_LQ_Cells.rds")

LQ_Cells <- readRDS(file = "20241211_Healthy_LQ_Cells.rds", refhook = NULL)

FeaturePlot(All.merged, features = 'KRT5')+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

FeaturePlot(All.merged, features = 'percent.mt' , raster=FALSE)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold"))

FeaturePlot(All.merged, features = 'nCount_RNA', raster=FALSE)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold"))

FeaturePlot(All.merged, features = 'nFeature_RNA', raster=FALSE)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold"))

DimPlot(object = All.merged, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = T, 
        pt.size = 0.0001, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cluster')


# Remove low_quality cells ####

HQ_Cells <- subset(All.merged, subset = percent.mt < 20 & nFeature_RNA > 200 & nCount_RNA > 750 )

# Save as RDS file
saveRDS(HQ_Cells, file = "20241212_Healthy_FT_HQ_Cells.rds")

# Load RDS
HQ_Cells <- readRDS(file = "20241212_Healthy_FT_HQ_Cells.rds", refhook = NULL)


# Seurat workflow + Harmony integration (without doublets) ####

# Need to normalize and run PCA before running Harmony
HQ_Cells <- NormalizeData(object = HQ_Cells, assay = 'RNA')
HQ_Cells <- FindVariableFeatures(object = HQ_Cells, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
HQ_Cells <- ScaleData(object = HQ_Cells, assay = 'RNA')
HQ_Cells <- RunPCA(object = HQ_Cells, assay = 'RNA', features = VariableFeatures(object = HQ_Cells),
                   reduction.name = 'pca_RNA', reduction.key = 'pca_RNA_')

# Determine the 'dimensionality' of the dataset
ElbowPlot(HQ_Cells,
          reduction = 'pca_RNA',
          ndims = 50) 

# Calculate the number of PCs that contain some proportion (95%) of the variance
max_pcs <- npcs(seu = HQ_Cells, var.toal = 0.95, reduction = 'pca_RNA')

# Run Harmony
gc()
HQ_Cells <- RunHarmony(HQ_Cells,
                       group.by.vars = "Sample.ID",
                       reduction = "pca_RNA", 
                       assay = "RNA", 
                       plot_convergence = TRUE)

# Normal Seurat clustering
HQ_Cells <- FindNeighbors(object = HQ_Cells, 
                          reduction = "harmony", 
                          k.param = 70, 
                          dims = 1:max_pcs, 
                          graph.name = 'harmony_snn')

HQ_Cells <- FindClusters(object = HQ_Cells, resolution = 0.7, graph.name = 'harmony_snn')

# Set Python UMAP via reticulate
umap.method = 'umap-learn'
metric = 'correlation'

HQ_Cells <- RunUMAP(object = HQ_Cells, reduction = "harmony", dims = 1:max_pcs)

# UMAP 
DimPlot(object = HQ_Cells, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.01, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = HQ_Cells, 
        reduction = 'umap', 
        group.by = "Source", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.01, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = HQ_Cells, 
        reduction = 'umap', 
        group.by = "Classify_Doublets", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.01, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')

# Save If Looks Good ##

saveRDS(HQ_Cells, file = "20241212_FT_Healthy_Cells_Final.rds")


# Load RDS
HQ_Cells <- readRDS(file = "20241212_FT_Healthy_Cells_Final.rds", refhook = NULL)


## Identify Cell Types ##

FTAllMarkers <- FindAllMarkers(HQ_Cells, 
                               max.cells.per.ident = 5000)

write.csv(FTAllMarkers, file = "20241212_FT_All_Markers.csv")


features <- c("VIM","COL1A1","GATA2",                               # Fibroblast
              "DES","PALLD","ACTA2","MYH11","POSTN",                # MyoFibroblast
              # Smooth Muscle
              "CSPG4","MCAM","NES",                                 # Pericyte
              "ICAM1","PCAM1","SELE","VWF","KDR",                   # Blood Endothelial
              "LYVE1","PPROX1","PDPN","FLT4",                       # Lymph Endothelial
              "EPCAM","KRT8","FOXJ1","OVGP1",                       # Epithelial
              "PTPRC",                                              # Hematopoietic 
              "CD8A","CD4","CD3E",                                  # T Cell
              "KLRC1","RUNX3",                                      # T/NK Cell
              "KLRD1","TNFRSF8","NCAM1",                            # NK Cell
              "CD38","JCHAIN",                                      # B Cell
              "FCGR3A","AIF1","CD68","CSF1R","ITGAX","MSLN")               # Macrophage

DotPlot(object = HQ_Cells,                    # Seurat object
        assay = 'RNA',                        # Name of assay to use.  Default is the active assay
        features = features,                  # List of features (select one from above or create a new one)
        # Colors to be used in the gradient
        col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
        col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
        dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
        dot.scale = 6,                        # Scale the size of the points
        group.by = NULL,              # How the cells are going to be grouped
        split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
        scale = TRUE,                         # Whether the data is scaled
        scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
        scale.min = NA,                       # Set lower limit for scaling
        scale.max = NA                        # Set upper limit for scaling
)+    labs(x = NULL,                              # x-axis label
           y = NULL)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 12 , face = "italic"))+
  theme(axis.text.y = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))


# Stacked Bar Graph 2 #

table <- table(HQ_Cells@active.ident ,
               HQ_Cells@meta.data$Sample.ID)    # Create a table of counts

IDs = unique(HQ_Cells@meta.data$Sample.ID)

df <- data.frame(table) 

table2 <- table(HQ_Cells@meta.data$Sample.ID)

library(forcats)

ggplot(data = df2,                # Dataset to use for plot.  Needs to be a data.frame  
       aes(x = Var1,              # Variable to plot on the x-axis
           y = Freq,              # Variable to plot on the y-axis
           fill = factor(Var2,    # Variable to fill the bars
                         levels = IDs),
          )) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", IDs,
                    values = c(natparks.pals(name="Arches2",n=12),
                               natparks.pals(name="Glacier",n=10),
                               natparks.pals(name="WindCave",n=12),
                               natparks.pals(name='Denali', n=20, direction = -1))) +
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis


## AGE TABLE ##

table <- table(HQ_Cells@active.ident ,
               HQ_Cells@meta.data$Age)    # Create a table of counts

IDs = sort(unique(HQ_Cells@meta.data$Age))


df <- data.frame(table) 

table2 <- table(HQ_Cells@meta.data$Age)

colors_blue <- c("#00008B", "#0000CD","#4169E1",  "#1E90FF", "#6495ED", "#87CEEB", "#87CEFA", "#B0E0E6")
colors_purple <- c("#E6E6FA", "#D8BFD8", "#DA70D6", "#BA55D3", "#9932CC","purple")
colors_red <- c( "#FA8072", "red","#B22222")

colors <- c(colors_blue , colors_purple , colors_red)

ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
       aes(x = Var1,              # Variable to plot on the x-axis
           y = Freq,              # Variable to plot on the y-axis
           fill = factor(Var2,    # Variable to fill the bars
                         levels = as.character(IDs)))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", IDs,
                    values = c(colors)) +
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis


# Stacked Bar Graph Location #

table <- table(HQ_Cells@active.ident ,
               HQ_Cells@meta.data$Location)    # Create a table of counts

IDs = unique(HQ_Cells@meta.data$Location)

df <- data.frame(table) 

table2 <- table(HQ_Cells@meta.data$Location)


ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
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
                    values = c('red','maroon','blue','purple','orange','grey'))+
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis


ggsave(filename = "20250402_krt5_age_human_UMAP.pdf" , plot = epi_human , width = 15, height = 12, dpi =600)



## Feature Plots ##

FeaturePlot(HQ_Cells, features = 'TP53', raster= F)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))



features <- c("VIM","COL1A1","GATA2",                               # Fibroblast
              "DES","PALLD","ACTA2","MYH11","POSTN",                # MyoFibroblast
              # Smooth Muscle
              "CSPG4","MCAM","NES",                                 # Pericyte
              "ICAM1","PCAM1","SELE","VWF","KDR",                   # Blood Endothelial
              "LYVE1","PPROX1","PDPN","FLT4",                       # Lymph Endothelial
              "EPCAM","KRT8","FOXJ1","OVGP1",                       # Epithelial
              "PTPRC",                                              # Hematopoietic 
              "CD8A","CD4","CD3E",                                  # T Cell
              "KLRC1","RUNX3",                                      # T/NK Cell
              "KLRD1","TNFRSF8","NCAM1",                            # NK Cell
              "CD38","JCHAIN",                                      # B Cell
              "FCGR3A","AIF1","CD68","CSF1R","ITGAX", "MKI67")               # Macrophage

DotPlot(object = HQ_Cells,                    # Seurat object
        assay = 'RNA',                        # Name of assay to use.  Default is the active assay
        features = features,                  # List of features (select one from above or create a new one)
        # Colors to be used in the gradient
        col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
        col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
        dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
        dot.scale = 6,                        # Scale the size of the points
        group.by = NULL,              # How the cells are going to be grouped
        split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
        scale = TRUE,                         # Whether the data is scaled
        scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
        scale.min = NA,                       # Set lower limit for scaling
        scale.max = NA                        # Set upper limit for scaling
)+    labs(x = NULL,                              # x-axis label
           y = NULL)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 12 , face = "italic"))+
  theme(axis.text.y = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))




#### RENAME ####


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




table <- data.frame(FT_Named@meta.data$Sample,FT_Named@active.ident) 

write.csv(table, file = "20250115_all_human_labels_for_SAMap.csv")



features <- c("VIM","COL1A1","GATA2",                               # Fibroblast
              "DES","PALLD","ACTA2","MYH11","POSTN",                # MyoFibroblast
              # Smooth Muscle
              "CSPG4","MCAM","NES",                                 # Pericyte
              "ICAM1","SELE","VWF","KDR",                   # Blood Endothelial
              "LYVE1","PDPN","FLT4",                       # Lymph Endothelial
              "EPCAM","KRT8","FOXJ1","OVGP1",                       # Epithelial
              "PTPRC",                                              # Hematopoietic 
              "CD8A","CD4","CD3E",                                  # T Cell
              "KLRC1","RUNX3",                                      # T/NK Cell
              "KLRD1","TNFRSF8","NCAM1",                            # NK Cell
              "CD38","JCHAIN",                                      # B Cell
              "FCGR3A","AIF1","CD68","CSF1R","ITGAX", "MKI67")               # Macrophage


features <- c( "SPP1","FKBP5",
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


all_dp <- DotPlot(object = FT_Named,                    # Seurat object
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
  




ggsave(filename = "FIG1f_all_distal_human_dotplot.pdf", plot = all_dp, width = 10, height = 10, dpi = 600)










#### Distal Epithelial Subset ####


FeaturePlot(FT_Named, features = c("EPCAM"), pt.size = 0.75 , raster = F)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

Epi_Filter <- subset(FT_Named, 
                     idents = c("Secretory Epithelial 1", 
                                "Secretory Epithelial 2",
                                "Secretory Epithelial 3",
                                "Ciliated Epithelial"))

Idents(Epi_Filter) <- "Location"
Distal_Epi_Filter <- subset(Epi_Filter, idents = "Fimbria")



# Need to Normalize and Run PCA before Running Harmony ##

Distal_Epi_Filter <- NormalizeData(object = Distal_Epi_Filter, assay = 'RNA')
Distal_Epi_Filter <- FindVariableFeatures(object = Distal_Epi_Filter, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
Distal_Epi_Filter <- ScaleData(object = Distal_Epi_Filter, assay = 'RNA')
Distal_Epi_Filter <- RunPCA(object = Distal_Epi_Filter, assay = 'RNA', features = VariableFeatures(object = Distal_Epi_Filter),
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

pcs <- npcs(seu = Distal_Epi_Filter, var.toal = 0.95, reduction = 'pca_RNA')

## Run Harmony ##

gc()

Distal_Epi_Filter <- RunHarmony(Distal_Epi_Filter, 
                         group.by.vars = c("Sample.ID"), 
                         reduction = "pca_RNA", 
                         assay = "RNA", 
                         plot_convergence = TRUE)

## Normal Seurat clustering ##

Distal_Epi_Filter <- FindNeighbors(object = Distal_Epi_Filter, 
                            reduction = "harmony",
                            k.param = 70, 
                            dims = 1:pcs, 
                            graph.name = 'harmony_snn')

Distal_Epi_Filter <- FindClusters(object = Distal_Epi_Filter, 
                           resolution = 0.4, 
                           graph.name = 'harmony_snn')


## Set Python UMAP via reticulate ##

umap.method = 'umap-learn'
metric = 'correlation'

Distal_Epi_Filter <- RunUMAP(object = Distal_Epi_Filter, 
                      reduction = "harmony", dims = 1:pcs)




saveRDS(Distal_Epi_Filter, file = "20241216_Distal_Epi_Cells.rds")


Distal_Epi_Filter <- readRDS(file = "20241216_Distal_Epi_Cells.rds", refhook = NULL)

## UMAP ##

DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "Source", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "Sample.ID", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')

FeaturePlot(Distal_Epi_Filter, features = c("CD44"), pt.size = 0.75)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

DistalEpiMarkers <- FindAllMarkers(Distal_Epi_Filter, 
                               max.cells.per.ident = 5000)

write.csv(DistalEpiMarkers, file = "20241213_Distal_Epi_All_Markers.csv")


table <- data.frame(Distal_Epi_Filter@meta.data$Sample,Distal_Epi_Filter@active.ident) 

write.csv(table, file = "20250115_epi_human_labels_for_SAMap.csv")


## AGE TABLE ##

table <- table(Distal_Epi_Filter@active.ident ,
               Distal_Epi_Filter@meta.data$Age)    # Create a table of counts

IDs = sort(unique(Distal_Epi_Filter@meta.data$Age))


df <- data.frame(table) 

table2 <- table(Distal_Epi_Filter@meta.data$Age)

colors_blue <- c("#0000CD", "#1E90FF", "#6495ED", "#87CEFA")
colors_purple <- c("#D8BFD8", "#DA70D6", "#BA55D3", "#9932CC","purple")
colors_red <- c( "#FA8072", "red","#B22222")

colors <- c(colors_blue , colors_purple , colors_red)

ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
       aes(x = Var1,              # Variable to plot on the x-axis
           y = Freq,              # Variable to plot on the y-axis
           fill = factor(Var2,    # Variable to fill the bars
                         levels = as.character(IDs)))) + # Order of the stacked bars
  theme_classic() +               # ggplot2 theme
  # Bar plot
  geom_bar(position = 'fill',     # Position of bars.  Fill means the bars are stacked.
           stat = "identity",     # Height of bars represent values in the data
           size = 1) +            # Size of bars
  # Color scheme
  scale_fill_manual("Location", IDs,
                    values = c(colors)) +
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis


# Stacked Bar Graph Location #

table <- table(HQ_Cells@active.ident ,
               HQ_Cells@meta.data$Location)    # Create a table of counts

IDs = unique(HQ_Cells@meta.data$Location)

df <- data.frame(table) 

table2 <- table(HQ_Cells@meta.data$Location)


ggplot(data = df,                # Dataset to use for plot.  Needs to be a data.frame  
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
                    values = c('red','maroon','blue','purple','orange','grey'))+
  labs(x = NULL,                     # x-axis label
       y = "Fraction of Cells") +    # y-axis label
  theme(text = element_text(size = 15),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 60, hjust = 1, size = 11),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1))              # Text color and horizontal adjustment on y-axis


#### Filter 2 ####


FeaturePlot(Distal_Epi_Filter, features = c("LMNTD1"), pt.size = 0.75 , raster = F)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

Distal_Epi_Filter2 <- subset(Distal_Epi_Filter, 
                     idents = c("0", 
                                "1",
                                "2",
                                "3",
                                "5",
                                "7"))

# Need to Normalize and Run PCA before Running Harmony ##

Distal_Epi_Filter2 <- NormalizeData(object = Distal_Epi_Filter2, assay = 'RNA')
Distal_Epi_Filter2 <- FindVariableFeatures(object = Distal_Epi_Filter2, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
Distal_Epi_Filter2 <- ScaleData(object = Distal_Epi_Filter2, assay = 'RNA')
Distal_Epi_Filter2 <- RunPCA(object = Distal_Epi_Filter2, assay = 'RNA', features = VariableFeatures(object = Distal_Epi_Filter2),
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

pcs <- npcs(seu = Distal_Epi_Filter, var.toal = 0.95, reduction = 'pca_RNA')

## Run Harmony ##

gc()

Distal_Epi_Filter2 <- RunHarmony(Distal_Epi_Filter2, 
                                group.by.vars = c("Sample.ID"), 
                                reduction = "pca_RNA", 
                                assay = "RNA", 
                                plot_convergence = TRUE)

## Normal Seurat clustering ##

Distal_Epi_Filter <- FindNeighbors(object = Distal_Epi_Filter, 
                                   reduction = "harmony",
                                   k.param = 70, 
                                   dims = 1:pcs, 
                                   graph.name = 'harmony_snn')

Distal_Epi_Filter <- FindClusters(object = Distal_Epi_Filter, 
                                  resolution = 0.7, 
                                  graph.name = 'harmony_snn')


## Set Python UMAP via reticulate ##

umap.method = 'umap-learn'
metric = 'correlation'

Distal_Epi_Filter <- RunUMAP(object = Distal_Epi_Filter, 
                             reduction = "harmony", dims = 1:pcs)


DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')




saveRDS(Distal_Epi_Filter, file = "20250123_Distal_Epi_Cells.rds")








DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "seurat_clusters", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "Source", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F) + 
  labs(title = 'Colored by Cluster')

DimPlot(object = Distal_Epi_Filter, 
        reduction = 'umap', 
        group.by = "Sample.ID", 
        repel = TRUE, 
        label = TRUE, 
        pt.size = 0.5, 
        label.size = 5,
        raster = F,
        shuffle = T) + 
  labs(title = 'Colored by Cluster')

FeaturePlot(Distal_Epi_Filter2, features = c("KRT5"), pt.size = 0.75)+
  scale_color_viridis_c(option="F",begin=.4,end=0.99, direction = -1)+
  theme(plot.title = element_text(size = 32,face = "bold.italic"))

DistalEpiMarkers <- FindAllMarkers(Distal_Epi_Filter, 
                                   max.cells.per.ident = 5000)

write.csv(DistalEpiMarkers, file = "20250123_Distal_Epi_All_Markers3.csv")




stained_features <- c("CD44","CDKN2A","CTNNB1","ITGA6","KRT5","KRT6B","KRT6A",
                      "KRT7","LEF1","MKI67","PAX2","PAX8","PROM1","TP53","WT1", "MSLN", "UPK3B" , "LRRN4")


features <- c('UPK1A', "SPDEF", "OVGP1", "GSTM2", "SELENOP", "MSLN", "SLC1A3",
             "ITGA6", "PAX8",'ECRG4', 'SOX5', 'PDE4B', 'LCN2','KLF6',
             'TP53' , 'TP73','KRT5','FOXA2','PROM1','CLSTN2','SPEF2','DNAH12','FOXJ1', 'FAM166C' , 'CFAP126',
             'FAM183A' , 'LEF1','MAPK8' , 'EGR1')

DotPlot(object = Distal_Epi_Filter,                    # Seurat object
        assay = 'RNA',                        # Name of assay to use.  Default is the active assay
        features = features,                  # List of features (select one from above or create a new one)
        # Colors to be used in the gradient
        col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
        col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
        dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
        dot.scale = 6,                        # Scale the size of the points
        group.by = NULL,              # How the cells are going to be grouped
        split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
        scale = TRUE,                         # Whether the data is scaled
        scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
        scale.min = NA,                       # Set lower limit for scaling
        scale.max = NA                        # Set upper limit for scaling
)+    labs(x = NULL,                              # x-axis label
           y = NULL)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+ 
  theme(axis.text.x = element_text(size = 12 , face = "italic"))+
  theme(axis.text.y = element_text(size = 12))+
  theme(legend.title = element_text(size = 12))




DistalEpiMarkers2 <- FindAllMarkers(Distal_Epi_Filter2, 
                                   max.cells.per.ident = 5000)

write.csv(DistalEpiMarkers2, file = "20241217_Distal_Epi_All_Markers2.csv")





Distal_Epi_Filter <- readRDS(file = "20241216_Distal_Epi_Cells2.rds", refhook = NULL)



























Distal_Epi_Filter3 <- subset(Distal_Epi_Filter2, 
                             idents = c("0", 
                                        "1",
                                        "3",
                                        "4",
                                        "7"))



Idents(Distal_Epi_Filter3) <- "Classify_Doublets"
Distal_Epi_Filter4 <- subset(Distal_Epi_Filter3, idents = "Singlet")

Idents(Distal_Epi_Filter4) <- "seurat_clusters"



Distal_Epi_Filter5 <- subset(Distal_Epi_Filter, 
                             idents = c("0", 
                                        "2",
                                        "3",
                                        "4",
                                        "7"))

# Need to Normalize and Run PCA before Running Harmony ##

Distal_Epi_Filter2 <- NormalizeData(object = Distal_Epi_Filter2, assay = 'RNA')
Distal_Epi_Filter2 <- FindVariableFeatures(object = Distal_Epi_Filter2, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
Distal_Epi_Filter2 <- ScaleData(object = Distal_Epi_Filter2, assay = 'RNA')
Distal_Epi_Filter2 <- RunPCA(object = Distal_Epi_Filter2, assay = 'RNA', features = VariableFeatures(object = Distal_Epi_Filter2),
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

pcs <- npcs(seu = Distal_Epi_Filter2, var.toal = 0.95, reduction = 'pca_RNA')

## Run Harmony ##

gc()

Distal_Epi_Filter2 <- RunHarmony(Distal_Epi_Filter2, 
                                 group.by.vars = c("Sample.ID"), 
                                 reduction = "pca_RNA", 
                                 assay = "RNA", 
                                 plot_convergence = TRUE)

## Normal Seurat clustering ##

Distal_Epi_Filter2 <- FindNeighbors(object = Distal_Epi_Filter2, 
                                    reduction = "harmony",
                                    k.param = 70, 
                                    dims = 1:pcs, 
                                    graph.name = 'harmony_snn')

Distal_Epi_Filter2 <- FindClusters(object = Distal_Epi_Filter2, 
                                   resolution = 0.4, 
                                   graph.name = 'harmony_snn')


## Set Python UMAP via reticulate ##

umap.method = 'umap-learn'
metric = 'correlation'

Distal_Epi_Filter2 <- RunUMAP(object = Distal_Epi_Filter2, 
                              reduction = "harmony", dims = 1:pcs)



