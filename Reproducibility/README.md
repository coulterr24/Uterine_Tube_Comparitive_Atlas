# Data Download, Alignment, Processing, and Comparison Pipeline 

### Single Cell RNA-Sequencing (scRNA-seq) Data Download
In addition to our own collected mouse [GSE252786](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252786) and human [GSE324855](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi) samples available for download. We downloaded human data from the following sources:
   +  Dinh et al [GSE151214](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151214)
   +  Ulrich et al [GSE178101](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178101)
   +  Lengyel et al [EGAC00001003114](https://ega-archive.org/datasets/EGAD00001010076) <br>

To download large datasets at once, using the Eurpean Nucleotide Archive (ENA) enables the 'wget' function for accession numbers from BioProject IDs found in the Gene Expression Omnibus (GEO). Pseudocode for downloading a large project with multiple scRNA-seq samples is provided as `ENA_Downloader.sh`. <br><br>

### Raw Data Alignment
**CellRanger** (v.7.1.0, 10x Genomics): [tutorial](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct) <br>
Within a linux terminal and with CellRanger downloaded with the proper file paths, run the count command to process raw scRNA-seq reads <br>
`cellranger count --id=run_count_1kpbmcs \
   --fastqs=/path_to/uterine_tube1_fastqs \
   --sample=uterine_tube1 \
   --transcriptome=/path_to/tanscriptome` <br> <br>

Note: This is computationally intensive for most standard computers. If you are on a local PC, I would advise starting from the processed gene count matrices or the prepared Seurat Objects.


### Setting Up R and RStudio
Install R (v.4.2.0) and Rstudio (https://rstudio-education.github.io/hopr/starting.html): approximately 15-30 min. <br> <br>

Install packages found in the `RSessionInfo.txt` with proper versions to properly recreate the figures: approximately 90 min. <br> <br>

### Setting Up python and jupyter notebook
Install python (v.3.10.8) and jupyter notebook (https://jupyter.org/install): approximately 15-30 min. <br> <br>

Install packages found in the `pySessionInfo.txt` with proper versions to properly recreate the figures: approximately 30 min. <br> <br>

### Preprocessing Mouse and Human Objects
Refer to code in `coulterr24/MouseTE_scRNA/tree/main/Data_Preprocessing/` to reproduce preproccessing results for mice. <br>
Refer to code in `coulterr24/Uterine_Tube_Comparitive_Atlas/tree/main/Seurat_Preprocessing/` to reproduce preproccessing results for human. <br> <br>

#### Mouse Objects for All Cells
Use `MouseTE_scRNA_Metadata.csv` from `coulterr24/MouseTE_scRNA/tree/main/Data_Preprocessing/`to identify sample file directories and additional meta data. <br>
`SC_Mouse_Distal_Data_Preparation.R` and `SC_Mouse_Distal_Data_Preparation.R` will go step by step for removing ambient RNA signals with **SoupX** (v.1.6.1), identifying doublets with **DoubletFinder** (v.2.0.3), and preparing integrated/batch-corrected objects with **Seurat** (v.4.1.1) and **harmony** (v.0.1.0).
#### Human Objects for All Cells
Use `Human_FT_Healthy_meta.csv` found in `coulterr24/Uterine_Tube_Comparitive_Atlas/tree/main/Seurat_Preprocessing/` to identify sample file directories and additional meta data. <br>
`All_Healthy_Human_FT_PreProcess.R` and `Distal_Epi_Healthy_Human_PreProcess_FT.R` will go step by step for removing ambient RNA signals with **SoupX** (v.1.6.2), identifying doublets with **DoubletFinder** (v.2.0.3), and preparing integrated/batch-corrected objects with **Seurat** (v.4.3.0) and **harmony** (v.0.1.0).
#### Subsetting Epithelial Cells from Mouse Objects
Within `coulterr24/MouseTE_scRNA/tree/main/Data_Preprocessing/`, `SC_Mouse_Distal_Epithelial_Preparation.R` and `SC_Mouse_Proximal_Epithelial_Preparation.R` will go further to subset epithelial cells. <br>
Epithelial cells were identified by clusters that expressed *Krt8* and *Epcam* (common epithelial gene markers).<br>
After characterization of epithelial cell states, **phateR** (v.1.0.7) and **Monocle3** (v.1.2.7) were used in the distal dataset to infer differentiation trajectories withen the distal tubal epithelium. <br>
Final objects for all mouse datasets were created within the preprocessing scripts. These scripts are also available for download on [Dryad](https://doi.org/10.5061/dryad.4mw6m90hm).
#### Subsetting Epithelial Cells from Human Objects
Within  `coulterr24/Uterine_Tube_Comparitive_Atlas/tree/main/Seurat_Preprocessing/`, `SC_Mouse_Distal_Epithelial_Preparation.R` and `SC_Mouse_Proximal_Epithelial_Preparation.R` will go further to subset epithelial cells. <br>
Epithelial cells were identified by clusters that expressed *KRT8* and *EPCAM* (common epithelial gene markers).<br>
After characterization of epithelial cell states, **phateR** (v.1.0.7) and **Monocle3** (v.1.0.0) were used in the distal dataset to infer differentiation trajectories withen the distal tubal epithelium. <br>
Final objects for all human datasets were created within the preprocessing scripts. These scripts are also available for download on [Dryad](https://doi.org/10.5061/dryad.zs7h44jr5) [active upon publication].

### Characterization of Human and Mouse Objects in Seurat
Now that the objects are available, try recreating the figures from the manuscript with these datasets as the input. Full guides can be followed for each figure from the scripts found in `coulterr24/Uterine_Tube_Comparitive_Atlas/tree/main/R_Figure_Scripts/`. <br>

Each figure should take only about 2-10 minutes to recreate with the properly loaded objects and packages <br>

If there are any issues with certain packages, refer to 'RsessionInfo.txt' to ensure that versions match. <br>

### Generation of Human and Mouse Objects in SAMap
With both human and mouse Seurat objects available, we can move to python-based SAMap for cross-species mapping. Downloading and preparation of materials for SAMap was conducted as the pipeline states at (https://github.com/atarashansky/SAMap).

The initial BLAST mapping between mouse and human trascriptomes was computed with `tblastx` for mouse and human transcriptomes. This step can be computationally demanding taking about 3 days with larger annotated transcriptomes, and the output consists of two files each about 3 GB.

With **sam-algorithm** (v.1.0.2) and **samap** (v.1.0.12), the scRNA-seq data from humans and mice were prepared with there Seurat annotations based on the code provided in `coulterr24/Uterine_Tube_Comparitive_Atlas/tree/main/SAMap_Preprocessing_and_Figures`. This step takes about 3-4 hours to fully complete. To skip this step, prepared SAMap objects are avaiable for download at [Dryad](https://doi.org/10.5061/dryad.zs7h44jr5) [active upon publication].

The files in `coulterr24/Uterine_Tube_Comparitive_Atlas/tree/main/SAMap_Preprocessing_and_Figures` also enables figure generation. After downloading the prepared SAMap file, step-by-step code will allow for figure recreation within 10-30 minutes.

### Comparison between Human and Mouse Epithelial Cells based on SAMap
The preprocessing code includes generation of upregulated genes for each mouse and human epithelial cell state. These data were exported as csv files and are available in `coulterr24/Uterine_Tube_Comparitive_Atlas/tree/main/Correlation_Analysis`. Here, SAMap correlation scores were calculated for each gene and output as `20251205_Conservation_Reference_Genes_for MvH_DP.xlsx`. THese data can be interrogated for further correlated gene detection with promise of translatability.


