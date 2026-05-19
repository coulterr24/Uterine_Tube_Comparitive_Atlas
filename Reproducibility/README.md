# Data Processing and Analysis Pseudocode Pipeline

### Raw Data Alignment
**CellRanger** (v6.1.2, 10x Genomics): [tutorial](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct) <br>
Within a linux terminal and with CellRanger downloaded with the proper file paths, run the count command to process raw scRNA-seq reads <br>
`cellranger count --id=run_count_1kpbmcs \
   --fastqs=/path_to/uterine_tube1_fastqs \
   --sample=uterine_tube1 \
   --transcriptome=/path_to/mm39_ref_genome` <br> <br>

Note: This is computationally intensive for most standard computers. If you are on a local PC, I would advise starting from the processed gene count matrices or the prepared Seurat Objects

### Setting Up R and RStudio
Install R (v.4.1.1) and Rstudio (https://rstudio-education.github.io/hopr/starting.html): approximately 15-30 min <br> <br>

Install packages found in the 'sessionInfo.txt' with proper versions to properly recreate the figures: approximately 60 min <br> <br>

### Preprocessing Distal and Proximal Objects
Refer to code in `MouseTE_scRNA/Data_Preprocessing/` to reproduce preproccessing results <br> <br>
#### Distal and Proximal Objects for All Cells
Use `MouseTE_scRNA_Metadata.csv` to identify sample file directories and additional meta data <br>
`SC_Mouse_Distal_Data_Preparation.R` and `SC_Mouse_Distal_Data_Preparation.R` will go step by step for removing ambient RNA signals with **SoupX** (v.1.6.1), identifying doublets with **DoubletFinder** (v.2.0.3), and preparing integrated/batch-corrected objects with **Seurat** (v.4.1.1) and **harmony** (v.0.1.0).
#### Subsetting Epithelial Cells from Distal and Proximal Objects
`SC_Mouse_Distal_Epithelial_Preparation.R` and `SC_Mouse_Proximal_Epithelial_Preparation.R` will go further to subset epithelial cells for furhter analysis. <br>
Epithelial cells were identified by clusters that expressed *Krt8* and *Epcam* (common epithelial gene markers).<br>
After characterization of epithelial cell states, **phateR** (v.1.0.7) and **Monocle3** (v.1.2.7) were used in the distal dataset to infer differentiation trajectories withen the distal tubal epithelium. <br>
Final objects for all datasets were created within the preprocessing scripts. THese scripts are also available for download on [Dryad](https://doi.org/10.5061/dryad.4mw6m90hm) [active upon publication].

### Exploration of Final Distal and Proximal Objects
Now that the objects are available, try recreating the figures from the manuscript with these datasets as the input. Full guides can be followed for each figure from the scripts found in `MouseTE_scRNA/Figure_Scripts/`. <br>

Each figure should take only about 2-10 minutes to recreate with the properly loaded objects and packages <br>

If there are any issues with certain packages, refer to 'sessionInfo.txt' to ensure that versions match. <br>

Feel free to explore beyond the figures and delve into other cell types as well :)
