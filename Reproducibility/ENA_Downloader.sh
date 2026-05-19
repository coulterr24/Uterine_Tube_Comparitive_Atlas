#!/bin/bash

# CSV File with 'Experiment' and 'Run'
# CSV with No Headers and Runs in Col1 and Experiments in Col2

ENA_Download=ENA_Download.csv

# Read the CSV file into arrays
readarray -t sra_info < <(cut -d ',' -f 1,2 $ENA_Download)

# Loop through each study ID and study name
for i in "${!sra_info[@]}"; do
    # Extract the study ID and study name from the current array element
    Run=$(echo "${sra_info[$i]}" | cut -d ',' -f 1)
    Experiment=$(echo "${sra_info[$i]}" | cut -d ',' -f 2)
    
    # Set the output file name
    output_file=./MSKCC/${Experiment}/${Experiment}_S1_L001_R1_001.fastq.gz

    # Check if the output file already exists
    if [[ -f $output_file ]]; then
        # If the file exists, concatenate it with the new file
        echo "Downloading $study_name Read 1 as TMP1"
        #curl -o ./MSKCC/${Experiment}/TMP1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_1.fastq.gz
        wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_1.fastq.gz -O ./MSKCC/${Experiment}/TMP1.fastq.gz
        echo "Concatenating $output_file with $study_name Read 1"
        cat ./MSKCC/${Experiment}/TMP1.fastq.gz >> $output_file && rm ./MSKCC/${Experiment}/TMP1.fastq.gz
        echo "Downloading $study_name Read 2 as TMP2"
        #curl -o ./MSKCC/${Experiment}/TMP2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_2.fastq.gz
        wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_2.fastq.gz -O ./MSKCC/${Experiment}/TMP2.fastq.gz
        echo "Concatenating ./MSKCC/${Experiment}/${Run}_S1_L001_R2_001.fastq.gz with $study_name Read 2"
        cat ./MSKCC/${Experiment}/TMP2.fastq.gz >> ./MSKCC/${Experiment}/${Experiment}_S1_L001_R2_001.fastq.gz && rm ./MSKCC/${Experiment}/TMP2.fastq.gz 
    else
        mkdir ./MSKCC/${Experiment}/
        # If the file doesn't exist, download it
        echo "Downloading $study_name Read 1"
        #curl -o ./MSKCC/${Experiment}/${Run}_S1_L001_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_1.fastq.gz
        wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_1.fastq.gz -O ./MSKCC/${Experiment}/${Experiment}_S1_L001_R1_001.fastq.gz
        echo "Downloading $study_name Read 2 "
        #curl -o ./MSKCC/${Experiment}/${Run}_S1_L001_R1_002.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_2.fastq.gz
        wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${Run:0:6}/0${Run:(-2)}/${Run}/${Run}_2.fastq.gz -O ./MSKCC/${Experiment}/${Experiment}_S1_L001_R2_001.fastq.gz
    fi
done