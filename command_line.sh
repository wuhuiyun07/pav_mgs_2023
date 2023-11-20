#!/bin/bash
# transfer data to LONI
rsync --rsh=ssh --archive --stats --progress ./nextseq2000_fastqz \
        wuhuiyun@qbc.loni.org:/project/awlab/wuhuiyun/pav/rawdata
# snakemake fastp singularity
snakemake  -s ./workflow/rules/fastp.smk --use-conda --core 4 -np 

singularity shell -B /project/awlab/wuhuiyun fastp_0.23.3--h5f740d0_0.sif

singularity shell -B /project/awlab/wuhuiyun ./pav_sif/snakemake_7.32.4--hdfd78af_1.sif
snakemake  --use-conda -s ./workflow/rules/fastp.smk --core 8 -np

ls ./rawdata >filenames.txt # extract filenames from the raw data sample
filename="example.txt"
substring=$(echo "$filename" | awk -F'.' '{print $1}')
echo "Substring: $substring"


for file in *.fastq.gz; do 
    file="${file%.*}"
    file="${file##*.}"
    echo "${file:0:8}"
done