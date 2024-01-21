#!/bin/bash

# Specify the directory containing the individual folders
# main_directory="/path/to/main_directory"
main_directory="/project/awlab/wuhuiyun/pav_mgs_2023/results/assembly/"


# Loop through each folder
for folder in "$main_directory"/*/; do
    # Extract the sample ID from the folder name
    sample_id=$(basename "$folder")

    # Rename contigs.fasta to include the sample ID
    if [ -e "$folder/contigs.fasta" ]; then
        mv "$folder/contigs.fasta" "$folder/${sample_id}.contigs.fasta"
        echo "Renamed contigs.fasta in $folder"
    else
        echo "contigs.fasta not found in $folder"
    fi
done
