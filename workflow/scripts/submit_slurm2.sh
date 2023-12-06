#!/bin/bash
#SBATCH --job-name=pav-mgs 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%j.out-%N  
#SBATCH --error=slurm-%j.err-%N 
#SBATCH --account=loni_virus2023      # your account name
#SBATCH --partition=workq             # the partition

# export WORK_DIR=/project/awlab/wuhuiyun/pav_mgs_2023
# Load any required modules for your HPC.

module load snakemake
module load python
# Run snakemake
# snakemake --profile config/slurm --latency-wait 90 --use-singularity --use-conda --configfile config/test.yaml
# cd /project/awlab/wuhuiyun/pav_mgs_2023
export WORK_DIR=/project/awlab/wuhuiyun/pav_mgs_2023

snakemake -c 48 --use-singularity -j 48 -s workflow/rules/spades2.smk