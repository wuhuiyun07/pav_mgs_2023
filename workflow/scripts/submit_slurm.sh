#!/bin/bash
#SBATCH --job-name=pav-mgs 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=96:00:00
#SBATCH --output=log/hpc/slurm-%j.out-%N  
#SBATCH --error=log/hpc/slurm-%j.err-%N 
#SBATCH --account=loni_virus2023      # your account name
#SBATCH --partition=single             # the partition

export WORK_DIR=/project/awlab/wuhuiyun/pav_mgs_2023
# Load any required modules for your HPC.

# Run snakemake
snakemake --profile config/slurm --latency-wait 90 --use-singularity --use-conda --configfile config/test.yaml