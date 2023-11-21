#!/bin/bash
#SBATCH --job-name=pav-mgs # sbatch options here only affect the overall job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
#SBATCH --output=log/hpc/slurm-%j.out-%N
#SBATCH --erro=log/hpc/slurm-%j.err-%N 
#SBATCH --account=wuhuiyun      # your account name
#SBATCH --partition=single             # the partition

export WORK_DIR=/project/awlab/wuhuiyun/pav_mgs_2023
# Load any required modules for your HPC.

# Run snakemake
snakemake --profile config/slurm --latency-wait 90 --use-singularity --use-conda --configfile config/test.yaml