#!/bin/bash
#SBATCH --job-name=pav-mgs 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=72:00:00
#SBATCH --output=slurm-%j.out-%N  
#SBATCH --error=slurm-%j.err-%N 
#SBATCH --account=loni_virus2023      # your account name
#SBATCH --partition=single             # the partition

# export WORK_DIR=/project/awlab/wuhuiyun/pav_mgs_2023
# Load any required modules for your HPC.

export OMP_NUM_THREADS=48
singularity run -B /project resources/sifs/virsorter2.2.4.sif virsorter run --prep-for-dramv -w results/virsorter/16_3_S3 -i results/assembly/test/16_3_S3/contigs.fasta --min-length 1500 -j 4 all"
