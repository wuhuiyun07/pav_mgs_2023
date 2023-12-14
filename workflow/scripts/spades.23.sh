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

export OMP_NUM_THREADS=48

module load snakemake
module load python
module load spades/3.13.0/gcc-9.3.0
# Run snakemake
# snakemake --profile config/slurm --latency-wait 90 --use-singularity --use-conda --configfile config/test.yaml
# cd /project/awlab/wuhuiyun/pav_mgs_2023
export WORK_DIR=/project/awlab/wuhuiyun/pav_mgs_2023

# SAMPLES = "23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 Blank_8_15_23_S26".split()

 spades.py --version
        export OMP_NUM_THREADS=48
        spades.py --meta \
            --threads 48 \
            -1 results/trimmed/23_1_S11.R1.fastq.gz \
            -2 results/trimmed/23_1_S11.R2.fastq.gz \
            -o results/assembly/test/23_1_S11/

 spades.py --version
        export OMP_NUM_THREADS=48
        spades.py --meta \
            --threads 48 \
            -1 results/trimmed/23_2_S12.R1.fastq.gz \
            -2 results/trimmed/23_2_S12.R2.fastq.gz \
            -o results/assembly/test/23_2_S12/
 spades.py --version
        export OMP_NUM_THREADS=48
        spades.py --meta \
            --threads 48 \
            -1 results/trimmed/23_3_S13.R1.fastq.gz \
            -2 results/trimmed/23_3_S13.R2.fastq.gz \
            -o results/assembly/test/23_3_S13/
             spades.py --version
        export OMP_NUM_THREADS=48

spades.py --meta \
    --threads 48 \
    -1 results/trimmed/23_4_S14.R1.fastq.gz \
    -2 results/trimmed/23_4_S14.R2.fastq.gz \
    -o results/assembly/test/23_4_S14/

 spades.py --version
        export OMP_NUM_THREADS=48
        spades.py --meta \
            --threads 48 \
            -1 results/trimmed/23_5_S15.R1.fastq.gz \
            -2 results/trimmed/23_5_S15.R2.fastq.gz \
            -o results/assembly/test/23_5_S15/

    