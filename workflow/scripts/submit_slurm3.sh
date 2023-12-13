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

# snakemake --use-conda  \ 
#     # --configfile "$CONFIG" \
#     --latnecy-wait 60 \
#     --core 48 \
#     # --cluster "xenon -vvv scheduleer $SCH " \
#     --location local:// submit \
#     --working-directory . \
#     --jobs 48 \
#     --cpus-per-task 48\
#     -s workflow/rules/spades2.smk \
#     # --name smk.{rule} \ 
#     # --inherit-env \
#     --max-run-time 15 \
#     --working-directory . \
#     --stderr stderr-%j.log \
#     --stdout stdout-%j.log"


snakemake --use-conda  -c6 -j6 -s workflow/rules/spades2.smk

snakemake --use-singularity -c6 -j6 -s workflow/rules/spades2.smk