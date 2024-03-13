#!/usr/bin/env bash

# https://hecatomb.readthedocs.io/en/latest/install/

# download hecatomb by submitting a job on LONI, single node, 8 core, take about 1.5 hours
conda create -n hecatomb -c conda-forge -c bioconda hecatomb

## To use Hecatomb, activate your new conda env.
conda activate hecatomb
# Check that it's installed
hecatomb --help

# Download database
hecatomb install
rule download_db_file:
# output: /ddnA/project/awlab/wuhuiyun/.conda/envs/hecatomb/lib/python3.10/site-packages/hecatomb/snakemake/workflow/../databases/aa/virus_secondary_aa/sequenceDB_nodes.dmp

# HPC profile
## https://hecatomb.readthedocs.io/en/latest/profiles/


# run test dataset
hecatomb test

# pair end reads
hecatomb run --reads pav.test

##
