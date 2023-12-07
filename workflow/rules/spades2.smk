# rule run_metaspades:
#     output:
#         contigs="results/assembly/test2/{sample}.contigs.fasta",
#         scaffolds="results/assembly/test2/{sample}.scaffolds.fasta",
#         dir=directory("results/assembly/test2/{sample}_intermediate_files"),
#     benchmark:
#         "reports/assembly/{sample}.spades2.txt",
#     params:
#         k="auto",
#         extra="--only-assembler",
#     log:
#         "reports/assembly/{sample}.spades2.log",
#     threads: 48,
#     resources:
#         mem_mb=250000,
#         cpus_per_task=48,
#         time=60 * 24,
#     container:
#         "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0"
#     script:
#         "../scripts/spades_script.py"

import os
import json
import glob
import pandas as pd
from snakemake.utils import validate

samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)

rule all:
    input: expand("reports/assembly/{sample}.spades2.log", sample=SAMPLES)

wildcard_constraints:
    dataset="\d+"

rule run_metaspades:
    input:
        reads=["results/trimmed/test/{sample}.R1.fastq.gz", "results/trimmed/test/{sample}.R2.fastq.gz"],
    output:
        contigs="results/assembly/test2/{sample}.contigs.fasta",
        scaffolds="results/assembly/test2/{sample}.scaffolds.fasta",
        dir=directory("results/assembly/test2/{sample}_intermediate_files"),
    params:
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/{sample}.spades2.log",
    threads: 48,
    resources:
        mem_mb=250000,
        cpus_per_task=48,
        time=60 * 24,
    # container:
    #     "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0",
    conda:
        "../envs/spades.yml"
    shell:
        # "singularity exec docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0 "
        "spades.py --meta "
        "--threads {threads} "
        "--memory {resources.mem_mb} "  # Use resources.mem_mb instead of resources.mem
        "-o {output} "
        "-k {params.k} "
        "{input} "
        "{params.extra} "
        "> {log} 2>&1"
