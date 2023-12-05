
import pandas as pd
wildcard_constraints:
    dataset="\d+"

# samples: config/samples-template.tsv

samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)

rule all:
    input: expand("reports/assembly/{sample}.spades.txt", sample=SAMPLES)
container: "docker://continuumio/miniconda3:4.4.10"

rule run_metaspades:
    input:
        reads=["results/test/{sample}.R1.fastq.gz", "results/test/{sample}.R2.fastq.gz"],
    output:
        contigs="results/assembly/{sample}.contigs.fasta",
        scaffolds="results/assembly/{sample}.scaffolds.fasta",
        dir=directory("results/assembly/{sample}_intermediate_files"),
    benchmark:
        "reports/assembly/{sample}.spades.txt"
    params:
        # all parameters are optional
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/{sample}.spades.log",
    # container:
    #     "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0"
    threads: 8
    resources:
        mem_mem=250000,
        time=60 * 24,
    wrapper:
        "v2.11.1/bio/spades/metaspades"