import pandas as pd
wildcard_constraints:
    dataset="\d+"

# samples: config/samples-template.tsv

samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)

rule all:
    input: expand("reports/assembly/{sample}.spades2.txt", sample=SAMPLES)

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/spades/metaspades.html


rule run_metaspades:
    input:
        reads=["results/trimmed/test/{sample}.R1.fastq.gz", "results/trimmed/test/{sample}.R2.fastq.gz"],
    output:
        contigs="results/assembly/test2/{sample}.contigs.fasta",
        scaffolds="results/assembly/test2/{sample}.scaffolds.fasta",
        dir=directory("results/assembly/test2/{sample}_intermediate_files"),
    benchmark:
        "reports/assembly/{sample}.spades2.txt" 
    params:
        # all parameters are optional
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/{sample}.spades2.log",
    # container:
    #     "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0"
    threads: 48
    resources:
        mem_mb=250000,
        cpus_per_task=48,
        # mem_mem=250000,
        time=60 * 24
    # container: 
    #     "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0"
    # script:
    #     "../scripts/spades_script.py"
    wrapper:
        "v3.0.2/bio/spades/metaspades"
    


