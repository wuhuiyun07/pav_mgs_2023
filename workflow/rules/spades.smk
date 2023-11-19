container: "docker://continuumio/miniconda3:4.4.10"
rule run_metaspades:

    input:
        reads=["results/trimmed/{sample}_R1.fastq.gz", "results/trimmed/{sample}_R2.fastq.gz"],
    output:
        contigs="results/assembly/{sample}_contigs.fasta",
        scaffolds="results/assembly/{sample}_scaffolds.fasta",
        dir=directory("results/assembly/{sample}_intermediate_files"),
    benchmark:
        "reports/assembly/{sample}_spades.txt"
    params:
        # all parameters are optional
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/{sample}_spades.log",
    container:
        "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0"
    threads: 8
    resources:
        mem_mem=250000,
        time=60 * 24,
    # wrapper:
    #     "v2.11.1/bio/spades/metaspades"