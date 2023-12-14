sample = "24_4_S19"
rule run_metaspades:
    input:
        reads=["results/trimmed/{sample}.R1.fastq.gz", "results/trimmed/{sample}.R2.fastq.gz"],
    output:
        contigs="results/assembly/test/{sample}.contigs.fasta",
        scaffolds="results/assembly/test/{sample}.scaffolds.fasta",
        dir=directory("results/assembly/test/{sample}_intermediate_files"),
    params:
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/{sample}.spades2.log",
    shell:
        """
        spades.py --version
        export OMP_NUM_THREADS=48
        spades.py --meta \\
            -1 results/trimmed/{wildcards.sample}.R1.fastq.gz \\
            -2 results/trimmed/{wildcards.sample}.R2.fastq.gz \\
            -o results/assembly/test/{wildcards.sample}/
        """
