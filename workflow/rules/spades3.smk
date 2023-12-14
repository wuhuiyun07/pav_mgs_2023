import pandas as pd
wildcard_constraints:
    dataset="\d+"

# samples: config/samples-template.tsv

# samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
# SAMPLES = samples_df["sample_name"].tolist()
# print(SAMPLES)

rule run_metaspades:
    input:
        reads=["results/trimmed/24_4_S19.R1.fastq.gz", "results/trimmed/24_4_S19.R2.fastq.gz"],
    output:
        contigs="results/assembly/test/24_4_S19.contigs.fasta",
        scaffolds="results/assembly/test/24_4_S19.scaffolds.fasta",
        dir=directory("results/assembly/test/24_4_S19_intermediate_files"),
    params:
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/24_4_S19.spades2.log",
    shell:
        """
        spades.py --version
        export OMP_NUM_THREADS=48
        spades.py --meta\\   
            -1 results/trimmed/24_4_S19.R1.fastq.gz \\
            -2 results/trimmed/24_4_S19.R2.fastq.gz \\
            -o results/assembly/test/
        """