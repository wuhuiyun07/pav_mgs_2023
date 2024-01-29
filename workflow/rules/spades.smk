
# import pandas as pd
# wildcard_constraints:
#     dataset="\d+"

# # samples: config/samples-template.tsv

# SAMPLES = samples_df["sample_name"].tolist()
# print(SAMPLES)

SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()
print(SAMPLES)

rule all_spades:
    input: expand("reports/spades.3.13/spades.log", sample = SAMPLES)

rule metaspades:
    input:
        R1 = "results/trimmed/{sample}.R1.fastq.gz",
        R2 = "results/trimmed/{sample}.R2.fastq.gz"
    output:
        dir = directory("results/spades.3.13/{sample}"),
        contigs="results/spades.3.13/{sample}/contigs.fasta",
    params:
        k= "auto",
    log:
        "reports/spades.3.13/{sample}/spades.log",     
    container:
        "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0"
    threads: 48
    resources:
        mem_mem=250000,
        time=60 * 24,
    shell:
        r"""
        spades.py --meta  -o {output.dir} -t {threads}  --pe1-1 {input.R1} --pe1-2 {input.R2}  
        # mv results/assembly/{sample}/contigs.fasta results/assembly/{sample}.contigs.fasta
        """
    # wrapper:
    #     "v3.0.2/bio/spades/metaspades"
   


# code from github metagenome assembly
# https://github.com/metagenome-atlas/metagenome-assembly/blob/main/workflow/rules/spades.smk