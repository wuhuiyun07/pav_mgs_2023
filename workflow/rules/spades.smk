
# import pandas as pd
# wildcard_constraints:
#     dataset="\d+"

# # samples: config/samples-template.tsv

# samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
# SAMPLES = samples_df["sample_name"].tolist()
# print(SAMPLES)

SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()
print(SAMPLES)

# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/spades/metaspades.html

# container: "docker://continuumio/miniconda3:4.4.10"
# rule OMP:
#     shell:
#         "export OMP_NUM_THREADS=48"


rule all_spades:
    input: expand("reports/assembly/{sample}.spades.txt", sample=SAMPLES)

rule metaspades:
    input:
        R1 = expand("results/trimmed/{sample}.R1.fastq.gz", sample = SAMPLES),
        R2 = expand("results/trimmed/{sample}.R2.fastq.gz", sample = SAMPLES)
    output:
        contigs="results/assembly/{sample}.contigs.fasta",
        scaffolds="results/assembly/{sample}.scaffolds.fasta",
        dir=directory("results/assembly/{sample}_intermediate_files"),
    benchmark:
        "reports/assembly/{sample}.spades.txt" 
    params:
        # all parameters are optional
        k= "auto",
        extra= "--only-assembler",
    log:
        "reports/assembly/{sample}.spades.log",
    container:
        "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0"
    threads: 48
    resources:
        mem_mem=250000,
        time=60 * 24,
    shell:
        r"""mkdir {output.dir}
        spades.py --meta  -t {threads}   
        -o {output.dir} 
        -k {params.k}   
        --pe1-1 {input.R1}
        --pe1-2 {input.R2}  
        """
    # wrapper:
    #     "v3.0.2/bio/spades/metaspades"
   


