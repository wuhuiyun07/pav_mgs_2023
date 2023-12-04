# SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()

# SAMPLES, = glob_wildcards("rawdata/{sample}_L001_R1_001.fastq.gz")

# rule all:
#     input: expand("results/trimmed/{sample}.html", sample=SAMPLES)

# import pandas as pd

samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)

# SAMPLES, = glob_wildcards("rawdata/{sample}_L001_R1_001.fastq.gz")

# dataset = rawdata
# group = sample


rule all:
    input: expand("results/trimmed/{sample}.html", sample=SAMPLES)
print(SAMPLES)

rule fastp_pe:
    input:
        sample=["rawdata/{sample}_L001_R1_001.fastq.gz", "reads/pe/{sample}L001_R2_001.fastq.gz"]
    output:
        trimmed=["results/trimmed/{sample}.R1.fastq.gz", "results/trimmed/{sample}.R2.fastq.gz"],
        html="reports/fastp/{sample}.html",
        json="reports/fastp/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    threads: 2
    wrapper:
        "v3.0.2/bio/fastp"

