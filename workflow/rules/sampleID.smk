SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 ".split()
print(SAMPLES)

rule all:
    input:
        file = expand("results/annotation2/{sample}.csv", sample=SAMPLES),

rule add_sampleID:
    input:
        file="results/annotation/{sample}_screened.csv",
    output:
        individual="results/annotation2/{sample}.csv",
    script:
        "../scripts/process_data.R"