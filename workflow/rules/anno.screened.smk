SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 ".split()
print(SAMPLES)

rule all:
    input:
        file = expand("results/annotation/{sample}_screened.csv", sample=SAMPLES),
        # file.2 = "results/combined.csv"

rule screened_annotation:
    input:
        script = "workflow/scripts/anno.screen.R",
        vs2_file="results/vs2/{sample}.vs2.final-viral-score.tsv",
        checkV_file="results/checkV/{sample}.checkv.quality_summary.tsv",
        diamond_file="results/diamond/{sample}.tsv",
    output:
        individual="results/annotation/{sample}_screened.csv",
        # file.2 = "results/combined.csv"
    script:
        "../scripts/anno.screen.R"
