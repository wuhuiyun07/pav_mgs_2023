# rule multiqc:
#     output:
#         mqc_out = directory('results/multiqc_out'),
#         mqc_in  = directory('results/multiqc_in'),
#     input:
#         salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
#         kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
#         fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
#     container:
#         "docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
#     shell:
#         """mkdir {output.mqc_in}
#            ln -snr -t {output.mqc_in} {input}
#            multiqc {output.mqc_in} -o {output.mqc_out}
#         """

SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()

rule multiqc:
    input:
        fastp= expand("report/fastp/{sample}.html", sample = SAMPLES),
    output:
        mqc_out = directory('results/multiqc_out'),
        html=  "result/multiqc_report.html",
    container:
        "docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    shell:
        "multiqc {input.fastp} -o {output.mqc_out}"