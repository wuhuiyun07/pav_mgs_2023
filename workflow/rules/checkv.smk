SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 ".split()
print(SAMPLES)

rule all_checkV:
    input:
        expand("results/checkV/{sample}.checkv.quality_summary.tsv", sample = SAMPLES)

rule chg_file_name:
    output:
        checkV = "results/checkV/{sample}.checkv.quality_summary.tsv"
    input:
        checkV = "results/checkV/{sample}/quality_summary.tsv"
    shell:
        """ scp {input.checkV} {output.checkV} """

rule checkV_spades:
    input:
        fasta ="results/spades.3.15/{sample}/contigs.fasta"
    params:
        threads = "16",
        dir = "results/checkV/{sample}",
        db = "resources/checkv-db-v1.5/"
    container:
        "../sifs/checkv_0.9.0--pyhdfd78af_0.sif"
    output:
        dir = directory("results/checkV/{sample}")
    shell:
        """
        checkv end_to_end {input.fasta} {output.dir} -t {params.threads} -d {params.db}
        """


rule checkv_database:
    output:
        db = "resources/checkv-db-v1.5/"
    container:
        "../sif/checkv_0.9.0--pyhdfd78af_0.sif"
    shell:
        "checkv download_database resources/"



# #to have wildcards in the input of a rule but not in the output of the rule
# def table_inputs(folder, name, wildcards):
#     files=expand("%s{sample}%s" % (folder, name), sample=SAMPLE)
#     return files

# rule merge_checkV:
#     input:
#         i=table_inputs(folder="checkV/", name="/quality_summary.tsv", wildcards=SAMPLE)
#     output:
#         merge="checkV/merged_checkV_allcontigs.tsv"
#     script:
#         "Scripts/merge_checkV.py"


