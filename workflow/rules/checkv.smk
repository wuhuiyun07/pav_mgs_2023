# rule all_checkV:
#     input:
#         merged = "results/checkV/merged_checkV_allcontigs.tsv"

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
        dir ="results/checkV/{sample}"
    shell:
        """
        check end_to_end {input.fasta} {output.dir} -t {params.threads} -d {params.db}
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


