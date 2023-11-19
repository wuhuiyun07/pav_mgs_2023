
rule all_virsorter2:
    input:
        all="virsorter2/vs2_merged_file.txt"


rule virsorter2:
    input:
        "results/assembly/{sample}_contigs.fasta"
    output:
        "results/virsorter2/{sample}_final-viral-combined.fa",
        "results/virsorter2/{sample}_final-viral-score.tsv",
        "results/virsorter2/{sample}_final-viral-boundary.tsv"
    params:
        path="results/virsorter2/",
        nodes="16"
    container:
        "docker://jiarong/virsorter:latest"
    shell:
        "virsorter run --prep-for-dramv -w {params.path} -i {input} -j {params.nodes} all"

# #to have wildcards in the input of a rule but not in the output of the rule
# def table_inputs(folder, name, wildcards):
#     files=expand("%s{sample}%s" % (folder, name), sample=SAMPLE)
#     return files

rule merge_vs2:
    input:
        i=table_inputs(folder="virsorter2/", name="/final-viral-score.tsv", wildcards=SAMPLE)
    params:
        s=SAMPLE
    output:
        o="virsorter2/vs2_merged_file.txt"
    script:
        "Scripts/merge_virsorter2.py"


