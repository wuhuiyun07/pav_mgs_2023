
# rule all_virsorter2:
#     input:
#         all="virsorter2/vs2_merged_file.txt"


SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()
print(SAMPLES)

rule all_v2:
    input: 
        DNA.score = expand("results/vs2_DNA/{sample}/final-viral-score.tsv", sample=SAMPLES),
        RNA.score = expand("results/vs2_RNA/{sample}/final-viral-score.tsv", sample=SAMPLES)

rule vs2_DNA:
    input:
        "results/assembly/spades3.13/{sample}/contigs.fasta"
    output:
        fa = "results/vs2_DNA/{sample}/final-viral-combined.fa",
        score ="results/vs2_DNA/{sample}/final-viral-score.tsv",
        boundary = "results/vs2_DNA/{sample}/final-viral-boundary.tsv"
    params:
        path = "results/vs2_DNA/{sample}",
        nodes = "8"
    container:
        "../sifs/virsorter_2.2.4--pyhdfd78af_1.sif"
    shell:
        "virsorter run -w {params.path} -i {input} -j {params.nodes} all"
        
rule vs2_RNA:
    input:
        "results/assembly/spades3.13/{sample}/contigs.fasta"
    output:
        fa = "results/vs2_RNA/{sample}/final-viral-combined.fa",
        score ="results/vs2_RNA/{sample}/final-viral-score.tsv",
        boundary = "results/vs2_RNA/{sample}/final-viral-boundary.tsv"
    params:
        path = "results/vs2_RNA/{sample}",
        nodes = "8"
    container:
        "../sifs/virsorter_2.2.4--pyhdfd78af_1.sif"
    shell:
        "virsorter run -w {params.path} -i {input} --include-groups RNA -j {params.nodes} all" # include RNA virus group
    
    # snakemake -c 48 -s workflow/rules/virsorter2.smk --use-singularity results/vs2_RNA/16_4_S4/final-viral-score.tsv -p


rule vs2_db:
    output:
        db = "resources/vs2_db",
    container:
        "../sifs/virsorter_2.2.4--pyhdfd78af_1.sif"
    shell:
        "virsorter config --init-source --db-dir={output.db}"
        
# #to have wildcards in the input of a rule but not in the output of the rule
# def table_inputs(folder, name, wildcards):
#     files=expand("%s{sample}%s" % (folder, name), sample=SAMPLE)
#     return files

# rule merge_vs2:
#     input:
#         i=table_inputs(folder="virsorter2/", name="/final-viral-score.tsv", wildcards=SAMPLE)
#     params:
#         s=SAMPLE
#     output:
#         o="virsorter2/vs2_merged_file.txt"
#     script:
#         "Scripts/merge_virsorter2.py"


