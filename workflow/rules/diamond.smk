SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()
print(SAMPLES)

NA = ["vs2_DNA", "vs2_RNA"]
print(NA)


rule all_diamond:
    input: 
        expand("results/diamond_blastp/{vs2_na}/{sample}.tsv", vs2_na = NA, sample=SAMPLES)


rule diamond_vs2:
    output:
        tsv = "results/diamond_blastp/{vs2_na}/{sample}.tsv"
    input:    
        fa = "results/{vs2_na}/{sample}/final-viral-combined.fa",
        db = "resources/ncbi_db/protein/viral.1.protein.dmnd"
    container:
        "../sifs/diamond_latest.sif"
    threads: 8
    log:
        "reports/diamond_blastp/{sample}.{vs_na}.log",
    params:
        fmt = "6 qseqid sseqid pident length mismatch evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle",  # Additional arguments
    # wrapper:
    #     "v3.3.5-42-g895739f/bio/diamond/blastp"
    shell:
        "diamond blastx "
        " --threads {threads}"
        " -q {input.fa} "
        " -d {input.db} "
        " -o {output.tsv}"
        " --header"
        " --outfmt {params.fmt}"
        " --verbose"
        " --compress 1"
        " --log {log}"


rule ncbi_db:
    output: 
        protein ="resources/ncbi_db/protein/viral.1.protein.faa"
    shell:
        "wget -P resources/ncbi_db/https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
        "gunzip {output.protein}"


rule diamond_makedb:
    input:
        fname = "resources/ncbi_db/protein/{reference}.faa",
    output:
        "resources/ncbi_db/protein/{reference}.dmnd",
    log:
        "logs/diamond_makedb/{reference}.log"
    params:
        extra=""
    threads: 8
    container:
        "../sifs/diamond_latest.sif"
    shell:
        "diamond makedb"
        " --threads {threads}"
        " --in {input}"
        " --taxonmap"
        " --taxonnames"
        " --db {output}"
