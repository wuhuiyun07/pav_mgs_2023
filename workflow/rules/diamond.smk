SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()
print(SAMPLES)

NA = ["vs2_DNA", "vs2_RNA"]
print(NA)


rule all_diamond:
    input: 
        expand("results/diamond_blastp/{vs2_na}/{sample}.tsv", vs2_na = NA, sample=SAMPLES)


rule diamond_vs2_RNA:
    output:
        tsv = "results/diamond_blastp/vs_RNA/{sample}.tsv"
    input:    
        fa = "results/vs2_RNA/{sample}/final-viral-combined.fa",
        db = "resources/ncbi_db/protein/viral.1.protein.dmnd"
    container:
        "../sifs/diamond_latest.sif"
    threads: 8
    log:
        "reports/diamond_blastp/vs_RNA/{sample}.log",
    params:
        fmt = "6 qseqid sseqid pident length mismatch evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle",  # Additional arguments
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
        

rule ncbi_db:
    output: 
        protein ="resources/ncbi_db/protein/viral.1.protein.faa"
    shell:
        "wget -P resources/ncbi_db/protein/https://:viral/viral.1.protein.faa.gz"
        "wget -P resources/ncbi_db/protein/ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
        "gunzip {output.protein}"


rule diamond_makedb:
    input:
        refseq = "resources/ncbi_db/protein/{reference}.faa",
        nodes = "resources/ncbi_db/nodes.dmp",
        names = "resources/ncbi_db/names.dmp",
        map = "resources/ncbi_db/protein/prot.accession2taxid.gz"
        
    output:
        "resources/ncbi_db/protein/{reference}.dmnd",
    log:
        "logs/diamond_makedb/{reference}.log"
    threads: 8
    container:
        "../sifs/diamond_latest.sif"
    shell:
        "diamond makedb"
        " --threads {threads}"
        " --in {input.refseq}"
        " --db {output}"
        " --taxonmap {input.map}"
        " --taxonnodes {input.nodes}"
        " --taxonnames {input.names}"
        
