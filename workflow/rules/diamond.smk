    
rule ncbi_db:
    output: 
        protein ="resources/ncbi_db/protein/viral.1.protein.faa"
    shell:
        "wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
        "gunzip viral.1.protein.faa.gz"
        


rule diamond_makedb:
    input:
        fname = "resources/ncbi_db/protein/{reference}.faa",
    output:
        fname = "resources/ncbi_db/protein/{reference}.dmnd"
    log:
        "logs/diamond_makedb/{reference}.log"
    params:
        extra=""
    threads: 8
    wrapper:
        "v3.3.5-42-g895739f/bio/diamond/makedb"




rule diamond_vs2:
    output:
        tsv = "results/diamond/{sample}.tsv"
    input:    
        fa = "results/vs2_{cond}/{sample}/final-viral-combined.fa"
    container:
        "../sifs/diamond_lastest.sif"
    shell:
        "diamond blastx "

    