    
rule ncbi_db:
    output: 
        protein ="resources/ncbi_db/protein/viral.1.protein.faa"
    shell:
        "wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
        "gunzip viral.1.protein.faa.gz"
        "gunzip viral.1.1.genomic.fna.gz"
        


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
        " --in {input.fname}"
        " --db {output.fname}"



rule diamond_vs2:
    output:
        tsv = "results/diamond_blastp/vs2_RNA/{sample}.tsv"
    input:    
        fname_fasta = "results/vs2_RNA/{sample}/final-viral-combined.fa",
        fname_db = "resources/ncbi_db/protein/viral.1.protein.dmnd"
    # container:
    #     "../sifs/diamond_lastest.sif"
    threads: 8
    log:
        "logs/diamond_blastp/{sample}.RNA.log",
    params:
        extra= "--header --compress 1",  # Additional arguments
    wrapper:
        "v3.3.5-42-g895739f/bio/diamond/blastp"
    # shell:
    #     "diamond blastp "
    #     "--threads {threads}"
    #     "-q {input.fa} "
    #     "-d {input.db} "
    #     "-o {output.tsv}"
    #     "--verbose"
    #     "--log"

    