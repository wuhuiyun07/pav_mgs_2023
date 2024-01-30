    
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
    wrapper:
        "v3.3.5-42-g895739f/bio/diamond/makedb"


from snakemake.shell import shell


extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

rule diamond_makedb_2:
    input:
        fname = "resources/ncbi_db/genome/{reference}.fna"
    output:
        fname = "resources/ncbi_db/genome/{reference}.dmnd"
    log:
        "logs/diamond_makedb/{reference}.log"
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
        tsv = "results/diamond/{sample}.tsv"
    input:    
        fa = "results/vs2_{cond}/{sample}/final-viral-combined.fa"
    container:
        "../sifs/diamond_lastest.sif"
    shell:
        "diamond blastx "

    