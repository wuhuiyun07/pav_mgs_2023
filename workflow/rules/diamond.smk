SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25".split()
print(SAMPLES)

rule all_diamond:
    input: expand("results/diamond_vs2/{sample}.diamond.tsv", sample = SAMPLES) 

rule chg_file_name:
    output:
        diamond = "results/diamond_vs2/{sample}.diamond.tsv",
        # vs2 = "results/vs2/{sample}.vs2.final-viral-score.tsv",
        # checkV = "results/checkV/{sample}.checkv.quality_summary.tsv"
    input:
        diamond = "results/diamond_vs2/{sample}.tsv",
        vs2 = "results/vs2/{sample}/final-viral-score.tsv",
        checkV = "results/checkV/{sample}/quality_summary.tsv"
    shell:
        r""" scp {input.diamond} {output.diamond} \"""
        # """ #scp {input.vs2} {output.vs2} \"""
        # """ scp {input.checkV} {output.checkV}  \"""


rule gunzip_diamond:
    output: 
        "results/diamond_vs2/{sample}.tsv"
    input:
        "results/diamond_vs2/{sample}.tsv.gz"
    shell:
        "gunzip {input} "


rule diamond_vs2:
    output:
        tsv = "results/diamond_vs2/{sample}.tsv.gz"
    input:    
        fa = "results/vs2/{sample}/final-viral-combined.fa",
        db = "resources/ncbi_db/protein/viral.1.protein.dmnd"
    container:
        "../sifs/diamond_latest.sif"
    threads: 8
    params:
        fmt = "6 qseqid sseqid pident length mismatch evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle",  # Additional arguments
    shell:
        r"""diamond blastx """
        """ --threads {threads}"""
        """ -q {input.fa} """
        """ -d {input.db} """
        """ -o {output.tsv} """
        """ --header """
        """ --outfmt {params.fmt} """
        """ --verbose """
        """ --compress 1 """


        
        

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
        
