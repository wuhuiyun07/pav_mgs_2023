import pandas as pd
wildcard_constraints:
    dataset="\d+"

# samples: config/samples-template.tsv

samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)

rule all:
    input:
        trimmed=["results/trimmed/{sample}.R1.fastq.gz", "results/trimmed/{sample}.R2.fastq.gz"],

# rule OMP:
#     shell:
#         "export OMP_NUM_THREADS=48"

rule run_metaspades:
    input:
        reads=["results/trimmed/{sample}.R1.fastq.gz", "results/trimmed/{sample}.R2.fastq.gz"],
    output:
        contigs="results/assembly/test/{sample}.contigs.fasta",
        scaffolds="results/assembly/test/{sample}.scaffolds.fasta",
        dir=directory("results/assembly/test/{sample}_intermediate_files"),
    params:
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/{sample}.spades2.log",
    threads: 48,
    resources:
        mem_mb=180000,
        cpus_per_task=48,
        time=60 * 24,
    # container:
    #     "docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0",
    conda:
        "../envs/spades.yml"   
    # shell:
    #     # "singularity exec docker://quay.io/biocontainers/spades:3.15.5--h95f258a_0 "
    #     "spades.py --meta "
    #     "--threads {threads} "
    #     "--memory {resources.mem_mb} "  # Use resources.mem_mb instead of resources.mem
    #     # "{input.reads} "
    #     "--pe1-1 {input.reads[0]} "  # Specify the first mate pair
    #     "--pe1-2 {input.reads[1]} "  # Specify the second mate pair
    #     "-o {output} "        
    #     "-k {params.k} "       
    #     "{params.extra} "
    #     "> {log} 2>&1"
    script:
        "../scripts/spades_script.py"
