rule run_metaspades:
    output:
        contigs="results/assembly/test2/{sample}.contigs.fasta",
        scaffolds="results/assembly/test2/{sample}.scaffolds.fasta",
        dir=directory("results/assembly/test2/{sample}_intermediate_files"),
    benchmark:
        "reports/assembly/{sample}.spades2.txt",
    params:
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/{sample}.spades2.log",
    threads: 48,
    resources:
        mem_mb=250000,
        cpus_per_task=48,
        time=60 * 24,
    script:
        "../scripts/spades_script.py"
