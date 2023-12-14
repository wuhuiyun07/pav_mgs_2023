rule run_metaspades:
    input:
        reads=["results/trimmed/24_4_S19.R1.fastq.gz", "results/trimmed/24_4_S1.R2.fastq.gz"]
    output:
        contigs="results/assembly/test/24_4_S19.contigs.fasta",
    params:
        k="auto",
        extra="--only-assembler",
    log:
        "reports/assembly/24_4_S19.spades2.log",
    shell:
        """
        spades.py --version
        export OMP_NUM_THREADS=48
        spades.py --meta \\
            -1 results/trimmed/24_4_S19.R1.fastq.gz \\
            -2 results/trimmed/24_4_S19.R2.fastq.gz \\
            -o results/assembly/test/24_4_S19/
        """
# bash:
        
#         spades.py --version
#         export OMP_NUM_THREADS=48
#         spades.py --meta \
#             -1 results/trimmed/24_4_S19.R1.fastq.gz \
#             -2 results/trimmed/24_4_S19.R2.fastq.gz \
#             -o results/assembly/test/24_4_S19/
        