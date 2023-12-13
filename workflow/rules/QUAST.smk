rule quast:
    input:
        fasta="genome.fasta",
        ref="genome.fasta",
        #gff="annotations.gff",
        #pe1="reads_R1.fastq",
        #pe2="reads_R2.fastq",
        #pe12="reads.fastq",
        #sv_bedpe="sv.bed",
    output:
        multiext("{sample}/report.", "html", "tex", "txt", "pdf", "tsv"),
        multiext("{sample}/transposed_report.", "tex", "txt", "tsv"),
        multiext(
            "{sample}/basic_stats/",
            "cumulative_plot.pdf",
            "GC_content_plot.pdf",
            "gc.icarus.txt",
            "genome_GC_content_plot.pdf",
            "NGx_plot.pdf",
            "Nx_plot.pdf",
        ),
        multiext(
            "{sample}/contigs_reports/",
            "all_alignments_genome.tsv",
            "contigs_report_genome.mis_contigs.info",
            "contigs_report_genome.stderr",
            "contigs_report_genome.stdout",
        ),
        "{sample}/contigs_reports/minimap_output/genome.coords_tmp",
        "{sample}/icarus.html",
        "{sample}/icarus_viewers/contig_size_viewer.html",
        "{sample}/quast.log",
    log:
        "logs/{sample}.quast.log",
    params:
        extra="--min-contig 5 --min-identity 95.0",
    wrapper:
        "v3.1.0/bio/quast"