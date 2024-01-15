REPLICATES = ["1", "2", "3"]
CONDITIONS = ["ref", "etoh60", "temp33"]
KMERS = ["19","21","25","27"]


rule all_max:
    input:
        expand("assem/{condition}_k{kmer}_max_contig.txt", condition = CONDITIONS,
							   kmer =KMERS)

rule cutadapt:
    output:
        read1 = "cutadapt/{sample}_1.fq",
        read2 = "cutadapt/{sample}_2.fq"
    input: 
        read1 = "reads/{sample}_1.fq",
        read2 = "reads/{sample}_2.fq"
    params:
        adapter = "AGATCGGAAGAGC"
    conda: "envs/cutadapt.yml"
    shell:
        """
        cutadapt -a {params.adapter} -A {params.adapter} \
	    -o {output.read1} -p {output.read2} \
           {input.read1} {input.read2}
        """

rule concatenate:
    output: 
        read1 = "concatenate/{condition}_1.fq",
        read2 = "concatenate/{condition}_2.fq"
    input:
        read1s = ["cutadapt/{condition}_1_1.fq", "cutadapt/{condition}_2_1.fq", "cutadapt/{condition}_3_1.fq"],
        read2s = ["cutadapt/{condition}_1_2.fq", "cutadapt/{condition}_2_2.fq", "cutadapt/{condition}_3_2.fq"],
    conda: "envs/cutadapt.yml"
    shell:
        r"""cat {input.read1s} > {output.read1} 
             cat {input.read2s} > {output.read2} 
         """

rule assembly:
    output: "assem/{sample}_k{kmer}_contigs.fa"
    input: 
        read1 = "concatenate/{sample}_1.fq",
        read2 = "concatenate/{sample}_2.fq"
    container:
         "docker://biocontainers/velvet:v1.2.10_cv3"
    params: 
        tmpdir = "velvet_tmp_{sample}_{kmer}"
    shell:
        r"""
        velveth {params.tmpdir} {wildcards.kmer} -shortPaired -fastq -separate {input.read1} {input.read2} 
        velvetg {params.tmpdir} 
        mv {params.tmpdir}/contigs.fa {output}
         """

rule max_contig:
    input: "assem/{assem}_contigs.fa"
    output: "assem/{assem}_max_contig.txt"
    container: "docker://quay.io/biocontainers/bbmap:38.96--h5c4e2a8_0"
    shell: 
        "stats.sh {input} | grep 'Max contig length:' > {output}"

