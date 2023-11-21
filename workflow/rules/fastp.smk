# RAWDATA = f'{config["path"]["root"]}/{config["folder"]["data"]}/{{IDs}}/{{IDs}}.fastq.gz'

# def get_input(wildcards):
#     unit = units.loc[wildcards.sample].loc[wildcards.unit]

#     if pd.isna(unit["fq1"]):
#         # SRA sample (always paired-end for now)
#         accession = unit["sra"]
#         return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

#     if unit["fq1"].endswith("gz"):
#         ending = ".gz"
#     else:
#         ending = ""

#     if pd.isna(unit["fq2"]):
#         # single end local sample
#         return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
#             S=unit.sample_name, U=unit.unit_name, E=ending
#         )
#     else:
#         # paired end local sample
#         return expand(
#             "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
#                 S=unit.sample_name, U=unit.unit_name, E=ending
#             ),
#             read=["fq1", "fq2"],
#         )
# rule all:
#     input:
#         expand("{sample}_{unit}.html",sample =config["samples"])

rule fastq: 
    input:
        get_input
    output:
        fastq1="results/trimmed/{sample}_{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}_{unit}_R2.fastq.gz",
        json="results/trimmed/{sample}_{unit}.json",
        html="results/trimmed/{sample}_{unit}.html",
    shell:
        """
        set +u;source activate {config[envs][metabagpipes]};set -u;

        mkdir -p $(dirname $(dirname {output}))
        mkdir -p $(dirname {output})

        fastp --thread {config[cores][fastp]} \
            -i {input} \
            -o {output} \
            -j $(dirname {output})/$(echo $(basename $(dirname {output}))).json \
            -h $(dirname {output})/$(echo $(basename $(dirname {output}))).html

        """


# rule fastp:
#     output: 
#     "report/{sample}.html"

# rule fastp_pe:
# # https://github.com/OpenGene/fastp/
#     input:
#         sample=["test_fastp/{sample}_L001_R1_001.fastq.gz", "test_fastp/{sample}_L001_R2_001.fastq.gz"]
#     output:
#         trimmed=["test_fastp/trimmed_{sample}_R1.fastq.gz", "test_fastp/trimmed_{sample}_R2.fastq.gz"],
#         html="report/{sample}.html",
#         json="report/{sample}.json"
#     log:
#         "reports/{sample}.log"
#     threads: 2,
#     # conda:
#     #     "workflow/envs/fastp.yml"
#     # shell:
#     #     "fastp -i test_fastp/{sample}_L001_R1_001.fastq.gz -I test_fastp/{sample}_L001_R2_001.fastq.gz -o test_fastp/trimmed_{sample}_R1.fastq.gz -O test_fastp/trimmed_{sample}_R2.fastq.gz"
#     wrapper:
#         "v2.13.0/bio/fastp"


# rule fastp_pe:
# # https://github.com/OpenGene/fastp/
#     input:
#         sample=["rawdata/{sample}_L001_R1_001.fastq.gz", "rawdata/{sample}_L001_R2_001.fastq.gz"]
#     output:
#         trimmed=["results/trimmed/{sample}_R1.fastq.gz", "results/trimmed/{sample}_R2.fastq.gz"],
#         html="report/fastp/{sample}.html",
#         json="report/fastp/{sample}.json"
#     log:
#         "reports/fastp/{sample}.log"
#     threads: 8
#     wrapper:
#         "v2.11.1/bio/fastp"

