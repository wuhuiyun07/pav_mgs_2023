DATA_READS = f'{config["path"]["root"]}/{config["folder"]["data"]}/{{IDs}}/{{IDs}}.fastq.gz'

rule qfilter: 
    input:
        READS = DATA_READS
    output:
        f'{config["path"]["root"]}/{config["folder"]["qfiltered"]}/{{IDs}}/{{IDs}}.fastq.gz', 
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

