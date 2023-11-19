rule fastp:
    output: 
    "report/{sample}.html"

rule fastp_pe:
# https://github.com/OpenGene/fastp/
    input:
        sample=["test_fastp/{sample}_L001_R1_001.fastq.gz", "test_fastp/{sample}_L001_R2_001.fastq.gz"]
    output:
        trimmed=["test_fastp/trimmed_{sample}_R1.fastq.gz", "test_fastp/trimmed_{sample}_R2.fastq.gz"],
        html="report/{sample}.html",
        json="report/{sample}.json"
    log:
        "reports/{sample}.log"
    threads: 2,
    # conda:
    #     "workflow/envs/fastp.yml"
    # shell:
    #     "fastp -i test_fastp/{sample}_L001_R1_001.fastq.gz -I test_fastp/{sample}_L001_R2_001.fastq.gz -o test_fastp/trimmed_{sample}_R1.fastq.gz -O test_fastp/trimmed_{sample}_R2.fastq.gz"
    wrapper:
        "v2.13.0/bio/fastp"

# rule all:
#     input: 
#     "report/29_2_S22.html"
#     "report/29_5_S25.html"
#     "repot/16_1_S1.html" 
#     "repot/16_2_S2.html" 
#     "repot/16_3_S3.html" 
#     "repot/16_4_S4.html" 
#     "repot/16_5_S5.html" 
#     "repot/22_1_S6.html" 
#     "repot/22_2_S7.html" 
#     "repot/22_3_S8.html" 
#     "repot/22_4_S9.html" 
#     "repot/22_5_S10.html" 
#     "repot/23_1_S11.html" 
#     "repot/23_2_S12.html" 
#     "repot/23_3_S13.html" 
#     "repot/23_4_S14.html" 
#     "repot/23_5_S15.html" 
#     "repot/24_1_S16.html" 
#     "repot/24_2_S17.html" 
#     "repot/24_3_S18.html" 
#     "repot/24_4_S19.html" 
#     "repot/24_5_S20.html" 
#     "repot/29_1_S21.html" 
#     "repot/29_2_S22.html" 
#     "repot/29_3_S23.html" 
#     "repot/29_4_S24.html" 
#     "repot/29_5_S25.html"
#     "repot/Blank_8_15_23_S26.html" 

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

