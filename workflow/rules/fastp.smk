# SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()

# SAMPLES, = glob_wildcards("rawdata/{sample}_L001_R1_001.fastq.gz")

# rule all:
#     input: expand("results/trimmed/{sample}.html", sample=SAMPLES)

# import pandas as pd
# wildcard_constraints:
#     dataset="\d+"

# samples: config/samples-template.tsv

# samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
# SAMPLES = samples_df["sample_name"].tolist()
# print(SAMPLES)

SAMPLES = "16_1_S1 16_2_S2 16_3_S3 16_4_S4 16_5_S5 22_1_S6 22_2_S7 22_3_S8 22_4_S9 22_5_S10 23_1_S11 23_2_S12 23_3_S13 23_4_S14 23_5_S15 24_1_S16 24_2_S17 24_3_S18 24_4_S19 24_5_S20 29_1_S21 29_2_S22 29_3_S23 29_4_S24 29_5_S25 Blank_8_15_23_S26".split()
print(SAMPLES)



# rule qfilterVis:
#     input: 
#         pdf = "results/trimmed/qfilterVis.pdf"
#     output: 
#         text = "results/trimmed/qfilterVis.stats",
#         plot = "results/trimmed/qfilterVis.pdf"
#     script:
#         """
#         export OMP_NUM_THREADS=48
#         module load r/4.3.2/gcc-9.3.0
#         # cd reports/fastp

#         echo -e "\nGenerating quality filtering results file qfilter.stats: ... "
#         for folder in report/fastp/;do
#             for file in $folder*.json;do
#                 ID=$(echo $file|sed 's|/.*$||g')
#                 readsBF=$(head -n 25 $file|grep total_reads|cut -d ':' -f2|sed 's/,//g'|head -n 1)
#                 readsAF=$(head -n 25 $file|grep total_reads|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
#                 basesBF=$(head -n 25 $file|grep total_bases|cut -d ':' -f2|sed 's/,//g'|head -n 1)
#                 basesAF=$(head -n 25 $file|grep total_bases|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
#                 q20BF=$(head -n 25 $file|grep q20_rate|cut -d ':' -f2|sed 's/,//g'|head -n 1)
#                 q20AF=$(head -n 25 $file|grep q20_rate|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
#                 q30BF=$(head -n 25 $file|grep q30_rate|cut -d ':' -f2|sed 's/,//g'|head -n 1)
#                 q30AF=$(head -n 25 $file|grep q30_rate|cut -d ':' -f2|sed 's/,//g'|tail -n 1)
#                 percent=$(awk -v RBF="$readsBF" -v RAF="$readsAF" 'BEGIN{{print RAF/RBF}}' )
#                 echo "$ID $readsBF $readsAF $basesBF $basesAF $percent $q20BF $q20AF $q30BF $q30AF" >> qfilter.stats
#                 echo "Sample $ID retained $percent * 100 % of reads ... "
#             done
#         done

#         echo "Done summarizing quality filtering results ... \nMoving to /stats/ folder and running plotting script ... "
#         # mv qfilter.stats {config[path][root]}/{config[folder][stats]}
#         # cd {config[path][root]}/{config[folder][stats]}

#         Rscript workflow/scripts/qfilterVis.R
#         echo "Done. "
#         rm Rplots.pdf
#         """

rule multiqc:
    input:
        fastp= expand("reports/trimmed/{sample}.json", sample = SAMPLES),
    output:
        mqc_out = directory('results/multiqc_out'),
        mqc_in  = directory('results/multiqc_in'),
    container:
        "docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
    shell:
        """mkdir {output.mqc_in}
           ln -snr -t {output.mqc_in} {input}
           multiqc {output.mqc_in} -o {output.mqc_out}
        """


rule fastp:
    output: 
        fq1 ="results/trimmed/{sample}.R1.fastq.gz",
        fq2 ="results/trimmed/{sample}.R2.fastq.gz",
        html = "reports/trimmed/{sample}.html",
        json="reports/trimmed/{sample}.json"
    input:  
        R1 = ["rawdata/{sample}_L001_R1_001.fastq.gz"],
        R2 = ["rawdata/{sample}_L001_R2_001.fastq.gz"],
    threads: 4
    container: 
        "../sifs/fastp_0.23.3--h5f740d0_0.sif"
    shell:
        "fastp -i {input.R1} -I {input.R2}  -o {output.fq1} -O {output.fq2} --json {output.json} --html {output.html}"
    
# rule fastp:
#     input:
#         sample=["rawdata/{sample}_L001_R1_001.fastq.gz", "rawdata/{sample}_L001_R2_001.fastq.gz"]
#     output:
#         trimmed1="results/trimmed/{sample}.R1.fastq.gz", 
#         trimmed2="results/trimmed/{sample}.R2.fastq.gz",
#         html="reports/trimmed/{sample}.html",
#         json="reports/trimmed/{sample}.json"
#     # conda:
#     #     "envs/fastp.yml"
#     threads: 4
#     wrapper:
#         "v3.0.2/bio/fastp"