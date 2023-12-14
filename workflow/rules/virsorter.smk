# fetch testing data
# wget -O test.fa https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa
# # run classification with 4 threads (-j) and test-out as output diretory (-w)
# virsorter run -w test.out -i test.fa --min-length 1500 -j 4 all
# ls test.out

# singularity exec [options] <container> <command>
singularity exec -B /project/awlab/wuhuiyun/pav_mgs_2023 \
./workflow/envs/virsorter_2.2.4--pyhdfd78af_1.sif \
virsorter run -w /results/virsorter/16_5_S5.out -i /results/assembly/16_5_S5.contigs.fasta --min-length 1500 -j 4 all

samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)

rule fastp_html:
    input: expand("reports/fastp/{sample}.html", sample=SAMPLES)

rule virsorter2:
    input:
        "Assemblies/{sample}_contigs.fasta"
    output:
        out1="VirSorter2/{sample}/final-viral-score.tsv",
        out2="VirSorter2/{sample}/final-viral-boundary.tsv"
    params:
        path="results/virSorter2/",
        nodes="16"
    container:
        "envs/virsorter_2.2.4--pyhdfd78af_1.sif"
    shell:
        "virsorter run --prep-for-dramv -w {params.path} -i {input} -j {params.nodes} all"