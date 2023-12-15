# fetch testing data
# wget -O test.fa https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa
# # run classification with 4 threads (-j) and test-out as output diretory (-w)
# virsorter run -w test.out -i test.fa --min-length 1500 -j 4 all
# ls test.out

# singularity exec [options] <container> <command>
singularity run -B /project/awlab/wuhuiyun/pav_mgs_2023 \
./resources/sifs/virsorter2.sif \ 
virsorter --version 
# run setup
# rm -rf db
virsorter run -w test.out -i test.fa --min-length 1500 -j 4 all

singularity run -B /project/  ./awlab/wuhuiyun/pav_mgs_2023/.test/virsorter-test/resources/sifs/virsorter2.sif run -w test.out -i test.fa --min-length 1500 -j 4 all

singularity run -B /project/awlab/wuhuiyun/pav_mgs_2023 \
./resources/sifs/virsorter2.sif \
virsorter run -w ./results/virsorter/16_5_S5.out -i ./results/assembly/16_5_S5.contigs.fasta --min-length 300 -j 4 all

singularity run -B /project virsorter2.sif --version
singularity run -B /project virsorter2.sif setup -d db -j 4
singularity run -B /project virsorter2.sif run -w test.out -i test.fa --min-length 1500 -j 4 all
singularity run -B /project virsorter2.sif run -w test2.out -i test.fa --min-length 1500 -j 4 all

singularity run -B /project virsorter2.sif run -w test2.out -i test.fa --min-length 1500 -j 4 all

#use another virsorter sif
singularity run -B /project virsorter2.2.4.sif --version

singularity run -B /project virsorter2.2.4.sif  -w test2.out -i test.fa --min-length 1500 -j 4 all
virsorter  run -w test2.out -i test.fa --min-length 1500 -j 4 all # used singularity shell
# working diretory= pav_mgs_2023
singularity run -B /project resources/sifs/virsorter2.2.4.sif virsorter run -w results/virsorter/16_1_S1 -i results/assembly/test/16_1_S1/contigs.fasta --min-length 1500 -j 4 all
singularity run -B /project resources/sifs/virsorter2.2.4.sif virsorter run -w results/virsorter/16_5_S5 -i results/assembly/test/16_5_S5/contigs.fasta --min-length 1500 -j 4 all
# singularity shell
singularity shell -B /project resources/sifs/virsorter2.2.4.sif 

singularity run -B /project resources/sifs/virsorter2.2.4.sif virsorter run --prep-for-dramv -w results/virsorter/16_4_S4 -i results/assembly/test/16_4_S4/contigs.fasta --min-length 1500 -j 4 all"











samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)

rule all_virsorter2:
    input:
        all="VirSorter2/vs2_merged_file.txt"

# download virsorter database in current directory
rule virsorter_db:
    output:
        "db-vs2/Done_all_setup"
    conda:
        "Envs/virsorter2.yaml"
    shell:
        "virsorter setup -d db-vs2 -j 4"

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