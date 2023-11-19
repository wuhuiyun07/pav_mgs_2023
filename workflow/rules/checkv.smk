rule all_checkV:
    input:
        merged="checkV/merged_checkV_allcontigs.tsv"

rule checkV_contamination:
    input:
        fasta="Assemblies/{sample}_3000.fasta"
    params:
        threads="4",
        directory="/nfs/turbo/cee-kwigg/hegartyb/SnakemakeAssemblies3000/checkV/{sample}",
        database="/home/hegartyb/checkV-db-v0.6/"
    container:
        "docker://quay.io/biocontainers/checkV:0.9.0--pyhdfd78af_0"
    output:
        out="checkV/{sample}/contamination.tsv"
    shell:
        """
        checkV contamination {input.fasta} {params.directory} -t {params.threads} -d {params.database}
        """

rule checkV_completeness:
    input:
        fasta="ViralSeqs/{sample}_3000_viralonly.fa",
        contamination="checkV/{sample}/contamination.tsv"
    output:
        out="checkV/{sample}/completeness.tsv"
    params:
        threads="4",
        directory="/nfs/turbo/cee-kwigg/hegartyb/SnakemakeAssemblies3000/checkV/{sample}",
        database="/home/hegartyb/checkV-db-v0.6/"
    container:
        "docker://quay.io/biocontainers/checkV:0.9.0--pyhdfd78af_0"
    shell:
        """
        checkV completeness {input.fasta} {params.directory} -t {params.threads} -d {params.database}
        """

rule checkV_genomes:
    input:
        fasta="ViralSeqs/{sample}.fa",
        complete="checkV/{sample}/completeness.tsv"
    params:
        threads="16",
        directory="/nfs/turbo/cee-kwigg/hegartyb/SnakemakeAssemblies3000/checkV/{sample}",
        database="/home/hegartyb/checkV-db-v0.6/"
    container:
        "docker://quay.io/biocontainers/checkV:0.9.0--pyhdfd78af_0"
    output:
        out="checkV/{sample}/complete_genomes.tsv"
    shell:
        """
        checkV complete_genomes {input.fasta} {params.directory}
        """

rule checkV_quality:
    input:
        fasta="ViralSeqs/{sample}_3000_viralonly.fa",
        genome="checkV/{sample}/complete_genomes.tsv"
    params:
        threads="16",
        directory="/nfs/turbo/cee-kwigg/hegartyb/SnakemakeAssemblies3000/checkV/{sample}",
        database="/home/hegartyb/checkV-db-v0.6/"
    container:
        "docker://quay.io/biocontainers/checkV:0.9.0--pyhdfd78af_0"
    output:
        out="checkV/{sample}/quality_summary.tsv"
    shell:
        """
        checkV quality_summary {input.fasta} {params.directory}
        """

#to have wildcards in the input of a rule but not in the output of the rule
def table_inputs(folder, name, wildcards):
    files=expand("%s{sample}%s" % (folder, name), sample=SAMPLE)
    return files

rule merge_checkV:
    input:
        i=table_inputs(folder="checkV/", name="/quality_summary.tsv", wildcards=SAMPLE)
    output:
        merge="checkV/merged_checkV_allcontigs.tsv"
    script:
        "Scripts/merge_checkV.py"


