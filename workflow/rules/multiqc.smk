# rule multiqc:
#     output:
#         mqc_out = directory('results/multiqc_out'),
#         mqc_in  = directory('results/multiqc_in'),
#     input:
#         salmon =   expand("salmon.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
#         kallisto = expand("kallisto.{cond}_{rep}", cond=CONDITIONS, rep=REPLICATES),
#         fastqc =   expand("reads.{cond}_{rep}_{end}_fastqc.zip", cond=CONDITIONS, rep=REPLICATES, end=["1","2"]),
#     container:
#         "docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
#     shell:
#         """mkdir {output.mqc_in}
#            ln -snr -t {output.mqc_in} {input}
#            multiqc {output.mqc_in} -o {output.mqc_out}
#         """

rule multiqc:
    input:
        fastp= "report/fastp/{sample}.html",
    output:
        mqc_out = directory('results/multiqc_out'),
        html=  "result/multiqc_report.html",
    container:
        "docker://quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0"
  shell:
        "multiqc {input.fastp} -o {output.mqc_out}"