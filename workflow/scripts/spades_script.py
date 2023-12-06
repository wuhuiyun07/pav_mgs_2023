import os
import shutil
from snakemake.shell import shell

# Infer output directory
if hasattr(snakemake.output, "dir"):
    output_dir = snakemake.output.dir
else:
    output_file = (
        snakemake.output.contigs
        if hasattr(snakemake.output, "contigs")
        else snakemake.output[0]
    )
    output_dir = os.path.split(output_file)[0]

# Parse parameters
extra = snakemake.params.get("extra", "")
kmers = snakemake.params.get("k", "'auto'")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

mem_gb = snakemake.resources.mem_mb // 1000 if hasattr(snakemake.resources, "mem_mb") else ""
memory_requirements = f"--memory {mem_gb}" if mem_gb else ""

# Check if params.txt exists
if not os.path.exists(os.path.join(output_dir, "params.txt")):
    # Parse short reads
    reads = snakemake.input.reads if hasattr(snakemake.input, "reads") else snakemake.input
    assert len(reads) > 1, "Metaspades requires a paired-end library with at least 2 fastq files."

    input_arg = f"--pe1-1 {reads[0]} --pe1-2 {reads[1]}"

    if len(reads) >= 3:
        input_arg += f" --pe1-m {reads[2]}"

        if len(reads) >= 4:
            input_arg += f" --pe1-s {reads[3]}"

    # Parse long reads
    for longread_name in ["pacbio", "nanopore"]:
        if hasattr(snakemake.input, longread_name):
            input_arg += f" --{longread_name} {getattr(snakemake.input, longread_name)}"

    # Execute MetaSPAdes
    shell(
        f"spades.py --meta "
        f"--threads {snakemake.threads} "
        f"{memory_requirements} "
        f"-o {output_dir} "
        f"-k {kmers} "
        f"{input_arg} "
        f"{extra} "
        f"> {snakemake.log[0]} 2>&1"
    )
else:
    # Restart Spades if params.txt exists
    shell(f"echo '\n\nRestart Spades...' >> {log[0]}")
    shell(f"rm -f {output_dir}/pipeline_state/stage_*_copy_files 2>> {log}")
    shell(
        f"spades.py --meta "
        f"--restart-from last "
        f"--threads {snakemake.threads} "
        f"{memory_requirements} "
        f"-o {output_dir} "
        f">> {snakemake.log[0]} 2>&1"
    )

    # Rename/move output files
    Output_key_mapping = {
        "contigs": "contigs.fasta",
        "scaffolds": "scaffolds.fasta",
        "graph": "assembly_graph_with_scaffolds.gfa",
    }

    for key, value in Output_key_mapping.items():
        if hasattr(snakemake.output, key):
            file_produced = os.path.join(output_dir, value)
            file_renamed = getattr(snakemake.output, key)
            if file_produced != file_renamed:
                shutil.move(file_produced, file_renamed)
