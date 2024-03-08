#!/usr/bin/env python

import pandas as pd
import os
import shutil
from snakemake.shell import shell

# Read input files
vs2_file = snakemake.input["vs2_file"]
checkV_file = snakemake.input["checkV_file"]
diamond_file = snakemake.input["diamond_file"]

vs2 = pd.read_csv(vs2_file, sep="\t")
checkV = pd.read_csv(checkV_file, sep="\t")
diamond = pd.read_csv(diamond_file, sep="\t", skiprows=3, header=None)

# Rename column names for diamond data
diamond.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle"]

# Filter checkV data
checkV_screened = checkV[(checkV["checkv_quality"].isin(["Low-quality", "Medium-quality", "High-quality"])) & (checkV["contig_length"] > 1000)]

# Filter vs2 data
vs2_screened = vs2[(vs2["max_score"] > 0.5) & (vs2["length"] > 1000)]

# Separate seqname column in vs2_screened
vs2_screened[["contig_id", "gene"]] = vs2_screened["seqname"].str.split("||", expand=True)

# Join vs2 and checkV data
contigs_for_diamond = pd.merge(vs2_screened, checkV_screened, on="contig_id", how="inner")

# Filter and sort diamond data
diamond_combined = diamond[(diamond["pident"] > 30) & (diamond["bitscore"] > 50) & (diamond["length"] > 30)]
diamond_combined = diamond_combined.sort_values(by=["bitscore", "evalue", "length", "pident"], ascending=[False, True, False, False])

# Separate qseqid column in diamond_combined
diamond_combined[["contig_id", "gene"]] = diamond_combined["qseqid"].str.split("||", expand=True)

# Perform left join and write output
annotation = pd.merge(contigs_for_diamond, diamond_combined, on="contig_id", how="left").drop_duplicates(subset="contig_id", keep="first")
annotation.to_csv(snakemake.output["annotation"], sep="\t", index=False)
