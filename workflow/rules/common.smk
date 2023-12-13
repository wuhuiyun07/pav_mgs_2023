import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate


samples_df = pd.read_csv("config/samples-template.tsv", sep="\t")
SAMPLES = samples_df["sample_name"].tolist()
print(SAMPLES)


def get_final_output():
    final_output = expand("results/assembly/{sample}.fasta")	
    return final_output