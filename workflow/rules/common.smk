import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

# validate(config, schema="../config.schema.yaml")

sample = (pd.read_csv(config["sample"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

print(sample)




# def get_final_output():
#     final_output = expand(
# 	"results/trimmed/{samples}.tsv")	
#     return final_output