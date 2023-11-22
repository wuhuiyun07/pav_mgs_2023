import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate


samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str}, encoding='latin-1')
    .set_index("sample_name", drop=False)
    .sort_index()
    
)




def get_final_output():
    final_output = expand(
	"results/trimmed/{samples}.html")	
    return final_output