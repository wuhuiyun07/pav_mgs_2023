#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

def read_data(vs2_path, checkV_path, diamond_path, output_path):
    # Read data files
    vs2 = pd.read_csv(vs2_path, sep='\t')
    checkV = pd.read_csv(checkV_path, sep='\t')
    diamond = pd.read_csv(diamond_path, sep='\t', skiprows=3, header=None)
    
    # Rename column names for diamond data
    diamond.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "evalue", 
                       "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", 
                       "sphylums", "stitle"]
    
    # Filter vs2 data
    vs2_screened = vs2[(vs2['max_score'] > 0.5) & (vs2['length'] > 1000)]
    vs2_screened[['contig_id', 'gene']] = vs2_screened['seqname'].str.split(r'\|\|', expand=True)
    
    # Filter checkV data
    checkV_screened = checkV[(checkV['checkv_quality'].isin(["Low-quality", "Medium-quality", "High-quality"])) & 
                             (checkV['contig_length'] > 1000)]
    
    # Join vs2 and checkV data
    contigs_for_diamond = pd.merge(vs2_screened, checkV_screened, on="contig_id", how="inner")
    
    # Filter and sort diamond data
    diamond_combined = diamond[(diamond[2] > 30) & (diamond[6] > 50) & (diamond[3] > 30)]
    diamond_combined.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "evalue", 
                                "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", 
                                "sphylums", "stitle"]
    
    # Perform left join
    screen_result = pd.merge(contigs_for_diamond, diamond_combined, on="contig_id", how="left")
    screen_result.drop_duplicates(subset='contig_id', keep="first", inplace=True)
    screen_result.sort_values(by=['checkv_quality', 'max_score_group', 'skingdoms', 'bitscore'], ascending=[True, True, False, False], inplace=True)
    
    # Save result to file
    screen_result.to_csv(output_path, index=False)

# Example usage:
vs2_path = "results/vs2/24_4_S19.vs2.final-viral-score.tsv"
checkV_path = "results/checkV/24_4_S19.checkv.quality_summary.tsv"
diamond_path = "results/diamond/24_4_S19.tsv.gz"
output_path = "results/screen_24_4_S19.csv"

read_data(vs2_path, checkV_path, diamond_path, output_path)
