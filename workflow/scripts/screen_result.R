#!/usr/bin/env Rscript
# trail #7
library(readr)
library(dplyr)
library(tidyverse)

read_data <- function(vs2_path, checkV_path, diamond_path, output_path) {
  vs2 <- read_tsv(vs2_path, col_names = TRUE)
  checkV <- read_tsv(checkV_path, col_names = TRUE)
  diamond <- read_tsv(diamond_path, skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)
  
  # Rename column names for diamond data
  colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")
  
  # Filter vs2 data
  vs2_screened <- vs2 %>%
    filter(max_score > 0.5,
           length > 1000) %>%
    separate(seqname, into = c("contig_id", "gene"), sep = "\\|\\|")
  
  # Filter checkV data
  checkV_screened <- checkV %>%
    filter(checkv_quality %in% c("Low-quality", "Medium-quality", "High-quality"),
           contig_length > 1000)
  
  # Join vs2 and checkV data
  contigs_for_diamond <- full_join(vs2_screened, checkV_screened, by = "contig_id")
  
  # Filter and sort diamond data
  diamond_combined <- diamond %>%
    filter(pident > 30,
           bitscore > 50,
           length > 30) %>%
    arrange(desc(bitscore), evalue, desc(length), desc(pident)) %>%
    separate(qseqid, into = c("contig_id", "gene"), sep = "\\|\\|")
  
  # Perform left join
  screen_result <- left_join(contigs_for_diamond, diamond_combined, by = "contig_id") %>%
    distinct(contig_id, .keep_all = TRUE) %>% 
    arrange(checkv_quality, max_score_group, desc(skingdoms), desc(bitscore)) 
  
  # Save result to file
  write_csv(screen_result, output_path)
}

# example:
# read_data("results/vs2/24_4_S19.vs2.final-viral-score.tsv","results/checkV/24_4_S19.checkv.quality_summary.tsv","results/diamond/24_4_S19.tsv.gz", "results/24_4_S19_screened.csv")

vs2_path <- snakemake@input$vs2_path
checkV_path <- snakemake@input$checkV_path
diamond_path <- snakemake@input$diamond_path
output_path <- snakemake@output[1]


read_data(vs2_path, checkV_path, diamond_path, output_path)