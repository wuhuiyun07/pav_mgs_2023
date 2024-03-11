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

vs2_paths <- list.files("results/vs2", full.names = TRUE)
checkV_paths <- list.files("results/checkV", full.names = TRUE)
diamond_paths <- list.files("results/diamond", full.names = TRUE)
output_paths <- paste0("results/annotation/result_", 1:length(vs2_paths), ".csv")

# Iterate over input and output paths and call read_data function
for (i in 1:length(vs2_paths)) {
  read_data(vs2_paths[i], checkV_paths[i], diamond_paths[i], output_paths[i])
}
# example:
# read_data("results/vs2/24_4_S19.vs2.final-viral-score.tsv","results/checkV/24_4_S19.checkv.quality_summary.tsv","results/diamond/24_4_S19.tsv.gz", "results/24_4_S19_screened.csv")

# vs2_path <- snakemake@input$vs2_path
# checkV_path <- snakemake@input$checkV_path
# diamond_path <- snakemake@input$diamond_path
# output_path <- snakemake@output[1]


# read_data(vs2_path, checkV_path, diamond_path, output_path)


# # trail 8 
# SAMPLES <- c("16_1_S1", "16_2_S2", "16_3_S3", "16_4_S4", "16_5_S5", "22_1_S6", "22_2_S7", "22_3_S8", "22_4_S9", "22_5_S10", "23_1_S11", "23_2_S12", "23_3_S13", "23_4_S14", "23_5_S15", "24_1_S16", "24_2_S17", "24_3_S18", "24_4_S19", "24_5_S20", "29_1_S21", "29_2_S22", "29_3_S23", "29_4_S24", "29_5_S25" )


# read_data <- function() {
#   vs2_files <- list.files("results/vs2", full.names = TRUE)
#   checkV_files <- list.files("results/checkV", full.names = TRUE)
#   diamond_files <- list.files("results/diamond", full.names = TRUE)
  
#   for (i in 1:length(vs2_files)) {
#     x <- vs2_files[i]
#     y <- checkV_files[i]
#     z <- diamond_files[i]
    
#     vs2 <- read_tsv(x, col_names = TRUE)
#     checkV <- read_tsv(y, col_names = TRUE)
#     diamond <- read_tsv(z, skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)
    
#     # Rename column names for diamond data
#     colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")
    
#     # Filter vs2 data
#     vs2_screened <- vs2 %>%
#       filter(max_score > 0.5,
#              length > 1000) %>%
#       separate(seqname, into = c("contig_id", "gene"), sep = "\\|\\|")
    
#     # Filter checkV data
#     checkV_screened <- checkV %>%
#       filter(checkv_quality %in% c("Low-quality", "Medium-quality", "High-quality"),
#              contig_length > 1000)
    
#     # Join vs2 and checkV data
#     contigs_for_diamond <- full_join(vs2_screened, checkV_screened, by = "contig_id")
    
#     # Filter and sort diamond data
#     diamond_combined <- diamond %>%
#       filter(pident > 30,
#              bitscore > 50,
#              length > 30) %>%
#       arrange(desc(bitscore), evalue, desc(length), desc(pident)) %>%
#       separate(qseqid, into = c("contig_id", "gene"), sep = "\\|\\|")
    
#     # Perform left join
#     screen_result <- left_join(contigs_for_diamond, diamond_combined, by = "contig_id") %>%
#       distinct(contig_id, .keep_all = TRUE) %>% 
#       arrange(checkv_quality, max_score_group, desc(skingdoms), desc(bitscore)) 
    
#     # Save results to file
#     filename <- paste0("result_", i, ".csv")
#     write_csv(screen_result, filename)
#   }
# }