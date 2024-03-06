# #!/usr/bin/env Rscript

# library(tidyverse)
# library(dplyr)
# library(tibble)
# library(ggplot2)

# # setwd("/Users/huiyunwu/Desktop/Virus_particle/pav_mgs_2023")
# getwd()

# vs2_file <- input[1]
# checkV_file <- input[2]
# diamond_file <- input[3]

# output_file <- output[0]

# # test_dat <- read_csv(snakemake@input[["test"]])
# vs2 <- read_tsv(vs2_file)  #results/vs2/{sample}.vs2.final-viral-score.tsv
# checkV <- read_tsv(checkV_file) #results/checkV/{sample}.checkv.quality_summary.tsv
# diamond<-read_tsv(diamond_file, skip=3, col_names=FALSE, show_col_types = FALSE), results/diamond_vs2/{sample}.diamond.tsv
# # read_tsv("results/diamond/16_4_S4.tsv", skip=3, col_names=FALSE, show_col_types = FALSE)
# colnames(diamond)<- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")
# diamond


# checkV_screened <- checkV %>% 
#   filter(checkv_quality == "Low-quality" |
#          checkv_quality == "Medium-quality"|
#          checkv_quality == "High-quality" ) %>%
#   filter(contig_length > 1000)


# vs2_screened <- vs2 %>% 
#   filter(max_score > 0.5) %>% 
#   filter(length > 1000)

# vs2_screened_sep <- separate(vs2_screened, seqname, into = c("contig_id", "gene"), sep = "\\|\\|" )

# contigs_for_diamond <- full_join(vs2_screened_sep, checkV_screened, by = c ("contig_id" = "contig_id"))


# diamond_combined <- diamond %>%
#   arrange(., desc(bitscore)) %>% 
#   arrange(., evalue) %>% 
#   arrange(., desc(length)) %>% 
#   arrange(., desc(pident)) %>% 
#   filter(pident > 30 &
#          bitscore >50 &
#          length > 30)

# diamond_combined_sep <- separate(diamond_combined, qseqid, into = c("contig_id", "gene"), sep = "\\|\\|" )


# annotation<-left_join(contigs_for_diamond, diamond_combined_sep, 
#                       by = c ("contig_id" = "contig_id")) %>% 
#             distinct(contig_id, .keep_all = TRUE) %>% 
#             write_tsv(output_file)


# # %>% write_csv(snakemake@output[["csv"]])



#!/usr/bin/env Rscript

library(tidyverse)

# Set working directory if needed
# setwd("/Users/huiyunwu/Desktop/Virus_particle/pav_mgs_2023")

# Print working directory
getwd()

# Input and output files
vs2_file <- snakemake@input[["vs2_file"]]
checkV_file <- snakemake@input[["checkV_file"]]
diamond_file <- snakemake@input[["diamond_file"]]
output_file <- snakemake@output[["annotation"]]

# Read input files
vs2 <- read_tsv("results/vs2/{sample}.vs2.final-viral-score.tsv")
checkV <- read_tsv("results/checkV/{sample}.checkv.quality_summary.tsv")
diamond <- read_tsv("results/diamond_vs2/{sample}.diamond.tsv", skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)

# Rename column names for diamond data
colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")

# Filter checkV data
checkV_screened <- checkV %>%
  filter(checkv_quality %in% c("Low-quality", "Medium-quality", "High-quality")) %>%
  filter(contig_length > 1000)

# Filter vs2 data
vs2_screened <- vs2 %>%
  filter(max_score > 0.5) %>%
  filter(length > 1000) %>%
  separate(seqname, into = c("contig_id", "gene"), sep = "\\|\\|")

# Join vs2 and checkV data
contigs_for_diamond <- full_join(vs2_screened, checkV_screened, by = "contig_id")

# Filter and sort diamond data
diamond_combined <- diamond %>%
  filter(pident > 30 & bitscore > 50 & length > 30) %>%
  arrange(desc(bitscore), evalue, desc(length), desc(pident)) %>%
  separate(qseqid, into = c("contig_id", "gene"), sep = "\\|\\|")

# Perform left join and write output
annotation <- left_join(contigs_for_diamond, diamond_combined, by = "contig_id") %>%
  distinct(contig_id, .keep_all = TRUE) %>%
  write_tsv(output_file)



