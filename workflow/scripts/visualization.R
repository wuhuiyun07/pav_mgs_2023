library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)

# setwd("/Users/huiyunwu/Desktop/Virus_particle/pav_mgs_2023")
getwd()

test_dat <- read_csv(snakemake@input[["test"]])

diamond<-read_tsv("results/diamond/16_4_S4.tsv", skip=3, col_names=FALSE, show_col_types = FALSE)
# read_tsv("results/diamond/16_4_S4.tsv", skip=3, col_names=FALSE, show_col_types = FALSE)
colnames(diamond)<- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")
diamond

vs2 <- read_tsv(file = "results/vs2/16_4_S4.vs2.final-viral-score.tsv")
checkV <- read_tsv(file = "results/checkV/16_4_S4.checkv.quality_summary.tsv")


checkV_screened <- checkV %>% 
  filter(checkv_quality == "Low-quality" |
         checkv_quality == "Medium-quality"|
         checkv_quality == "High-quality" ) %>%
  filter(contig_length > 1000)


vs2_screened <- vs2 %>% 
  filter(max_score > 0.5) %>% 
  filter(length > 1000)

vs2_screened_sep <- separate(vs2_screened, seqname, into = c("contig_id", "gene"), sep = "\\|\\|" )

contigs_for_diamond <- full_join(vs2_screened_sep, checkV_screened, by = c ("contig_id" = "contig_id"))


diamond_combined <- diamond %>%
  arrange(., desc(bitscore)) %>% 
  arrange(., evalue) %>% 
  arrange(., desc(length)) %>% 
  arrange(., desc(pident)) %>% 
  filter(pident > 30 &
         bitscore >50 &
         length > 30)

diamond_combined_sep <- separate(diamond_combined, qseqid, into = c("contig_id", "gene"), sep = "\\|\\|" )


annotation<-left_join(contigs_for_diamond, diamond_combined_sep, 
                      by = c ("contig_id" = "contig_id")) %>% 
            distinct(contig_id, .keep_all = TRUE)

write_csv(annotation, file = "results/annotation.csv")


# %>% write_csv(snakemake@output[["csv"]])

