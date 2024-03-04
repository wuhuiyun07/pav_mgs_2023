library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)

setwd("/Users/huiyunwu/Desktop/Virus_particle/pav_mgs_2023")
getwd()

diamond<-read_tsv("results/diamond/16_4_S4.tsv", skip=3, col_names=FALSE, show_col_types = FALSE)
# read_tsv("results/diamond/16_4_S4.tsv", skip=3, col_names=FALSE, show_col_types = FALSE)
colnames(diamond)<- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")

mammal<-read_tsv("/Users/huiyunwu/Desktop/Mentees/Katie/ONR10623barcode03.txt", col_names = FALSE, show_col_types = FALSE)
# colnames(mammal)<- c("NA","qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")



diamond_screened <- diamond %>% 
  filter(pident > 30 &
         bitscore >50 &
         length > 30)


vs2 <- read_tsv(file = "results/vs2/22_5_S10.vs2.final-viral-score.tsv")
checkV <- read_tsv(file = "results/checkV/22_5_S10.checkv.quality_summary.tsv")

checkV_screened <- checkV %>% 
  filter(checkv_quality == "Low-quality" |
         checkv_quality == "Medium-quality"|
         checkv_quality == "High-quality" ) %>%
  filter(contig_length > 1500)

vs2_screened <- vs2 %>% 
  filter(max_score > 0.5) %>% 
  filter(length > 1500)
  

