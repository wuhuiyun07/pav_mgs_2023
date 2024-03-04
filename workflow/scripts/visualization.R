library(tidyverse)
library(dplyr)
library(tibble)
library(ggplot2)

setwd("/Users/huiyunwu/Desktop/Virus_particle/pav_mgs_2023")
getwd()

# file_list <- list.files(pattern = "results/vs2/test/*.tsv")
# vs2_combined <- map_dfr(file_list, ~read_tsv(.x))
# vs2 <-
#   list.files( pattern = "results/vs2/test/*.tsv") %>%
#   map_dfr(~read_tsv(.))

# checkV

# diamond <- read_csv(file="results/diamond/16_5_S5.tsv", col_names = FALSE)
# write_tsv(diamond, file= "results/diamond/16_5_S5_v2.tsv")
# diamond <- read_csv(file="results/diamond/16_5_S5_v2.tsv",col_names = FALSE, separate("\t"))



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
  

