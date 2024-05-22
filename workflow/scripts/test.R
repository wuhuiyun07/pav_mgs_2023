library(dplyr)
library(tidyverse)


getwd()

sample<- read_tsv(snakemake@input[["vs2_file"]]) %>% 
        mutate(sampleID = "{sample}") %>% 
        select(sampleID, everything()) %>%
        write_csv(snakemake@output[[1]])

# vs2 <- read_tsv(snakemake@input[["vs2_file"]], col_names = TRUE, show_col_types = FALSE)

# # test <- read_tsv("results/vs2/{sample}.vs2.final-viral-score.tsv")
#         # %>% mutate(sampleID = "{sample}")%>%
#         #   select(sampleID, everything())
#           %>% write_csv("results/test/{sample}_test.csv")
