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
library(ggplot2)
# library(glue)
# library(lubridate)

# Set working directory if needed
# setwd("/Users/huiyunwu/Desktop/Virus_particle/pav_mgs_2023")

# Print working directory
getwd()

vs2_file <- list.files("results/vs2", full.name = TRUE)
checkV_file <- list.files("results/checkV", full.name = TRUE)
diamond_file <- list.files("results/diamond_vs2", full.name = TRUE)

# Input and output files
# vs2_file <- snakemake@input[["vs2_file"]]
# checkV_file <- snakemake@input[["checkV_file"]]
# diamond_file <- snakemake@input[["diamond_file"]]
# output_file <- snakemake@output[["annotation"]]


# Read input files
vs2 <- read_tsv(vs2_file, col_names = TRUE)
checkV <- read_tsv(checkV_file, col_names = TRUE)
diamond <- read_tsv(diamond_file, skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)

# Rename column names for diamond data
colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")


# Filter vs2 data
vs2_screened <- vs2 %>%
  filter(max_score > 0.5) %>%
  filter(length > 1000) %>%
  separate(seqname, into = c("contig_id", "gene"), sep = "\\|\\|")


# Filter checkV data
checkV_screened <- checkV %>%
  filter(checkv_quality %in% c("Low-quality", "Medium-quality", "High-quality")) %>%
  filter(contig_length > 1000)%>%
  arrange(checkv_quality,miuvig_quality, desc(completeness))


# Join vs2 and checkV data
contigs_for_diamond <- full_join(vs2_screened, checkV_screened, by = "contig_id")

# Filter and sort diamond data
diamond_combined <- diamond %>%
  filter(pident > 30 & bitscore > 50 & length > 30) %>%
  arrange(desc(bitscore), evalue, desc(length), desc(pident)) %>%
  separate(qseqid, into = c("contig_id", "gene"), sep = "\\|\\|")

# print(colnames(annotation))
# 1] "contig_id"           "gene.x"              "dsDNAphage"          "ssDNA"               "NCLDV"               "RNA"                 "lavidaviridae"       "max_score"          
#  [9] "max_score_group"     "length.x"            "hallmark"            "viral"               "cellular"            "contig_length"       "provirus"            "proviral_length"    
# [17] "gene_count"          "viral_genes"         "host_genes"          "checkv_quality"      "miuvig_quality"      "completeness"        "completeness_method" "contamination"      
# [25] "kmer_freq"           "warnings"            "gene.y"              "sseqid"              "pident"              "length.y"            "mismatch"            "evalue"             
# [33] "bitscore"            "staxids"             "sscinames"           "sskingdoms"          "skingdoms"           "sphylums"            "stitle"  
# Perform left join and write output
screen_result <- left_join(contigs_for_diamond, diamond_combined, by = "contig_id") %>%
  distinct(contig_id, .keep_all = TRUE) %>% 
  select(contig_id, gene.y, contig_length, bitscore, max_score, max_score_group, viral_genes,checkv_quality, completeness,miuvig_quality,skingdoms,sphylums,sscinames) %>%
  arrange(checkv_quality, miuvig_quality,desc(completeness),max_score_group, desc(skingdoms), desc(bitscore),desc(contig_length))   %>%
  write_tsv("results/screen_result.tsv")

# contig_length
length<- screen_result %>% 
  count(contig_id, contig_length, completeness) %>%
  arrange(desc(contig_length))

ggplot(screen_result, aes(x=contig_length)) + 
  geom_histogram()

ggplot(screen_result, aes(x=completeness)) + 
  geom_histogram()

ggplot(screen_result, aes(x= viral_genes))+
  geom_histogram()


ggplot(screen_result %>% count(max_score_group), aes(x= max_score_group, y= n))+
  geom_bar(stat = "identity")

ggplot(screen_result %>% count(skingdoms), aes(x= skingdoms, y= n))+
  geom_bar(stat = "identity")

ggplot(screen_result %>% count(sphylums), aes(x= sphylums , y= n))+
  geom_bar(stat = "identity")

ggplot(screen_result %>% count(checkv_quality), aes(x= checkv_quality, y= n))+
  geom_bar(stat = "identity")



  count(max_score_group,skingdoms, miuvig_quality, checkv_quality) #%>%


  count(gene.y, miuvig_quality, checkv_quality) #%>%
  count(miuvig_quality) 
  count(checkv_quality)
  # %>% count(checkv_quality)




### FUNCTION
screen_data <- function(vs2_file, checkV_file, diamond_file) {
  library(tidyverse)

  # Read input files
  vs2 <- read_tsv(vs2_file, col_names = TRUE)
  checkV <- read_tsv(checkV_file, col_names = TRUE)
  diamond <- read_tsv(diamond_file, skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)

  # Rename column names for diamond data
  colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")

  # Filter vs2 data
  vs2_screened <- vs2 %>%
    filter(max_score > 0.5) %>%
    filter(length > 1000) %>%
    separate(seqname, into = c("contig_id", "gene"), sep = "\\|\\|")

  # Filter checkV data
  checkV_screened <- checkV %>%
    filter(checkv_quality %in% c("Low-quality", "Medium-quality", "High-quality")) %>%
    filter(contig_length > 1000)

  # Join vs2 and checkV data
  contigs_for_diamond <- full_join(vs2_screened, checkV_screened, by = "contig_id")

  # Filter and sort diamond data
  diamond_combined <- diamond %>%
    filter(pident > 30 & bitscore > 50 & length > 30) %>%
    arrange(desc(bitscore), evalue, desc(length), desc(pident)) %>%
    separate(qseqid, into = c("contig_id", "gene"), sep = "\\|\\|")

  # Perform left join
  screen_result <- left_join(contigs_for_diamond, diamond_combined, by = "contig_id") %>%
    distinct(contig_id, .keep_all = TRUE) %>% 
    arrange(checkv_quality, max_score_group, desc(skingdoms), desc(bitscore)) 

  return(screen_result)
}

# Example usage:
# result <- screen_data("results/vs2_file.tsv", "results/checkV_file.tsv", "results/diamond_file.tsv")

