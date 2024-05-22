library(dplyr)
library(tidyverse)

sample<- read_tsv(snakemake@input[[2]])%>%
            write_tsv(snakemake@output[[1]])


# vs2_file <- list.files("results/vs2", full.name = TRUE)
# checkV_file <- list.files("results/checkV", full.name = TRUE)
# diamond_file<- list.files("results/diamond", full.name = TRUE)


vs2 <- read_tsv(snakemake@input[["vs2_file"]], col_names = TRUE)
checkV <- read_tsv(snakemake@input[["checkV_file"]], col_names = TRUE)
diamond <- read_tsv(snakemake@input[["diamond_file"]], skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)

# diamond <- read_tsv("results/diamond/24_5_S20.tsv", skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)
# Rename column names for diamond data
colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")

# test <- read_tsv("results/vs2/29_5_S25.vs2.final-viral-score.tsv")%>%
#           mutate(sampleID = "29_5_S25")%>%
#           select(sampleID, everything())

head(vs2)
head(checkV)
head(diamond)

# Filter vs2 data
vs2_screened <- vs2 %>%
  filter(max_score > 0.5) %>%
  filter(length > 1000) %>%
  separate(seqname, into = c("contig_id", "gene"), sep = "\\|\\|")

head(vs2_screened)

# Filter checkV data
checkV_screened <- checkV %>%
  filter(checkv_quality %in% c("Low-quality", "Medium-quality", "High-quality")) %>%
  filter(contig_length > 1000)%>%
  arrange(checkv_quality,miuvig_quality, desc(completeness))

# print(checkV_screened)
# Join vs2 and checkV data
contigs_for_diamond <- full_join(vs2_screened, checkV_screened, by = "contig_id")
# print(contigs_for_diamond)

# Filter and sort diamond data
diamond_combined <- diamond %>%
  filter(pident > 30 & bitscore > 50 & length > 30) %>%
  arrange(desc(bitscore), evalue, desc(length), desc(pident)) %>%
  separate(qseqid, into = c("contig_id", "gene"), sep = "\\|\\|")
head(diamond_combined)

# print(colnames(annotation))
# 1] "contig_id"           "gene.x"              "dsDNAphage"          "ssDNA"               "NCLDV"               "RNA"                 "lavidaviridae"       "max_score"          
#  [9] "max_score_group"     "length.x"            "hallmark"            "viral"               "cellular"            "contig_length"       "provirus"            "proviral_length"    
# [17] "gene_count"          "viral_genes"         "host_genes"          "checkv_quality"      "miuvig_quality"      "completeness"        "completeness_method" "contamination"      
# [25] "kmer_freq"           "warnings"            "gene.y"              "sseqid"              "pident"              "length.y"            "mismatch"            "evalue"             
# [33] "bitscore"            "staxids"             "sscinames"           "sskingdoms"          "skingdoms"           "sphylums"            "stitle"  
# Perform left join and write output
screen_result <- left_join(contigs_for_diamond, diamond_combined, by = "contig_id") %>%
  distinct(contig_id, .keep_all = TRUE) %>% 
  # select(contig_id, staxids, pident, gene.y, contig_length, bitscore, max_score, max_score_group, viral_genes,checkv_quality, completeness,miuvig_quality,skingdoms,sphylums,sscinames) %>%
  arrange(checkv_quality, miuvig_quality,desc(completeness),max_score_group, desc(skingdoms), desc(bitscore),desc(contig_length))  %>%
  write_csv(snakemake@output[["individual"]])
head(screen_result)

getwd()

# combined_file <- list.files("results/annotation", full.name = TRUE)

# combined<-read_csv(combined_file)%>%
#           write_csv("results/combined2.csv")
# colnames(combined)
#  [1] "contig_id"           "gene.x"              "dsDNAphage"          "ssDNA"               "NCLDV"              
#  [6] "RNA"                 "lavidaviridae"       "max_score"           "max_score_group"     "length.x"           
# [11] "hallmark"            "viral"               "cellular"            "contig_length"       "provirus"           
# [16] "proviral_length"     "gene_count"          "viral_genes"         "host_genes"          "checkv_quality"     
# [21] "miuvig_quality"      "completeness"        "completeness_method" "contamination"       "kmer_freq"          
# [26] "warnings"            "gene.y"              "sseqid"              "pident"              "length.y"           
# [31] "mismatch"            "evalue"              "bitscore"            "staxids"             "sscinames"          
# [36] "sskingdoms"          "skingdoms"           "sphylums"            "stitle"             
