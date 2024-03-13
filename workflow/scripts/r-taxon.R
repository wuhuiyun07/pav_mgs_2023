

######R code#####
# packages available in r-tidy conda environment
library(taxonomizr)
# library(tidyverse)
# library(tibble)

# diamond<-read_tsv(snakemake@input[["diamond_file"]])
# diamond<-read.tsv("results/diamond/24_5_S20.tsv.gz",skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)
# colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")

# taxID<-diamond_output$staxids
#change to string format 

# awk -F'\t' 'NR==5' ONR10623barcode03.tsv #select fifth row
# awk -F'\t' '{print $9}' ONR10623barcode03.tsv >> taxaID.csv #select ninth column

# awk -F'\t' '{print $8}' results/diamond/24_5_S20.tsv >> results/diamond/24_5_S20.taxaID.csv #select eighth column for staxids

taxaID <- read.table("results/diamond/24_4_S19.tsv", sep = "\t", header = FALSE)[,8]

taxa<-getTaxonomy(taxaID,'../r-taxon/accessionTaxa.sql')
colnames(taxa)[1] <- "staxids"
# print(taxa)

write.csv(taxa, snakemake@output[["taxon_file"]])
# taxonomy  <- write.csv(taxa, "results/diamond/24_4_S19.taxa.csv")


