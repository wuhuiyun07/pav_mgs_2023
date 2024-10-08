#!/usr/bin/env Rscript

# conda env export > workflow/rules/taxonomizr.yml 

library(taxonomizr)

# diamond<-read_tsv(snakemake@input[["diamond_file"]])
# diamond<-read.tsv("results/diamond/24_5_S20.tsv.gz",skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)
# colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")

# taxID<-diamond_output$staxids
#change to string format 

# awk -F'\t' 'NR==5' ONR10623barcode03.tsv #select fifth row
# awk -F'\t' '{print $9}' ONR10623barcode03.tsv >> taxaID.csv #select ninth column

# awk -F'\t' '{print $8}' results/diamond/24_5_S20.tsv >> results/diamond/24_5_S20.taxaID.csv #select eighth column for staxids
getwd()
# taxaID <- read.table("results/diamond/24_4_S19.tsv", sep = "\t", header = FALSE)[,8]
# taxaID <- read.table(snakemake@input[[diamond_file]], sep = "\t", header = FALSE)[,8]
# taxaID <- read.csv("results/combined.csv", header = TRUE)[,2]
#results <- read.table("diamond.Wu_24_5_S20.csv",sep = "\t")
results <- read.table("results/diamond_tmp/Wu_23_5_S15.tsv", sep = "\t")
dim(results)
taxaid<-as.character(results[,8])
taxa <-getTaxonomy(taxaid, '/data/lab/hwu/resources/accessionTaxa.sql')
write.csv(taxa, 'results/taxa/Wu_23_5_S15.taxa.csv')


# colnames(taxaID)
# taxa<-getTaxonomy(taxaID,'../r-taxon/accessionTaxa.sql')
# colnames(taxa)[1] <- "staxids"
# print(taxa)

# write.csv(taxa, snakemake@output[["taxon_file"]])
# taxonomy  <- write.csv(taxa, "results/diamond/24_4_S19.taxa.csv")
write.csv(taxa,"results/bigtable.csv")





1   # scripts/my_analysis.R 
  1 #!/usr/bin/env Rscript 
  2  
  3 library(taxonomizr) 
  4  
  5 getwd() 
  6 results <- read.table("./results/diamond_tmp/{sample}.tsv", sep = "\t") 
  7 dim(results) 
  8 taxaid<-as.character(results[,8]) 
  9 taxa <-getTaxonomy(taxaid, '/data/lab/hwu/resources/accessionTaxa.sql') 
 10 write.csv(taxa, 'results/taxa/{sample}/taxa-blast-tmp.csv') 
 11  
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
~                                                                                                                                    
-- VISUAL --                                                                                             12        1,6           All