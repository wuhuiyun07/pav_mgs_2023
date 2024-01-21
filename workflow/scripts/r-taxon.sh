
R 
######R code#####
# packages available in r-tidy conda environment
library(taxonomizr)
library(tidyverse)
library(tibble)
library(dplyr)
library(ggplot2)
library(dplyr)

diamond_output<-read_tsv("./t-taxon/dimo-tsv/ONR10623barcode03.tsv")
taxID<-diamond_output[,9]
#change to string format 

awk -F'\t' 'NR==5' ONR10623barcode03.tsv #select fifth row
awk -F'\t' '{print $9}' ONR10623barcode03.tsv >> taxaID.csv #select ninth column

echo "$(cat taxaID.csv)"
#sed add double quotes to each line
sed 's/.*/"&"/' taxaID.csv > taxaID_string.txt




#note this will require a lot of hard drive space, bandwidth and time to process all the data from NCBI
# prepareDatabase('accessionTaxa.sql') #done
# use Accession number
# blastAccessions<-c("Z17430.1","Z17429.1","X62402.1") 
# ids<-accessionToTaxa(blastAccessions,'accessionTaxa.sql')
# getTaxonomy(ids,'accessionTaxa.sql')

# light weight database
# prepareDatabase(getAccessions=FALSE)
# prepareDatabase(types='prot')

# use taxaID
taxaID<-read.csv("taxaID_string.csv", quote = "\"")
paste("\"", taxaID, "\"", sep = "")

data <- read.table("taxaID_string.txt", quote = "\"", sep = "\t", header = TRUE)

taxaID<-read.csv("taxaID_string.csv", quote = T, row.names = F)


test <- read.csv("taxaID.csv", quote = '"')
test <- read.csv("taxaID.csv", as.is = T)
test <- read.csv("taxaID.csv", as.is = T, quote = '"')

quoted_taxaID <- paste('"', taxaID, '"', sep = '')
quoted_vector <- paste('"', numeric_vector, '"', sep = '')

taxaID<-read.csv("taxaID_string.csv")
taxaID<-c("1283077", "2163636","2785027")
taxa<-getTaxonomy(taxaID,'accessionTaxa.sql')
print(taxa)

taxonomy  <- save.csv(taxa, results/anonotation/taxonomy.csv)
