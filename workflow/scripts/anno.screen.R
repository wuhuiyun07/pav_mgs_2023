library(dplyr)
library(tidyverse)

sample<- read_tsv(snakemake@input[[2]])%>%
            write_tsv(snakemake@output[[1]])



# vs2 <- read_tsv("results/vs2/Wu_23_4_S14/final-viral-score.tsv", col_names = TRUE)
# checkV <- read_tsv("results/checkV-megahit/Wu_23_4_S14/quality_summary.tsv", col_names = TRUE)
# diamond <- read_tsv("results/diamond_tmp/Wu_23_4_S14.tsv", skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)

vs2 <- read_tsv(snakemake@input[["vs2"]], col_names = TRUE)
checkV <- read_tsv(snakemake@input[["checkv"]], col_names = TRUE)
diamond <- read_tsv(snakemake@input[["diamond"]], skip = 3, col_names = FALSE, col_types = cols(), show_col_types = FALSE)


# Rename column names for diamond data
colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "evalue", "bitscore", "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle")


head(vs2)
head(checkV)
head(diamond)

# Filter vs2 data
vs2_screened <- vs2 %>%
  filter(max_score > 0.5) %>%
  filter(length > 2000) %>%
  separate(seqname, into = c("contig_id", "gene"), sep = "\\|\\|")

head(vs2_screened)

# Filter checkV data
checkV_screened <- checkV %>%
  filter(checkv_quality %in% c("Low-quality", "Medium-quality", "High-quality")) %>%
  filter(contig_length > 3000)%>%
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
  select(contig_id, staxids, pident, gene.y, contig_length, bitscore, max_score, max_score_group, viral_genes,checkv_quality, completeness,miuvig_quality,skingdoms,sphylums,sscinames) %>%
  arrange(checkv_quality, miuvig_quality,desc(completeness),max_score_group, desc(skingdoms), desc(bitscore),desc(contig_length)) 
head(screen_result)

# write_csv(screen_result, "results/screen/Wu_23_4_S14.csv")
write_csv(screen_result, snakemake@output[["highscore"]])

### add sample ID to the screen_result###
# Specify the full file path
# file_path <- "results/screen/Wu_26_4_S24/high.score.csv"
file_path <- snakamake@output[["highscore"]] 

# Load the CSV file
data <- read.csv(file_path)

# Extract the folder name
folder_name <- basename(dirname(file_path))

# Add a new column 'SampleID' as the first column
data <- cbind(SampleID = folder_name, data)

# Save the updated data to a new CSV file
write.csv(data, file = snakemake@output[["highscorewID"]], row.names = FALSE)


## add sampleID in highscore file
## add tzxonomizr result in highscore file
## add filter type, collection location in highscore file

# Load necessary libraries
library(dplyr)

# Set the directory path containing your CSV files
folder_path <- "results/screen/Wu_29_2_S22/"

# Get the folder name
folder_name <- basename(folder_path)

# Read your CSV file
csv_file <- "high.score.csv" # Your CSV file name
data <- read.csv(file.path(folder_path, csv_file))

# Add the folder name as a new column
data <- data %>%
  mutate(sampleID = folder_name)

# Move the FolderName column to the first position
data <- data %>%
  select(sampleID, everything())

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
