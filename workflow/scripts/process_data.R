
# Define function to process each input file
process_file <- function(input_file, output_file) {
    # Read input file
   data <- read.csv(input_file)

    # Extract the base filename without the directory path and extension
    base_filename <- basename(input_file) 

    # Extract the prefix from the base filename
    prefix <- sub("_screened.csv", "", base_filename)

    # Add prefix as the last column with the column name "sampleID"
    data$sampleID <- prefix
    data <- data[, c("sampleID", setdiff(names(data), "sampleID"))]

    # Write the modified data to the output file
    write.csv(data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}


input_file <- snakemake@input[["file"]]
output_file <- snakemake@output[["individual"]]

# Process the input file
process_file(input_file, output_file)

# data<-read.csv("results/annotation/16_2_S2_screened.csv")
# dim(data)
# colnames(data)
