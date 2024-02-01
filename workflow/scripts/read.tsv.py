# python
import pandas as pd
# Specify the file path to your TSV file
file_path = "path_to_your_file/file_name.tsv"

# Read the TSV file into a pandas DataFrame
data = pd.read_csv(file_path, sep='\t')

# Display the DataFrame (optional)
print(data)

# bash: Using cat and awk:
cat your_file.tsv | awk '{print $1, $2, $3}' 
cat results/16_4_S4.tsv | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}'
# This command will display the content of the TSV file, and awk will print the columns you specify.

# Using cut:
cut -f 1,2,3 your_file.tsv
# The cut command allows you to extract specific columns from a file based on a delimiter, which is tab by default.

# using awk with printf for better formatting:
awk -F'\t' '{printf "Column 1: %s, Column 2: %s, Column 3: %s\n", $1, $2, $3}' your_file.tsv

# using column command:
column -t -s $'\t' your_file.tsv
