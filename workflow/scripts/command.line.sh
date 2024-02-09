rename 's/ /_/g' *
# Replace all spaces with underscores (_) in file names in the current directory
rename 's/_L001_R2_001/.R2/g' *

# To replace "old_string" with "new_string" in filenames
find . -type f -name '_L001_R2_001' -exec sh -c 'mv "$1" "${1/_L001_R2_001/.R2}"' _ {} \;

sed -i 's/_L001_R2_001/.R2/g' 


directory="/results/trimmed"

# Loop through the files in the directory
for file in *_L001_R2_001*; do
    mv "$file" "$(echo "$file" | sed 's/_L001_R2_001/.R2/')"
done

# Replace '.txt' extension with '.md' for all files in the current directory
rename 's/\_L001_R1_001.fastq.gz$/R2.fastq.gz/' *_L001_R1_001.fastq.g

# export manual
# in singularity shell
diamond help > diamond_manual.txt
