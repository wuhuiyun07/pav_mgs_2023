directory="/results/trimmed"

for file in "$directory"/*; do
    if [ -f "$file" ]; then
        new_file=$(echo "$file" | sed "s/_L001_R2_001/.R2/")
        mv "$file" "$new_file"
    fi
done
