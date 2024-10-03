#download basespace cli
#https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
~/bs download project -i 414267855 -o /work/wuhuiyun/pavB2/rawdata --extension=fastq.gz
find /work/wuhuiyun/pavB2/rawdata -type f -name "*.fastq.gz" -exec mv {} /work/wuhuiyun/pavB2/fastq \;

~/bs download project -n 06212024_RNASeq_Dr.Wu_25x -o /work/wuhuiyun/pavB2/rawdata --extension=fastq.gz
bs download project -n 06212024_RNASeq_Dr.Wu_25x -o . --extension=fastq.gz
bs download project -n 06272024_RNASeq_Dr.Wu_26x_2024-06-26T18_32_11_1c0bccf -o . --extension=fastq.gz
# The project ID is often included in the URL of the project details page. 
# It typically follows the format basespace.illumina.com/projects/{project_id}.
bs list projects # find project id from command
# pavB4
bs download project -i 423009914 -o . --extension=fastq.gz
# pavB3
bs download project -i 422795373 -o . --extension=fastq.gz
find . -name "*.fastq.gz" -exec mv {} destination_folder/ \; # moving files resursivly from subdirectories