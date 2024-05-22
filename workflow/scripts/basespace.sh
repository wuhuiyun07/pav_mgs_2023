#download basespace cli
#https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview
~/bs download project -i 414267855 -o /work/wuhuiyun/pavB2/rawdata --extension=fastq.gz
find /work/wuhuiyun/pavB2/rawdata -type f -name "*.fastq.gz" -exec mv {} /work/wuhuiyun/pavB2/fastq \;
