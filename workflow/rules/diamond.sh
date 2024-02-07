
module load anaconda3/2020.07
conda update diamond
Source activate diamond-env #created a virtual environment to run diamond


ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip #Download NCBI taxonomic names zip file, you have to unzip this for the viral database to work
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz #download NCBI accession #s with taxonomic IDs from NCBI

##make viral database to run in Diamond with taxonomic names from NCBI, This took 4888s using 40 threads on HPC to create

diamond makedb --in /lustre/project/taw/kvigil/Reference/refseq_viral_proteins/viral.1.protein.faa.gz -d viralprotein --taxonmap /lustre/project/taw/kvigil/Reference/prot.accession2taxid.FULL.gz --taxonnames  /lustre/project/taw/kvigil/Reference/names.dmp --taxonnodes
/lustre/project/taw/kvigil/Reference/nodes.dmp

## execute diamond with viral database that has all the taxonomic names downloaded
diamond blastx -d /lustre/project/taw/kvigil/Reference/viralprotein.dmnd -q consensus.fasta -o ONR030223oysterspooledblastxdenovo.tsv --ultra-sensitive -f 6 qseqid sseqid pident length mismatch evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle


