
# fastqc
fastqc -t 12 -o out_path sample1_1.fq sample1_2.fq

#trimmomatic
trimmomatic PE ../raw_data/sludge/AS-Cu-1_1.fastq ../raw_data/sludge/AS-Cu-1_2.fastq -baseout ./sludge/AS-Cu-1.fastq ILLUMINACLIP:/home/xiarong/miniconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50 -phred33


#fastuniq
fastuniq -i AS-Cu-1_input.fastuniq -o AS-Cu-1_uniq.R1.fq -p AS-Cu-1_uniq.R2.fq

#megahit
megahit -t 40 -1 ../1.3_fastuniq_result/sludge/AS-Cu-1_uniq.R1.fq -2 ../1.3_fastuniq_result/sludge/AS-Cu-1_uniq.R2.fq -o ./AS-Cu-1_assembly --k-min 35 --k-max 95 --k-step 20 --min-contig-len 500 -m 0.1

#Virsorter2
virsorter run -w AS-Cu-1_vs2 -i /96t/xiarong/MP/2.1_megahit_assembly_result/sludge/final.contig/AS-Cu-1_final.contigs.fa --min-length 5000 -j 4 all

#Deepvirfinder
python dvf.py -i /96t/xiarong/MP/2.1_megahit_assembly_result/sludge/final.contig/AS-Cu-1_final.contigs.fa -o .  -l 5000 -c 40

#checkv
checkv end_to_end input_file.fna output_fiel -t 40

#mmseqs
mmseqs easy-cluster all_virus_contig.fa virus_contig_cluster --min-seq-id 0.95 --cov-mode 0 -c 0.85 tmp --threads 15

#coverm
coverm contig --coupled ./sludge/AS-Cu-1_uniq.R1.fq ./sludge/AS-Cu-1_uniq.R2.fq --reference virus_uniq_10kb.fasta -m rpkm -o coverm_result1 -t 20

#genomad
genomad annotate virus_uniq_10kb.fasta genomad_output  genomad_db

#vibrant
VIBRANT_run.py -i ../4_virus_abundance/virus_uniq_10kb.fasta  -t 15 -virome

#Dramv
virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i virus_uniq_10kb.fasta -w ./ --min-length 1500 --min-score 0.5 -j 50 all 
DRAM-v.py annotate -i final-viral-combined-for-dramv.fa -v viral-affi-contigs-for-dramv.tab -o ./dramv_out --skip_trnascan --threads 20 --min_contig_size 1000

#metaWRAP
metaWRAP binning --metabat2 --maxbin2 --concoct -a /cardboard-1_novirus_contig.fa -o output_dir cardboard-1_clean_1.fastq cardboard-1_clean_2.fastq -t 40

#dRep
dRep dereplicate drep_result -g all_bin/*.fa -p 20 -comp 70 -con 5 -sa 0.95

#gtdbtk
gtdbtk classify_wf --genome_dir /96t/xiarong/MP/8_bining/drep_bin/drep_result/dereplicated_genomes/ --out_dir ./ --extension fa --prefix bin --cpus 8

#eggnog-mapper
python emapper.py --data_dir /home/public/miniconda/eggnog/ -i bin_gene_cluster_result_rep_seq.fasta --output eggnog_out -m diamond --cpu 12

###host-prediction
#CRISPR
crt-mod -i /96t/xiarong/MP/8_bining/bin_ARG/bin_all_contig.fa fasta -o /96t/xiarong/MP/10_host_prediction/CRISPR/bin_spacer.fasta -f --ignore-empty --min-nr=1 --min-rl=20 --max-rl=50 --min-sl=20 --max-sl=60 --search-wl=7 -a 40

#tRNA
tRNAscan-SE -B -o tRNA.out -f rRNA.ss -m tRNA.stats -a tRNA.fa --detail --thread 40 virus_uniq_10kb.fasta

#homo
blastn -query /96t/xiarong/MP/8_bining/bin_ARG/bin_all_contig.fa -db ../CRISPR/phage_database/phage -out homo.xls -evalue 1e-5 -outfmt 6


#MetaCHIP
cut -f1,2 /96t/xiarong/MP/8_bining/bin_classfication_1208/classify/bin.bac120.summary.tsv > GTDB_classification_1208.tsv
MetaCHIP PI -p MP_1208 -r pcofg -t 40 -o ./MetaCHIP_PI_22 -i /96t/xiarong/MP/8_bining/drep_bin/all_bin_RN_22 -x fa -taxon GTDB_classification_1208.tsv
MetaCHIP BP -p MP_1208 -r pcofg -t 40 -o ./MetaCHIP_PI_22

#sortmerna
sortmerna --ref rRNAdatabase/smr_v4.3_default_db.fasta --reads tM-pbat-Bu3_uniq.R1.fq --reads tM-pbat-Bu3_uniq.R2.fq -workdir /96t/xiarong/MP/17_metatranscriptome/sortmerna_result/run2 --out2 --aligned /96t/xiarong/MP/17_metatranscriptome/sortmerna_result/tM-pbat-Bu3/aligned.fasta --sam --num_alignments 1 --fastx --other /96t/xiarong/MP/17_metatranscriptome/sortmerna_result/tM-pbat-Bu3/other.fasta --threads 15






