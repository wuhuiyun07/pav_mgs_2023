# fetch testing data
wget -O test.fa https://raw.githubusercontent.com/jiarong/VirSorter2/master/test/8seq.fa
# run classification with 4 threads (-j) and test-out as output diretory (-w)
virsorter run -w test.out -i test.fa --min-length 1500 -j 4 all
ls test.out