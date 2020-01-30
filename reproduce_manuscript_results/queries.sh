# Refseq transcripts from UCSC's hg38
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/refMrna.fa.gz
gunzip refMrna.fa.gz
# an example of sequence batch
seqkit sample -n 10000 refMrna.fa > sample_1000_transcripts.fa
# query with Reindeer
./Reindeer --query --count -k 21 -q sample_1000_transcripts.fa -l out_reindeer_2585_counts #out_reindeer_2585_counts is the directory where Reindeer's index was written
