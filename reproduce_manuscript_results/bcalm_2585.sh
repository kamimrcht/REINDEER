### 1/ download each dataset, and launch bcalm with appropriate filters using file lists from HowDe
mkdir -p bcalm_results
mkdir -p tmp
> sample_all.in
#get fastq.gz and launch bcalm on each file
while read -r filename threshold; do
    prefetch "$filename"
    parallel-fastq-dump --sra-id "$filename" --threads 20 --split-spot --gzip --tmpdir tmp
    loc=$(realpath $filename)
    echo "$loc.fastq.gz" > $filename.in
    bcalm -in $filename.in -kmer-size 21 -abundance-min "$threshold" -out-dir bcalm2 -out-tmp bcalm2 -nb-cores 20
    mv $filename.unitigs.fa bcalm_results
    ### 1/ download each dataset, and launch bcalm with appropriate filters using file lists from HowDe
rm -f ../tmp/sra/*.sra

mkdir -p bcalm_results
mkdir -p tmp
> sample_all.in
#get fastq.gz and launch bcalm on each file
while read -r filename threshold; do
    #echo $filename
    prefetch "$filename"
    parallel-fastq-dump --sra-id "$filename" --threads 20 --split-spot --gzip --tmpdir tmp
    loc=$(realpath $filename)
    #loc=$(realpath $filename)
    echo "$loc.fastq.gz" > $filename.in
    bcalm -in $filename.in -kmer-size 21 -abundance-min "$threshold" -out-dir bcalm2 -out-tmp bcalm2 -nb-cores 20
    mv $filename.unitigs.fa bcalm_results
    ls $PWD/bcalm_results/$filename.unitigs.fa >> fof_reindeer.lst
done < file_list



### cleaning
rm -f *.glue*
rm -f bcalm2/*

