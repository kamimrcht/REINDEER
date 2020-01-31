Reproduce the results presented in REINDEER's manuscript
========================================================

# Pre-requisites
* [BCALM2](https://github.com/GATB/bcalm)
* [parallel-fastq-dump](https://github.com/rvalieris/parallel-fastq-dump)
* [seqkit](https://bioinf.shenwei.me/seqkit/download/)

# Install REINDEER
`git clone --recursive https://github.com/kamimrcht/REINDEER.git`

`cd REINDEER`

`sh install.sh`


# Download datasets
We dowloaded datasets from SRA and computed a de Bruijn graph per dataset with BCALM2 using this script:

` sh bcalm_2585.sh `

After its execution, we have 2585 de Bruijn graphs corresponding to each dataset and a file of file (`fof_reindeer.lst`) recording paths to each of those de Bruijn graphs.


# Index construction
REINDEER is launched on all those graphs. It builds an index that will record abundances in each dataset :

`sh index_construction.sh`

Index construction step's output file is `out_reindeer_2585_counts`. 

Results shown in Table 3 of the manuscript were generated using this script. Each line enables to construct the index using different options (quantized counts, logs, presence/absence only).

In order to reproduce similar results to Table 4, we selected samples of 10, 100, 500 and 1000 datasets from `fof_reindeer.lst` and launched Reindeer (as in line 2 of `index_construction.sh`, changing `-f [new_fof]` parameter, with `[new_fof]` being the sample file of file.


# Query
We used the previously built index to perform queries of transcripts.
We got transcripts sequences by downloading Refseq transcripts from UCSC's hg38, and created batches with seqkit, then queried those batches files with REINDEER using this script (an example for a file of 1000 transcripts is given):

`sh queries.sh`

Results will be available in `results_query_1000_on_2585/query_results/out_query_Reindeer0.out`

# Query of oncogenes

We used the metadata associated with each sample as a proxy to determine if the sample is a cancer sample or not. Unfortunately as those metadata are not defined in a normalised way, we can only perform rudimentary full-text search. Such a search will also retrieve control samples that belong to a cancer experiment. Also, some experiments are small RNA-seq sequencing runs, in which our selected oncogenes won't be found.
In spite of those limitations our expectations are that the ''cancer-related'' group of datasets will contain a higher proportion of true RNA-seq cancer datasets than the ''non cancer-related'' group.

We systematically searched the term ''cancer'' or ''carcinogen'' in the metadata of the 2,585 RNA-seq and found it in 1,185 samples. Most of the time (1,112 samples) those samples' metadata also contained the term ''breast''. Thus 43% of the samples are cancer-related. We expect to find an over expression of oncogenes in those 1,185 samples compared to the ''`non cancer-related'' ones.
The list of the corresponding samples is accesible here: [cancer-related samples](data/cancer_dataset), [non cancer-related samples](data/no_cancer_dataset).

We selected four oncogenes well known for their role in breast cancers (ERBB2, FOXM1, MYC, and PIK3CA) as well as three tumor suppressor genes (BRCA1, PTEN, TP53) [Perera et l. 2012, Song et al. 2017].
We launched REINDEER on the longest transcript of those genes. The transcript was split in 100bp long sequences, to be able to take into account matches on a few exons. The sequences can be found here: [ERBB2](data/erbb2-split.fa), [FOXM1](data/foxm1.fa), [MYC](data/myc-split.fa), [PIK3CA](data/pik3ca-split.fa), [BRCA1](data/brca1-split.fa), [PTEN](data/pten-split.fa), [TP53](data/tp53-split.fa).
We required at least 78% of the k-mers to be found in the split sequences. Among the sequences we kept the maximal count returned by REINDEER.

