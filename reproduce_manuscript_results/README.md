Reproduce the results presented in REINDEER's manuscript
========================================================

# Pre-requisites
* BCALM2
* parallel-fastq-dump
* seqkit

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

