# REINDEER: efficient indexing of k-mer presence and abundance in sequencing datasets

   * [Motivation](#motivation)
   * [Installation](#installation)
      * [Requirements](#requirements)
      * [Quick start](#quick-start)
   * [Index construction](#index-construction)
      * [Starting with read files (raw or gzipped fasta/fastq)](#starting-with-read-files-raw-or-gzipped-fastafastq)
      * [Starting with De Bruijn graph files (raw or gzipped fasta files)](#starting-with-de-bruijn-graph-files-raw-or-gzipped-fasta-files)
   * [Query fasta files](#query-fasta-files)
   * [Index and query k-mers presence/absence only (not abundances)](#index-and-query-k-mers-presenceabsence-only-not-abundances)
   * [Output](#output)
      * [k-mer abundances](#k-mer-abundances)
      * [k-mer presence/absence](#k-mer-presenceabsence)
   * [Beta options](#beta-options)
      * [query the index on the disk instead of loading the index in-ram](#query-the-index-on-the-disk-instead-of-loading-the-index-in-ram)
      * [log counts/quantized counts](#log-countsquantized-counts)
      * [input paired-end reads (to bcalm)](#input-paired-end-reads-to-bcalm)
   * [Reproduce the manuscript's results](#reproduce-the-manuscripts-results)
   * [Citation](#citation)




# Motivation

REINDEER builds a data-structure that indexes k-mers and their abundances in a collection of datasets (raw RNA-seq or metagenomic reads for instance).
Then, a sequence (FASTA) can be queried for its presence and abundance in each indexed dataset.
While other tools (e.g. SBT, BIGSI) were also designed for large-scale k-mer presence/absence queries, retrieving abundances was so far unsupported (except for single datasets, e.g. using some k-mer counters like KMC, Jellyfish). REINDEER combines fast queries, small index size, and low memory footprint during indexing and queries. We showed it allows to index 2585 RNA-seq datasets (~4 billions k-mers) using less than 60GB of RAM and a final index size lower than 60GB on the disk.
Then, a REINDEER index can either be queried on disk (experimental feature, low RAM usage) or be loaded in RAM for faster queries.

<img src="./Images/reindeer.png" alt="drawing" width="400"/>

Note on presence/absence queries: REINDEER supports this type of queries, although other data structures are more fit for this task. See for instance:

* [HowDeSBT](https://www.ncbi.nlm.nih.gov/pubmed/31504157)
* [Vari-Merge](https://academic.oup.com/bioinformatics/article/35/14/i51/5529124)
* [BIGSI](https://www.ncbi.nlm.nih.gov/pubmed/30718882)
* [BiFrost](https://www.biorxiv.org/content/10.1101/695338v2)

to name a few.

# Installation

## Requirements
* GCC >= 4.8
* CMAKE >  3.10.0

To install, first clone the project:

`git clone --recursive https://github.com/kamimrcht/REINDEER.git`

Then:

`cd REINDEER`

`sh install.sh`

or

`make`

Test can be run:

`sh test.sh`

### Compilation tips

If REINDEER gives a `Error: no such instruction` during compilation, try replacing `-march=native -mtune=native` by `-msse4` in the file `makefile`. If this did not work, pleaes file an issue.

## Quick start
Have a look at the file of file format in `test/fof_unitigs.txt`.
Then build the index:

`./Reindeer --index -f test/fof_unitigs.txt -o quick_out`

and query:

`./Reindeer --query -q test/query_test.fa -l quick_out -o quick_query`

Results should be in `quick_query/query_results/out_query_Reindeer0.out`.

Help:

`./Reindeer --help`

# Index construction

## Starting with read files (raw or gzipped fasta/fastq)

Make sure you have installed Bcalm by doing `sh install`.
Let's assume you work with two files, `reads_1.fastq` and `reads_2.fastq`.
The first thing needed to is to create a file of file (fof) that record the path to the reads.
An example file can be found here: `test/fof.txt`
Then this file of file is provided for the index construction:

`./Reindeer --index -f test/fof.txt --bcalm`

By default, the output will be written in `output_reindeer`.
If you want to change the output, you can use the `-o` option:

`./Reindeer --index -f test/fof.txt --bcalm -o output_dir`

## Starting with De Bruijn graph files (raw or gzipped fasta files)

Let's assume you have already computed a De Bruijn graph for each dataset.
You can provide a file of file of each unitig file (fasta) instead of the read files.
This allows to skip the first Bcalm step.

`./Reindeer --index -f test/fof_unitigs.txt`


# Query fasta files

Simply provide the fasta query file (**single line**) to Reindeer using `-q`, along with the directory of index files that were generated during index construction (`output_reindeer` by default, can changed with `-o`), using `-l':

`./Reindeer --query -q test/query_test.fa -l output_reindeer`


# Index and query k-mers presence/absence only (not abundances)

By default, REINDEER records k-mers abundances in each input dataset.
In order to have k-mer presence/absence instead of abundance per indexed dataset, use `--nocount` option.

`./Reindeer --index -f test/fof_unitigs.txt --nocount`
`./Reindeer --query -l output_reindeer -q test/query_test.fa --nocount`



# Output

By default, query outputs are written in `output_reindeer/query_results/`

## k-mer abundances

The output of REINDEER looks like:

```    
>SRR10092187.6 HISEQ:815:HK2NNBCXY:1:1101:1224:2136 length=51 66 2:4 
```

Counts are reported in a particular format. In this example, `66` means that the queried k-mers correspond to a region which abundance is 66 in dataset 1. In the third column, we observe that two numbers were given, separated by a `:` It means that the sequence spans 2 unitigs in the graph of dataset 2, so REINDEER reported the counts in these 2 unitigs.

## k-mer presence/absence

When launched with `--nocount`, an example line from Reindeer's query is:

```    
>SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51 dataset1:100% dataset3:71%
```

A query sequence appears in the output if more than `t`% k-mers were found an indexed dataset.
Headers from the initial files are used to refer to query sequences (for instance `>SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51`).

Then, the indexed datasets it appears in are written, separated by spaces. The numbers (i.e., `1` and `3` in `dataset1:100% dataset3:71%`) refer to the rank (starting at 1) of the datasets in the inital fof file.

In this example, it means than `t`% k-mers from sequence `>SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51` were in the first and third sample written in the fof file.

Finally, we can see that 100% of the k-mers from the queried sequence appear in dataset `1`, and only 71% appear in dataset `3`.



# Beta options

## query the index on the disk instead of loading the index in-ram

First build the index:

`./Reindeer --index --disk-query -f fof_unitigs.txt`

Then query:

`./Reindeer --query --disk-query -l output_reindeer -q query.fa`


## log counts/quantized counts

`./Reindeer --index --log-count -f fof_unitigs.txt`

`./Reindeer --index --quantization -f test/fof_unitigs.txt`


## input paired-end reads (to bcalm)


`./Reindeer --index --paired-end --bcalm  -f fof.txt`

# Reproduce the manuscript's results

We provide a page with scripts to reproduce the results we show in our manuscript, link [here](https://github.com/kamimrcht/REINDEER/tree/master/reproduce_manuscript_results).


# Citation

Access to the preprint: [REINDEER: efficient indexing of k-mer presence and abundance
in sequencing datasets](https://www.biorxiv.org/content/10.1101/2020.03.29.014159v1.full.pdf)

Citation:
```Bibtex
@article{marchet2020reindeer,
  title={REINDEER: efficient indexing of k-mer presence and abundance in sequencing datasets},
  author={Marchet, Camille and Iqbal, Zamin and Gautheret, Daniel and Salson, Mika{\"e}l and Chikhi, Rayan},
  booktitle={ISMB},
  year={2020},
}
```
