# REINDEER
REINDEER  REad Index for abuNDancE quERy



# Installation

## Requirements
* GCC >= 4.8
* CMAKE >  3.10.0

To install, first clone the project:

`git clone --recursive https://github.com/kamimrcht/REINDEER.git`

Then:

`cd REINDEER`

`sh install.sh`

## Quick start
Have a look at the file of file format in `test/fof.txt`.
Then build the index:

`./Reindeer --index -f test/fof.txt --bcalm`

and query:

`./Reindeer --query -q test/query_test.fa -l output_reindeer`

Results should be in `output_reindeer/query_results`.

Help:

`./Reindeer --help`

# Index construction

## Starting with read files (fasta/fastq)

Let's assume you work with two files, `reads_1.fastq` and `reads_2.fastq`.
The first thing needed to is to create a file of file (fof) that record the path to the reads.
An example file can be found here: `test/fof.txt`
Then this file of file is provided for the index construction:

`./Reindeer --index -f test/fof.txt --bcalm`

By default, the output will be written in `output_reindeer`.
If you want to change the output, you can use the `-o` option:

`./Reindeer --index -f test/fof.txt --bcalm -o output_dir`

## Starting with De Bruijn graph files

Let's assume you have already computed a De Bruijn graph for each dataset.
You can provide a file of file of each unitig file (fasta) instead of the read files.
This allows to skip the first Bcalm step.

`./Reindeer --index -f test/fof_unitigs.txt`

Now if you have also previously built the main graph from all unitigs' k-mers, you can provide it as well. This allows to skip all Bcalm steps.

`./Reindeer --index -f test/fof_unitigs.txt -g test/output_bcalm/union_graph.unitigs.fa`

# Index query

* Use case 1: queries are in a fasta/fastq file.

Simply provide the fasta query file to Reindeer using `-q`, along with the directory of index files that were generated during index construction (`output_reindeer` by default), using `-l':

`./Reindeer --query -q test/query_test.fa -l output_reindeer`

* Use case 2: you'd like an interactive query, for instance for a web-service.

In this case do not provide the query file, and Reindeer will automatically request one.

Results are written in `out_query_BLight*.out`.

When a query file is provided, * is replaced by 0.

In interactive mode, if n queries are made, * equals 0 to n-1.

# Index and query k-mers abundances

In order to have k-mer abundances instead of presence/absence per indexed dataset, use `--count` option.

`./Reindeer --index -f test/fof.txt --bcalm --count`
` ./Reindeer --query -l output_reindeer -q test/query_test.fa --bcalm --count`



# Output

By default, query outputs are written in `output_reindeer/query_results/`

## k-mer presence/absence

Here is one example line from Reindeer's query:

>    >SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51 dataset1:100% dataset3:71%

A query sequence appears in the output if more than `t`% k-mers were found an indexed dataset.
Headers from the initial files are used to refer to query sequences (for instance `>SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51`).

Then, the indexed datasets it appears in are written, separated by spaces. The numbers (i.e., `1` and `2` in `dataset1:100% dataset3:71%`) refer to the rank (starting at 1) of the datasets in the inital fof file.

In this example, it means than `t`% k-mers from sequence `>SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51` were in the first and third sample written in the fof file.

Finally, we can see that 100% of the k-mers from the queried sequence appear in dataset `1`, and only 71% appear in dataset `3`.

## k-mer abundances

The output changes slightly:

>    >SRR10092187.6 HISEQ:815:HK2NNBCXY:1:1101:1224:2136 length=51 dataset1:66 dataset2:4

This time, instead of percentages of present k-mers, the query counts per datasets are indicated.
In this example, `dataset1:66` means that the queried k-mers correspond to a region which abundance is 66 in dataset 1.

