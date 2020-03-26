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

## k-mer presence/absence

Here is one example line from Reindeer's query:

>    >SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51 dataset1:100% dataset3:71%

A query sequence appears in the output if more than `t`% k-mers were found an indexed dataset.
Headers from the initial files are used to refer to query sequences (for instance `>SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51`).

Then, the indexed datasets it appears in are written, separated by spaces. The numbers (i.e., `1` and `3` in `dataset1:100% dataset3:71%`) refer to the rank (starting at 1) of the datasets in the inital fof file.

In this example, it means than `t`% k-mers from sequence `>SRR10092187.47 HISEQ:815:HK2NNBCXY:1:1101:1829:2108 length=51` were in the first and third sample written in the fof file.

Finally, we can see that 100% of the k-mers from the queried sequence appear in dataset `1`, and only 71% appear in dataset `3`.

## k-mer abundances

The output changes slightly:

>    >SRR10092187.6 HISEQ:815:HK2NNBCXY:1:1101:1224:2136 length=51 66 2:4

This time, instead of percentages of present k-mers, the query counts per datasets are indicated.
In this example, `66` means that the queried k-mers correspond to a region which abundance is 66 in dataset 1.

In the third column, we observe that two numbers were given, separated by a `:` It means that the sequence spans 2 unitigs in the graph of dataset 2, so REINDEER reported the counts in these 2 unitigs.


# Beta options

## log counts/quantized counts


`./Reindeer --index --log-count -f fof_unitigs.txt`

`./Reindeer --index --quantization -f test/fof_unitigs.txt `


## input paired-end reads (to bcalm)


`./Reindeer --index --paired-end --bcalm  -f fof.txt`

# Reproduce the manuscript's results

We provide a page with scripts to reproduce the results we show in our manuscript, link [here](https://github.com/kamimrcht/REINDEER/tree/master/reproduce_manuscript_results).



