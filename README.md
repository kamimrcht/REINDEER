# Blight

## de Bruijn graph index with light memory usage

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)


## Description:
Blight is an index structure able to index the kmers of a de Bruijn graph and to associate an unique identifier to each kmer.
The index can be used as a set that will return true when a query kmer is in the set false otherwise.
The index can also be used as an associative data structure that will return a unique identifier (between 1 and N where N is the number of indexed kmer) to each kmer in the set and -1 to kmer that are not indexed.
This index can therefore be used to associate information to kmers in a efficient and exact way with the use of a simple array of value.

## Usage:

./blight -g DBG_unitigs.fa -u query_sequences.fa -k 31

## Mandatoty options:

### -g graph file

The de Bruijn graph to index

Encoded in a fasta file where each represented a simple path of the graph

We recommend the use of BCALM2(https://github.com/GATB/bcalm) to construct such graph from a set of sequences

### -q query file

The sequences to query

Encoded in a fasta file

### -k k value used for graph


The kmer size used


## Performance options:

#### -m minimizer size

Size of minimizer used

Default value 7

Higher size usually lead to lower memory usage (below 12)

#### -n mphf number
Create 4^n mphf

Default value 5 (1024 MPHF)

More mean slower construction but better index

Must be <=m

High value lead to higher memory usage (>7)

#### -s file number
To use 4^s files

Default value 3 (64 files)

Higher mean lower memory usage during construction but more temporary files used

Must be <=n

#### -t core used
To use t threads for faster construction and query

Default value 1


#### -b bit saved per kmer
Will reduce the memory usage of b bit per kmer but query have to check 2^b kmers

Default value 6













