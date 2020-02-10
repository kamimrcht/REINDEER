# some code to determine whether Assumption 1 from the paper holds
# requires a KMC-counted dataset
# Rayan, Jan 2020

unitigs_file = "/opt/reindeer/datasets/reads/SRR1002076.unitigs.fa"
kmc_counts = "/opt/reindeer/datasets/reads/SRR1002076.countsk21.raw"
k=21

#sample_kmc=10000000
#sample_points=100000
sample_kmc=0
sample_points=0

tab = str.maketrans("ACGT","TGCA")
def revcomp(s):
    return s.translate(tab)[::-1]

def normalize(kmer):
    return min(kmer,revcomp(kmer))

counts_dict = dict()
for i,line in enumerate(open(kmc_counts)):
    kmer, count = line.strip().split()
    counts_dict[normalize(kmer)] = int(count)
    
    if sample_kmc != 0 and i > sample_kmc: break

print("loaded counts")

from Bio import SeqIO

points = []
fasta_sequences = SeqIO.parse(open(unitigs_file),'fasta')
for j,fasta in enumerate(fasta_sequences):
    name, sequence = fasta.id, str(fasta.seq)
    for field in fasta.description.split(' '):
        if field.startswith("km"):
            abundance=float(field[5:])
            break
    #print(fasta.description,abundance)

    for i in range(len(sequence)-k+1):
        kmer = normalize(sequence[i:i+k])
        if kmer not in counts_dict:
            #print("kmer not found",kmer)
            continue
        count = counts_dict[kmer]
        points += [(abundance,count)]    
    
    if sample_points != 0 and len(points)> sample_points: break

print(points)

import scipy.stats
print("spearman",scipy.stats.spearmanr(points))
print("pearson",scipy.stats.pearsonr(list(map(lambda x:x[0],points)),list(map(lambda x:x[1],points))))
