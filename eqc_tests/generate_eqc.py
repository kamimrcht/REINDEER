import math,itertools
import sys

n_eqc = int(sys.argv[1])

l = math.ceil(math.log(n_eqc,2))
print("l=",l)

# open datasets files
of = []
nb_kmers = []
filename_list = str(n_eqc) + "_eqc_list.txt"
filename_list_f = open(filename_list,"w")
for i in range(l):
	filename = str(n_eqc) + "_eqc_datasets/dataset%d.fa" % i
	filename_list_f.write(filename+"\n")
	of += [open(filename,"w")]
	nb_kmers += [0]
filename_list_f.close()
	

nb_eqc = 0
for eqc in itertools.product([0,1],repeat=l):
	nb_eqc += 1
	#print(eqc)
	#generate a deterministic kmer
	kmer = ""
	for c in eqc:
		kmer += "A" if c == 0 else "G"
	# place that kmer in a dataset
	for i,c in enumerate(eqc):
		if c == 1:
			of[i].write(">\n"+kmer+"\n")
			nb_kmers[i] += 1

print("generated",nb_eqc,"equivalence classes")

for i in range(l):
	print("write",nb_kmers[i],"kmers to dataset",i)
