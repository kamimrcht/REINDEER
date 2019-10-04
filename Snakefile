#rule target:
#	input:
#		"main_graph"
	
	
rule create_list_graph:
	output: 
		"list_graph"
	shell:
		"""
		rm -f list_graph
		> {output}
		"""

#rule bct_each_graph:
#	input:
#		"/home/camillemarchet/projects/BLIGHT/test_mantis/input_files_all/{sample}.fastq"
#	output:
#		"bct_{sample}"
#	shell:
#		"""
#		mkdir -p {output}
#		~/softwares/BCT/Bct.py -u {input} -t 4 -o {output}
#		#~/softwares/bcalm/build/bcalm -in {input} -kmer-size 31 -abundance-min 2 -nb-cores 4 -out {output}
#		rm -f *.h5 *glue* {output}/.dbg*
#		echo "{output}/dbg31.fa" >> list_graph
#		"""



rule bcalm_single_graphs:
	input:
		"/home/camillemarchet/projects/BLIGHT/{sample}.fastq"
	output:
		"bcalm_{sample}"
	shell:
		mkdir -p {output}
		~/softwares/BCT/Bct.py -u {input} -t 4 -o {output}
		~/softwares/bcalm/build/bcalm -in {input} -kmer-size 31 -abundance-min 2 -nb-cores 4 -out {output}
		rm -f *.h5 *glue* {output}/.dbg*
		echo "{output}/dbg31.fa" >> list_graph

rule bcalm_union_graph:
	input:
		"list_graph"
	output:
		"main_graph"
	shell:
		"""
		~/softwares/bcalm/build/bcalm -in {input} -nb-cores 4 -abundance-min 1 -out {output}
		rm -f *.h5 *glue*
		"""
		
rule reindeer:
	input:
		"main_graph"
	shell:
		"""
		/home/camillemarchet/projects/reindeer/Reindeer -o /home/camillemarchet/projects/reindeer/test/samples.list -g {input}.unitigs.fa -t 4 -q /home/camillemarchet/projects/reindeer/test/queries.fa -k 31 > log_reindeer ;
		"""




