SRC:=update_reindeer_info.cpp
BIN:=$(SRC:.cpp=)
SHELL:=/bin/bash

all: $(BIN)

$(BIN): $(SRC)

.PHONY: test

# do test
infofile:=reindeer_matrix_eqc_info
test: init
	for v in 1.1 1.0.2; do for mem in nodiskquery diskquery; do \
	if [ -d index-$$v-$$mem ]; then ./update_reindeer_info -f fof_unitigs.txt -i index-$$v-$$mem ;\
	diff index-$$v-$$mem/$(infofile).txt index-1.4-$$mem/$(infofile).txt > /dev/null ;\
	if [ $$? -eq 0 ]; then echo "Rename old $(info_file) file in index-$$v-$$mem: OK"; else echo "Error when renaming old binary index-$$v-$$mem/$(infofile)" ; fi ;\
	fi; done; done

init: clean
	ln -s ../test/fof_unitigs.txt
	ln -s ../test/output_bcalm
	../Reindeer --index -f fof_unitigs.txt --mem-query -o index-1.4-nodiskquery > log 2>&1
	../Reindeer --index -f fof_unitigs.txt -o index-1.4-diskquery > log 2>&1
	[ -s ../Reindeer-1.1 ] && ../Reindeer-1.1 --index -f fof_unitigs.txt -o index-1.1-nodiskquery > log 2>&1
	#[ -s ../Reindeer-1.1 ] && ../Reindeer-1.1 --index -f fof_unitigs.txt --disk_query -o index-1.1-diskquery > log 2>&1
	sif=/data/containers/reindeer.sif [ -e $$sif ] && singularity exec -B /scratch -B /data $$sif Reindeer --index -f fof_unitigs.txt -o index-1.0.2-nodiskquery > log 2>&1 || echo "no version 1.0.2"
	sif=/data/containers/reindeer.sif [ -e $$sif ] && singularity exec -B /scratch -B /data $$sif Reindeer --index -f fof_unitigs.txt --disk-query -o index-1.0.2-diskquery > log 2>&1 || echo "no version 1.0.2"
	#cp -r index-1.4 index-1.1
	#cp ../test/index-1.1/$(infofile) index-1.1/


clean:
	rm -f fof_unitigs.txt output_bcalm log $(infofile)*
	rm -rf index-1.1-* index-1.4* index-1.0.2*

