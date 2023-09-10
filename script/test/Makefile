infofile:=reindeer_matrix_eqc_info

SHELL:=/bin/bash

all: init
	../update_reindeer_info -f fof_unitigs.txt -i .
	@diff $(infofile).reindeer1.1 ../../test/index-1.1/$(infofile) > /dev/null ;\
	if [ $$? -eq 0 ]; then echo "Rename old $(info_file) file: OK"; else echo "Error when renaming old binary $(infofile)" ; fi
	@diff $(infofile) ../../test/index-1.4/$(infofile) > /dev/null ;\
	if [ $$? -eq 0 ]; then echo "Convert is OK"; make clean; else echo "Error when converting $(infofile)"; diff $(infofile) ../../test/index-1.4/$(infofile) ; fi

init: clean
	ln -s ../../test/fof_unitigs.txt
	ln -s ../../test/output_bcalm
	cp ../../test/index-1.1/$(infofile) .


clean:
	rm -f fof_unitigs.txt output_bcalm $(infofile)*

