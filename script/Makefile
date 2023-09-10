SRC:=update_reindeer_info.cpp
BIN:=$(SRC:.cpp=)
SHELL:=/bin/bash

all: $(BIN)

$(BIN): $(SRC)

.PHONY: test

# do test
infofile:=reindeer_matrix_eqc_info
test: init
	../update_info -f fof_unitigs.txt -i index-1.1
	@diff index-1.1/$(infofile).reindeer1.1 ../test/index-1.1/$(infofile) > /dev/null ;\
	if [ $$? -eq 0 ]; then echo "Rename old $(info_file) file: OK"; else echo "Error when renaming old binary $(infofile)" ; fi
	@diff index-1.1/$(infofile) index-1.4/$(infofile) > /dev/null ;\
	if [ $$? -eq 0 ]; then echo "Convert is OK"; make clean; else echo "Error when converting $(infofile)"; diff index-1.1/$(infofile) index-1.4/$(infofile) ; fi

init: clean
	ln -s ../test/fof_unitigs.txt
	ln -s ../test/output_bcalm
	../Reindeer --index -f fof_unitigs.txt -o index-1.4 > log 2>&1
	cp -r index-1.4 index-1.1
	cp ../test/index-1.1/$(infofile) index-1.1/


clean:
	rm -f fof_unitigs.txt output_bcalm log $(infofile)*
	rm -rf index-1.1

