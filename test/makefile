CXX=g++
CC=gcc

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS+=-O0 -DDEBUGi -fstrict-aliasing -msse4
	#WARNS= -Wall -Wextra -Wno-format -Werror=float-equal -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 -Werror=pointer-sign -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wno-unknown-pragmas -Wnarrowing -Wsuggest-final-methods  -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion -Wfatal-errors -Werror=old-style-declaration -Wduplicated-cond -Werror=write-strings -Werror=return-type -Werror=volatile-register-var -Wsuggest-final-types -Werror=missing-parameter-type -Werror=implicit-int
	DEBUG_SYMS=1
else
	CFLAGS+=-DNDEBUG -Ofast -flto -msse4 -fstrict-aliasing  -g
	CFLAGS2+=-DNDEBUG -Ofast -flto -msse4 -fstrict-aliasing -g
	#WARNS=-Wfatal-errors
endif

DEBUG_SYMS ?= 1
ifeq ($(DEBUG_SYMS), 1)
	CFLAGS+=-g
endif
CFLAGS2+= -w -Wall -std=gnu99 -DUSE_THREADS  -fstrict-aliasing -Iext $(DEFS)
CFLAGS+=-std=c++11 -pipe -lz -fopenmp -Iblight/lz4 ${WARNS}

EXEC=bench_monotig
SRC = bench_monotig.cpp ../src/utils.cpp ../src/monotig.cpp ../blight/blight.cpp ../blight/utils.cpp ../trle/trlec.c ../trle/trled.c ../blight/lz4/lz4frame.c ../blight/lz4/lz4.c ../blight/lz4/lz4hc.c ../blight/lz4/xxhash.c
OBJ = $(addsuffix .o, $(basename $(SRC)))

# Do test by default
# use benchmark to run "bench_monotig"

TESTDIR:=out_test_tmp
QUERY:=query_test

all: test_build test_query_def test_query_format
	@rm -fr $(TESTDIR);

test_build:
	@mkdir -p $(TESTDIR);
	@echo "$(YELLFG)Building REINDEER index $(NoColor)....................................................."
	../Reindeer --index -f fof_unitigs.txt -o $(TESTDIR) > $(TESTDIR)/log 2>&1;
	@eqc=$$(grep "Number of equivalence classes found" $(TESTDIR)/log | awk '{print $6}'); \
	if (( eqc==180 )); then echo "    -> $(GREEN)Construction test passed.$(NoColor)" ; else echo "    -> $(RED)Construction test error.$(NoColor)"; fi;

test_query_def:
	@echo "$(YELLFG)Querying index format; default $(NoColor).............................................."
	../Reindeer --query -l $(TESTDIR) -q $(QUERY).fa -o $(TESTDIR) >> $(TESTDIR)/log 2>&1;
	@if diff $(TESTDIR)/out_query_Reindeer_P40_$(QUERY)_0.out query_results/out_query_Reindeer_P40_$(QUERY)_0.out > /dev/null; \
		then echo "    -> $(GREEN)Query test passed.$(NoColor)";\
		else echo "    -> $(RED)Query test error.$(NoColor)"; fi;

test_query_format:
	@for f in sum average normalize; do \
	echo "$(YELLFG)Querying index format $$f $(NoColor)..................................................." ;\
	../Reindeer --query -l $(TESTDIR) -q $(QUERY).fa -o $(TESTDIR) --format $$f >> $(TESTDIR)/log 2>&1 ;\
	if diff $(TESTDIR)/out_query_Reindeer_P40_$(QUERY)_0.out query_results/out_query_Reindeer_P40_$(QUERY)_0-format_$$f.out > /dev/null; \
		then echo "    -> $(GREEN)Query test passed.$(NoColor)";\
		else echo "    -> $(RED)Query test error.$(NoColor)"; fi; done

benchmark: $(EXEC)

bench_monotig: $(OBJ)
	$(CXX) -o $@ $^ $(CFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS2)

.PHONY:= clean rebuild

clean:
	rm -rf $(OBJ) $(EXEC)

rebuild: clean $(EXEC)
# for color
# color for bash
RED:=\033[0;31m
GREEN:=\033[0;32m
NoColor:=\033[0m
YELLFG:=\033[1;43m
