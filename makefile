CC=g++

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS+=-O0 -DDEBUG
	WARNS= -Wall -Wextra -Wno-format -Werror=float-equal -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 -Werror=pointer-sign -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wno-unknown-pragmas -Wnarrowing -Wsuggest-final-methods  -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion -Wfatal-errors -Werror=old-style-declaration -Wduplicated-cond -Werror=write-strings -Werror=return-type -Werror=volatile-register-var -Wsuggest-final-types -Werror=missing-parameter-type -Werror=implicit-int
	DEBUG_SYMS=1
else
	CFLAGS+=-DNDEBUG -O3 -flto -march=native -mtune=native
	WARNS=
endif

DEBUG_SYMS ?= 1
ifeq ($(DEBUG_SYMS), 1)
	CFLAGS+=-g
endif

CFLAGS+=-std=c++11 -pipe -lz -fopenmp ${WARNS}
INC=blight.h bbhash.h common.h
EXEC=bench_blight split


all: $(EXEC)

Colored_De_Bruijn_graph_snippet: Colored_De_Bruijn_graph_snippet.o blight.o
	$(CC) -o $@ $^ $(CFLAGS)

interactive: interactive.o blight.o
	$(CC) -o $@ $^ $(CFLAGS)

Abundance_De_Bruijn_graph_snippet: Abundance_De_Bruijn_graph_snippet.o blight.o
	$(CC) -o $@ $^ $(CFLAGS)

bench_blight: bench_blight.o blight.o
	$(CC) -o $@ $^ $(CFLAGS)

bench_blight.o: bench_blight.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

split: split.o blight.o
	$(CC) -o $@ $^ $(CFLAGS)

split.o: split.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

Colored_De_Bruijn_graph_snippet.o: Colored_De_Bruijn_graph_snippet.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

interactive.o: interactive.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

Abundance_De_Bruijn_graph_snippet.o: Abundance_De_Bruijn_graph_snippet.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

blight.o: blight.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
