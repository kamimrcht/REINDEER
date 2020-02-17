CXX=g++
CC=gcc

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS+=-O0 -DDEBUG
	WARNS= -Wall -Wextra -Wno-format -Werror=float-equal -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 -Werror=pointer-sign -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wno-unknown-pragmas -Wnarrowing -Wsuggest-final-methods  -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion -Wfatal-errors -Werror=old-style-declaration -Wduplicated-cond -Werror=write-strings -Werror=return-type -Werror=volatile-register-var -Wsuggest-final-types -Werror=missing-parameter-type -Werror=implicit-int -g
	DEBUG_SYMS=1
else
	CFLAGS+=-DNDEBUG -Ofast -flto -march=native -mtune=native -g
	WARNS=
endif

CFLAGS2+= -w -Wall -std=gnu99 -DUSE_THREADS  -fstrict-aliasing -Iext $(DEFS)
CFLAGS+=-std=c++11 -pipe -lz -fopenmp -msse4 ${WARNS}
INC=blight.h bbhash.h common.h
EXEC= bench_minitig bench_blight


all: $(EXEC)

bench_blight: bench_blight.o blight.o utils.o
	$(CXX) -o $@ $^ $(CFLAGS)

bench_blight.o: bench_blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

bench_minitig: bench_minitig.o blight.h minitig.o    blight.o    utils.o trlec.o
	$(CXX) -o $@ $^ $(CFLAGS)

bench_minitig.o: bench_minitig.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

split: split.o blight.o
	$(CXX) -o $@ $^ $(CFLAGS)

split.o: split.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

blight.o: blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

minitig.o: minitig.cpp  $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

trlec.o: trlec.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

trled.o: trled.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
