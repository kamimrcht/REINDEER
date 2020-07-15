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
INC=../blight/blight.h ../blight/bbhash.h ../blight/common.h ../trle/trle.h ../trle/trle_.h ../trle/conf.h ../blight/utils.h ../blight/zstr.hpp ../blight/common.h ../blight/robin_hood.h ../src/monotig.hpp ../blight/lz4/lz4_stream.h ../src/utils.hpp
EXEC=bench_monotig
LZ4H=lz4frame.o lz4.o xxhash.o lz4hc.o


all: $(EXEC)

bench_monotig: bench_monotig.o blight.o utils_b.o trlec.o trled.o monotig.o utils.o $(LZ4H)
	$(CXX) -o $@ $^ $(CFLAGS)

utils.o: ../src/utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

bench_monotig.o: bench_monotig.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

monotig.o: ../src/monotig.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

blight.o: ../blight/blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)
	
utils_b.o: ../blight/utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

trlec.o: ../trle/trlec.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

trled.o: ../trle/trled.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

lz4frame.o: ../blight/lz4/lz4frame.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

lz4.o: ../blight/lz4/lz4.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

lz4hc.o: ../blight/lz4/lz4hc.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)
	
xxhash.o: ../blight/lz4/xxhash.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2) 

clean:
	rm -rf trlec.o trled.o utils.o main.o blight.o utils_b.o reindeer.o query.o build_index.o eq_classes.o launch_bcalm.o monotig.o matrix_operation.o bench_monotig.o lz4*.o xxhash.o
	rm -rf $(EXEC)

rebuild: clean $(EXEC)
