CXX=g++
CC=gcc

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS+=-O0 -DDEBUGi -fstrict-aliasing
	WARNS= -Wall -Wextra -Wno-format -Werror=float-equal -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 -Werror=pointer-sign -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wno-unknown-pragmas -Wnarrowing -Wsuggest-final-methods  -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion -Wfatal-errors -Werror=old-style-declaration -Wduplicated-cond -Werror=write-strings -Werror=return-type -Werror=volatile-register-var -Wsuggest-final-types -Werror=missing-parameter-type -Werror=implicit-int
	DEBUG_SYMS=1
else
	CFLAGS+=-DNDEBUG -Ofast -flto -march=native -mtune=native -fstrict-aliasing 
	CFLAGS2+=-DNDEBUG -Ofast -flto -march=native -mtune=native -fstrict-aliasing  
	#WARNS=-Wfatal-errors
endif

DEBUG_SYMS ?= 1
ifeq ($(DEBUG_SYMS), 1)
	CFLAGS+=-g
endif
CFLAGS2+= -w -Wall -std=gnu99 -DUSE_THREADS  -fstrict-aliasing -Iext $(DEFS)
CFLAGS+=-std=c++11 -pipe -lz -fopenmp ${WARNS}
INC=blight/blight.h blight/bbhash.h blight/common.h src/utils.hpp src/reindeer.hpp src/launch_bcalm.hpp trle/trle.h trle/trle_.h trle/conf.h src/eq_classes.hpp src/query.hpp src/build_index.hpp blight/utils.h src/eq_classes.hpp blight/zstr.hpp blight/common.h blight/robin_hood.h src/monotig.hpp src/matrix_operation.hpp
EXEC=Reindeer


all: $(EXEC)

Reindeer: main.o blight.o utils_b.o trlec.o trled.o monotig.o utils.o reindeer.o query.o build_index.o eq_classes.o launch_bcalm.o matrix_operation.o
	$(CXX) -o $@ $^ $(CFLAGS)

main.o: main.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

monotig.o: src/monotig.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

blight.o: blight/blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)
	
utils_b.o: blight/utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

trlec.o: trle/trlec.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

trled.o: trle/trled.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

utils.o: src/utils.cpp  $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)
	
reindeer.o: src/reindeer.cpp  $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)
	
query.o: src/query.cpp  $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)
	
eq_classes.o: src/eq_classes.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)
	
build_index.o: src/build_index.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)
	
launch_bcalm.o: src/launch_bcalm.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

matrix_operation.o: src/matrix_operation.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf trlec.o trled.o utils.o main.o blight.o utils_b.o reindeer.o query.o build_index.o eq_classes.o launch_bcalm.o monotig.o matrix_operation.o
	rm -rf $(EXEC)

rebuild: clean $(EXEC)
