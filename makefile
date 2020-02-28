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
INC=blight/blight.h blight/bbhash.h blight/common.h src/utils.hpp src/reindeer.hpp src/launch_bcalm.hpp trle/trle.h trle/trle_.h trle/conf.h src/eq_classes.hpp
EXEC=Reindeer


all: $(EXEC)

Reindeer: main.o blight.o utils.o trlec.o trled.o minitig.o
	$(CXX) -o $@ $^ $(CFLAGS)

main.o: main.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

minitig.o: src/minitig.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

blight.o: blight/blight.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)
	
utils.o: blight/utils.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

trlec.o: trle/trlec.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

trled.o: trle/trled.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

clean:
	rm -rf trlec.o trled.o utils.o main.o blight.o
	rm -rf $(EXEC)

rebuild: clean $(EXEC)
