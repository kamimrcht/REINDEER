CC=g++

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CFLAGS+=-O0 -DDEBUG
	WARNS= -Wall -Wextra -Wno-format -Werror=float-equal -Wuseless-cast -Wlogical-op -Wcast-align -Wtrampolines -Werror=enum-compare -Wstrict-aliasing=2 -Werror=parentheses -Wnull-dereference -Werror=restrict -Werror=logical-op -Wsync-nand -Werror=main -Wshift-overflow=2 -Werror=pointer-sign -Wcast-qual -Werror=array-bounds -Werror=char-subscripts -Wshadow -Werror=ignored-qualifiers -Werror=sequence-point -Werror=address -Wduplicated-branches -Wsign-compare -Wodr -Wno-unknown-pragmas -Wnarrowing -Wsuggest-final-methods  -Wformat-signedness -Wrestrict -Werror=aggressive-loop-optimizations -Werror=missing-braces -Werror=uninitialized -Wframe-larger-than=32768 -Werror=nonnull -Wno-unused-function -Werror=init-self -Werror=empty-body -Wdouble-promotion -Wfatal-errors -Werror=old-style-declaration -Wduplicated-cond -Werror=write-strings -Werror=return-type -Werror=volatile-register-var -Wsuggest-final-types -Werror=missing-parameter-type -Werror=implicit-int
	DEBUG_SYMS=1
else
	CFLAGS+=-DNDEBUG -O3 -flto -march=native -mtune=native
	WARNS=-Wfatal-errors
endif

DEBUG_SYMS ?= 1
ifeq ($(DEBUG_SYMS), 1)
	CFLAGS+=-g
endif

CFLAGS+=-std=c++11 -pipe -lz -fopenmp ${WARNS}
INC=blight/blight.h blight/bbhash.h blight/common.h src/utils.hpp
EXEC=Reindeer


all: $(EXEC)

Reindeer: reindeer.o blight.o
	$(CC) -o $@ $^ $(CFLAGS)

reindeer.o: reindeer.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)

blight.o: blight/blight.cpp $(INC)
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
