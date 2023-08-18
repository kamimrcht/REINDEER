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

all:

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
