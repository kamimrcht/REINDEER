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
CFLAGS+=-std=c++17 -pipe -lz -fopenmp -Iblight/lz4 ${WARNS}
INC=blight/blight.h blight/bbhash.h blight/common.h src/utils.hpp src/reindeer.hpp src/launch_bcalm.hpp trle/trle.h trle/trle_.h trle/conf.h src/eq_classes.hpp src/query.hpp src/build_index.hpp blight/utils.h src/eq_classes.hpp blight/zstr.hpp blight/common.h blight/robin_hood.h src/monotig.hpp src/matrix_operation.hpp blight/lz4/lz4_stream.h
EXEC=Reindeer reindeer_socket update_info
LZ4H=blight/lz4/lz4frame.o blight/lz4/lz4.o blight/lz4/xxhash.o blight/lz4/lz4hc.o

VERSION := $(shell git describe --tags --always)

all: $(EXEC)

Reindeer: main.o blight.o utils_b.o trlec.o trled.o monotig.o utils.o reindeer.o query.o build_index.o eq_classes.o launch_bcalm.o matrix_operation.o $(LZ4H)
	$(CXX) -o $@ $^ $(CFLAGS)

main.o: main.cpp $(INC)
	echo "#define VERSION \"$(VERSION)\"" > version.h
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

blight/lz4/lz4frame.o: blight/lz4/lz4frame.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

blight/lz4/lz4.o: blight/lz4/lz4.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

blight/lz4/lz4hc.o: blight/lz4/lz4hc.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)

blight/lz4/xxhash.o: blight/lz4/xxhash.c $(INC)
	$(CC) -o $@ -c $< $(CFLAGS2)


utils.o: src/utils.cpp  $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

reindeer.o: src/reindeer.cpp  $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

query.o: src/query.cpp  $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

eq_classes.o: src/eq_classes.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

build_index.o: src/build_index.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

launch_bcalm.o: src/launch_bcalm.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)

matrix_operation.o: src/matrix_operation.cpp $(INC)
	$(CXX) -o $@ -c $< $(CFLAGS)


reindeer_socket: reindeer_socket.o blight.o utils_b.o trlec.o trled.o utils.o reindeer.o query.o build_index.o eq_classes.o matrix_operation.o $(LZ4H)
	$(CXX) -o $@ $^ $(CFLAGS)

reindeer_socket.o: src/reindeer_socket.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

update_info: update_reindeer_info.o
	$(CXX) -o $@ $^ $(CFLAGS)

update_reindeer_info.o: script/update_reindeer_info.cpp
	$(CXX) -o $@ -c $< $(CFLAGS)

# do test
test:
	cd test && $(MAKE) all

.PHONY: test

clean:
	rm -f *.o $(EXEC)

rebuild: clean $(EXEC)
