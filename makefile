CC=g++
CFLAGS=   -Ofast -std=c++11 -flto   -pipe -funit-at-a-time  -Wfatal-errors -lz -fopenmp
LDFLAGS= -lpthread -lz -flto -fopenmp


EXEC=bench_blight Colored_De_Bruijn_graph_snippet Abundance_De_Bruijn_graph_snippet


all: $(EXEC)

Colored_De_Bruijn_graph_snippet: Colored_De_Bruijn_graph_snippet.o blight.o mmh.o
	$(CC) -o $@ $^ $(LDFLAGS)

Abundance_De_Bruijn_graph_snippet: Abundance_De_Bruijn_graph_snippet.o blight.o mmh.o
	$(CC) -o $@ $^ $(LDFLAGS)

bench_blight: bench_blight.o blight.o mmh.o
	$(CC) -o $@ $^ $(LDFLAGS)

bench_blight.o: bench_blight.cpp blight.h
	$(CC) -o $@ -c $< $(CFLAGS)

Colored_De_Bruijn_graph_snippet.o: Colored_De_Bruijn_graph_snippet.cpp blight.h
	$(CC) -o $@ -c $< $(CFLAGS)

Abundance_De_Bruijn_graph_snippet.o: Abundance_De_Bruijn_graph_snippet.cpp blight.h
	$(CC) -o $@ -c $< $(CFLAGS)

mmh.o: MurmurHash3.cpp MurmurHash3.h
	$(CC) -o $@ -c $< $(CFLAGS)

blight.o: blight.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
