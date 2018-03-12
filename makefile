CC=g++
CFLAGS=   -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -lz -fopenmp
LDFLAGS=-flto -lpthread -lz -fopenmp



EXEC=blight


all: $(EXEC)

blight: blight.o ksl.o
	$(CC) -o $@ $^ $(LDFLAGS)

blight.o: blight.cpp ksl.h
	$(CC) -o $@ -c $< $(CFLAGS)

ksl.o: ksl.cpp
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
