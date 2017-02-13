CC=g++
CFLAGS=-Wall  -Ofast -std=c++11  -flto -pipe -funit-at-a-time  -Wfatal-errors -fopenmp
LDFLAGS=-flto -lpthread -fopenmp



EXEC=blight

all: $(EXEC)



blight: blight.cpp
	$(CC) -o $@  $^ $(CFLAGS)

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CFLAGS)



clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)

