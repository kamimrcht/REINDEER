#!/bin/bash

mkdir bin;

cd bcalm2;
mkdir build;
cd build;
cmake ..;
make -j4;
mv bcalm ../..bin;

cd ../../src;
make -j4;
mv Reindeer ../bin;

rm -f *.o
