#!/bin/bash

mkdir bin;

cd bcalm2;
mkdir build;
cd build;
cmake ..;
make -j4;
mv bcalm ../..bin;

cd ../..;
make -j4;
mv Reindeer bin;
