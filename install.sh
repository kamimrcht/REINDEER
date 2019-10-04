#!/bin/bash

mkdir bin;

cd ../bcalm2;
mkdir build;
cd build;
cmake ..;
make;
mv bcalm ../..bin;

cd ../..;
make;
mv Reindeer bin;
