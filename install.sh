#!/bin/bash

rm -rf Reindeer bin;
mkdir bin;

( cd bcalm2 && rm -rf build && mkdir build &&  cd build && cmake .. && make -j4 )
mv bcalm2/build/bcalm bin;
make -j4;

rm -f *.o;
