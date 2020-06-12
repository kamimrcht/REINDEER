#!/bin/bash
set -e

#n=1000000000
#n=100000000
n=10000000

function log2 {
    local x=0
    for (( y=$1-1 ; $y > 0; y >>= 1 )) ; do
        let x=$x+1
    done
    echo $x
}
k=$(log2 n)
mkdir -p "$n"_eqc_datasets
python generate_eqc.py $n
rm -Rf "$n"_eqc_index/ && ../Reindeer  --index -f "$n"_eqc_list.txt -k $k -o "$n"_eqc_index -t 20
