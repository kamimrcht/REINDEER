mkdir out_test;
./Reindeer --count --index -f test/fof_unitigs.txt -o out_test > out_test/log 2>&1;
eqc=$(grep "Number of equivalence classes found" out_test/log | awk '{print $6}');
if (( eqc==180 )); then echo "Test passed." ; fi;
rm -r out_test;
