mkdir out_test;
echo "Building REINDEER index ....................................................."
./Reindeer --index -f test/fof_unitigs.txt -o out_test > out_test/log 2>&1;
eqc=$(grep "Number of equivalence classes found" out_test/log | awk '{print $6}');
if (( eqc==180 )); then echo "    -> Construction test passed." ; else echo "    -> Construction test error."; fi;

echo "Querying index..............................................................."
./Reindeer --query -l out_test -q  test/output_bcalm/SRR10092187_10k.unitigs.fa -o out_test >> out_test/log 2>&1;
DIFF=$(diff test/out_test/out_query_Reindeer0.out test/query_results/out_query_Reindeer0.out)
if [ "$DIFF" != "" ] ; then echo "    -> Query test error."; else echo "    -> Query test passed."; fi;

rm -r out_test;
