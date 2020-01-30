# Reindeer index construction - record counts
./Reindeer --index --count -f fof_reindeer.lst -k 21 -o out_reindeer_2585_counts
# Reindeer index construction - record colors
./Reindeer --index -f fof_reindeer.lst -k 21 -o out_reindeer_2585_colors
# Reindeer index construction - record quantized counts
./Reindeer --index --count --quantization -f fof_reindeer.lst -k 21 -o out_reindeer_2585_colors
# Reindeer index construction - record log counts
./Reindeer --index --count --log-count -f fof_reindeer.lst -k 21 -o out_reindeer_2585_colors
