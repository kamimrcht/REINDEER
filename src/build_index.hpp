#include <chrono>
#include "eq_classes.hpp"
#include "../blight/blight.h"
#include "../blight/utils.h"
#include "utils.hpp"
#ifndef BUI
#define BUI


using namespace std;

// read colors from bcalm headers
vector<uint8_t>  get_colors_minitigs(string& line);
// read counts from bcalm headers
vector<uint16_t> get_counts_minitigs(string& line);


// build color/count matrix and dump it
void build_matrix(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file);


// dispatch count vectors in files. Similar counts go in similar files
void write_matrix_in_bucket_files(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file, bool do_query_on_disk );


// color using minitig file: either build and dump the color matrix during the index construction, or load it during the query
void do_coloring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*>& compr_minitig_color,vector<unsigned>& compr_minitig_color_size, bool do_query_on_disk, string& nb_eq_class_file, long& eq_class_nb);
// load dumped index(+colors)
kmer_Set_Light* load_rle_index(uint k, string& color_load_file, string& color_dump_file, string& fof, bool record_counts, bool record_reads, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*> &compr_minitig_color,  vector<unsigned>& compr_minitig_color_sizes, bool do_query_on_disk, string& nb_eq_class_file, long& eq_class_nb);


// build index from new file
void build_index(uint k, uint m1,uint m2,uint m3, uint c, uint bit, string& color_load_file, string& color_dump_file, string& fof, bool record_counts, bool record_reads, uint64_t& color_number, kmer_Set_Light& ksl, uint nb_threads, bool exact, string& output, bool do_query_on_disk, bool quantize, bool do_log);



#endif
