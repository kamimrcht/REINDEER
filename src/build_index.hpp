#include <chrono>
#include "eq_classes.hpp"
#include "../blight/blight.h"
#include "../blight/utils.h"
#include "utils.hpp"
#include "../blight/lz4/lz4_stream.h"
#ifndef BUI
#define BUI


using namespace std;

// read colors from bcalm headers
vector<uint8_t>  get_colors_monotigs(string& line);
// read counts from bcalm headers
vector<uint16_t> get_counts_monotigs(string& line);


// build color/count matrix and dump it


// dispatch count vectors in files. Similar counts go in similar files
void write_matrix_in_bucket_files(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts,  uint k, uint nb_threads,  string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file, bool do_query_on_disk , bool quantize, bool log);


// color using monotig file: either build and dump the color matrix during the index construction, or load it during the query
void do_coloring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, uint k,  uint nb_threads,  string& output, vector<unsigned char*>& compr_monotig_color,vector<unsigned>& compr_monotig_color_size, bool do_query_on_disk,long& eq_class_nb, uint64_t& nb_colors,  bool quantize, bool log);

// load dumped index(+colors)

kmer_Set_Light* load_rle_index(uint k, string& color_load_file, string& color_dump_file, string& fof, bool record_counts, uint nb_threads, string& output, vector<unsigned char*> &compr_monotig_color,  vector<unsigned>& compr_monotig_color_sizes, bool do_query_on_disk, long& eq_class_nb, uint64_t& nb_colors,  bool quantize, bool log);

// build index from new file
void build_index(uint k, uint m1,uint m2,uint m3, uint bit, string& color_load_file, string& color_dump_file, string& fof, bool record_counts,  kmer_Set_Light* ksl, uint nb_threads,  string& output, bool do_query_on_disk, bool quantize, bool do_log, uint64_t nb_colors);


#endif
