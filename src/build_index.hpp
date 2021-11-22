#include <chrono>
#include "eq_classes.hpp"
#include "../blight/blight.h"
#include "../blight/utils.h"
#include "utils.hpp"
#include "../blight/lz4/lz4_stream.h"
#include "reindeer_index.hpp"
#ifndef BUI
#define BUI


using namespace std;

// read colors from bcalm headers
vector<uint8_t>  get_colors_monotigs(string& line);
// read counts from bcalm headers
vector<uint16_t> get_counts_monotigs(string& line);


// build color/count matrix and dump it


// dispatch count vectors in files. Similar counts go in similar files
void write_matrix_in_bucket_files(reindeer_index& index_values, kmer_Set_Light* ksl, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size);


// color using monotig file: either build and dump the color matrix during the index construction, or load it during the query
void do_coloring(reindeer_index& index_values, kmer_Set_Light* ksl, vector<unsigned char*>& compr_monotig_color,vector<unsigned>& compr_monotig_color_size);

// load dumped index(+colors)

kmer_Set_Light* load_rle_index(reindeer_index& index_values, vector<unsigned char*> &compr_monotig_color,  vector<unsigned>& compr_monotig_color_sizes);

// build index from new file
void build_index(reindeer_index& index_values, kmer_Set_Light* ksl);


#endif
