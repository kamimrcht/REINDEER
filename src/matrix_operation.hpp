

#include "utils.hpp"
#include "../trle/trle.h"

#ifndef MAT
#define MAT

using namespace std;


// dump rle vector on buffer
void dump_compressed_vector_buff(vector<uint16_t>& counts, int64_t monotig_id, string& buffer, unsigned char *in);
// load rle encoded matrix from disk (keep compressed in ram)
vector<unsigned char*> load_compressed_vectors(const string& input_file, vector<unsigned>&vector_sizes, uint64_t& color_number, uint64_t& monotig_number, long eq_class_nb);


void dump_compressed_vector_bucket_disk_query(vector<uint16_t>& counts, int64_t monotig_id, unsigned char *in,vector<ofstream*>& bucket_files,  vector<uint8_t>& colors, bool record_counts);
//~ void dump_compressed_vector_bucket(vector<uint16_t>& counts, int64_t monotig_id, unsigned char *in,  vector<ofstream*>& bucket_files, vector<uint8_t>& colors, bool record_counts );

void dump_compressed_vector_bucket(int64_t monotig_id, vector<ofstream*>& bucket_files, string& header );


void read_matrix_compressed_line(ifstream& in, int64_t& rank, char* comp, unsigned& comp_size);

#endif
