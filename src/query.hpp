

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <unordered_set>
#include <set>
#include <mutex>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include "../trle/trle.h"
#include "../blight/blight.h"
#include "utils.hpp"

#ifndef QQ
#define QQ

using namespace std;
using namespace chrono;




// decode rl encoded vector (color/count of ONE minitig)
unsigned char* decode_vector(unsigned char* minitig_counts, unsigned vector_size, uint16_t color_number, bool record_counts);

// convert char [] counts/colors to uint
vector<uint16_t> count_string_to_count_vector(unsigned char* count_char_minitig, unsigned size);

// convert char [] counts/colors to uint
vector<uint8_t> count_string_to_count_vector8(unsigned char* count_char_minitig, unsigned size);

// for ONE minitig, get a vector of its counts/colors (uint) from the rle encoded matrix
vector<uint16_t> get_count_minitig(unsigned char* minitig_counts, unsigned vector_size, uint16_t color_number, bool record_counts);



void get_position_vector_query_disk(vector<long>& position_in_file, string& position_file_name, uint64_t nb_minitig);
void get_matrix_line_query_disk(int64_t rank, unsigned char* color, unsigned& line_size, long position_in_matrix, ifstream& in );


long get_matrix_line_query(int64_t	rank, unsigned char* color, unsigned& line_size, vector<long>& position_in_file, vector<unsigned char*>& compr_minitig_color );

void get_colors_counts_query_eq_classes(vector<int64_t>& kmer_ids,   uint16_t color_number, vector<vector<uint16_t>>& query_counts,  vector<unsigned char*>& compr_minitig_color, vector<unsigned>&compr_minitig_color_size,vector<long>& position_in_file, bool record_counts, vector<vector<uint8_t>>& query_colors, bool do_query_on_disk, string& rd_file);

// for all queried k-mers, get the colors/counts in vector<vector<uint16_t>>& query_counts
void get_colors_counts(vector<int64_t>& kmer_ids, bool record_counts, uint16_t color_number, vector<int64_t>& kmers_colors, vector<vector<uint16_t>>& query_counts, vector<unsigned char*>& compr_minitig_color, vector<unsigned>&compr_minitig_color_size, vector<long>& position_in_file);

// compute a string that sums up the count(s) for each dataset
void write_count_output(bool record_counts, vector<vector<uint16_t>>& query_counts, uint16_t color_number , vector<string>& toW, vector<string>& color_counts);

void write_results_above_threshold(string& toWrite, vector<vector<uint16_t>>& query_counts, uint16_t color_number , vector<string>& toW, vector<string>& color_counts,  string& header, bool record_counts, uint threshold, string& line, uint k);

void write_output(vector<int64_t>& kmers_colors, string& toWrite, bool record_reads, bool record_counts, vector<vector<uint32_t>>& query_unitigID,vector<vector<uint32_t>>& query_unitigID_tmp,  uint16_t& color_number, string& header, string& line, uint k, uint threshold,  vector<vector<uint16_t>>& query_counts, vector<vector<uint8_t>>& query_colors);




void doQuery(string& input, string& name, kmer_Set_Light& ksl, uint64_t& color_number, uint k, bool record_counts, bool record_reads, uint threshold,vector<vector<uint32_t>>& query_unitigID, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color, vector<unsigned >& compr_minitig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb, uint64_t nb_minitig);

void query_by_file(uint& counter, string& entry, kmer_Set_Light& ksl, uint64_t& color_number,   uint k, bool record_counts, bool record_reads, uint threshold, vector<string>& bgreat_files, string& output, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color, vector<unsigned >& compr_minitig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb, uint64_t	 nb_minitig);


void perform_query(kmer_Set_Light& ksl, uint16_t& color_number,  uint k, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query, string& output, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color,  vector<unsigned>& compr_minitig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb);

#endif
