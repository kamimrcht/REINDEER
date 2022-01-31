

#include "../trle/trle.h"
#include "../blight/utils.h"
#include "../blight/robin_hood.h"
#include "matrix_operation.hpp"
#include "utils.hpp"

#ifndef EQ
#define EQ

using namespace std;


struct count_vector
{
	unsigned compressed_size;
	int64_t monotig_rank;
	string compressed;
};


struct compare_vec
{
	bool operator()(count_vector& vec1, count_vector& vec2)
	{
		return vec1.compressed < vec2.compressed;
	}
};


void sort_vectors(vector<count_vector>& matrix_lines);

// write final matrix of equivalence classes
void get_eq_classes(string& output, robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>>& bucket_class, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, ofstream* out);

void get_eq_classes_disk_query(string& output, robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>>& bucket_class, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, ofstream* out);

//sort count vectors by file, write one occurence per count in a new matrix file
void write_eq_class_matrix(string& output, vector<ofstream*>& all_files, uint64_t nb_unitigs,  bool do_query_on_disk, uint64_t nb_colors, ofstream* info);

#endif
