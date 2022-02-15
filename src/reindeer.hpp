
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
#include <set>
#include <mutex>
#include <ctype.h>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../blight/blight.h"
#include "../blight/zstr.hpp"
#include "../blight/strict_fstream.hpp"
#include "query.hpp"
#include "build_index.hpp"
#include "utils.hpp"

#ifndef REN
#define REN


using namespace std;
using namespace chrono;

template <class T>
class Reindeer_Index{
public:
// MPHF options
	uint m2;
	uint c;
	uint bit;
	uint ex;
	uint m1;
	uint m3;
	
	uint k; //size of k-mers
	
	// options
	uint threads; // number of threads
	bool record_counts; // if true, record counts else presence/absence //todo replace by color mode
	uint color_mode;
	bool do_query_on_disk; // if true the full index is dumped on the disk, else it is rebuilt and in ram
	
	// other counts modes
	bool quantize; // record quantization
	bool do_log; // record log
	
	// color variables
	uint64_t nb_colors; // number of samples
	uint64_t nb_monotig;
	long nb_eq_class;
	
	//files
	string reindeer_index_files;
	string color_dump_file;
	string monotig_files;
	string index_file;
	string matrix_eqc_info_file;
	string matrix_eqc_file;
	string matrix_eqc_position_file;
	bool dele_monotig_file;

	
	string color_load_file;
	string fof;
	string output;
	string dumped_index;
	string matrix_name;
	
	//constructor
	Reindeer_Index(uint pk, string& pfof, bool precord_counts, string& preindeer_index_files, uint pthreads,  bool pdo_query_on_disk, bool pquantize, bool pdo_log, uint pm1, uint pm3, bool dele_tmp);
	Reindeer_Index(string& output,string& output_query, uint threshold,  string& query, uint threads,  bool do_query_on_disk, bool dele_tmp);
	void print_Reindeer();


	void build_index(kmer_Set_Light* ksl);
	void read_info();
	void do_coloring(kmer_Set_Light* ksl,  vector<unsigned char*>& compr_monotig_color,vector<unsigned>& compr_monotig_color_size);
	kmer_Set_Light* load_rle_index(vector<unsigned char*> &compr_monotig_color,  vector<unsigned>& compr_monotig_color_sizes);
	void write_matrix_in_bucket_files(kmer_Set_Light* ksl, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size);
	void write_eq_class_matrix(vector<ofstream*>& all_files, ofstream* out_info);

	//~ Reindeer_Index()
	//~ {
		//~ high_resolution_clock::time_point t1 = high_resolution_clock::now();
		//~ kmer_Set_Light ksl(k,m1,m2,m3,threads,bit);
		//~ int systemRet;
		//~ // BUILD THE INDEX
		//~ nb_colors = get_color_number(fof);
		//~ build_index(k, m1, m2, m3, bit, color_load_file, color_dump_file, fof, record_counts, &ksl, threads,  output, do_query_on_disk, quantize, do_log, nb_colors);
		//~ high_resolution_clock::time_point t12 = high_resolution_clock::now();
		//~ duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
		//~ cout<<"Index building and Coloration done total: "<< time_span12.count() << " seconds."<<endl;
		//~ uint64_t mem(getMemorySelfMaxUsed());
		//~ cout<<"Max Memory used: "<< mem  << endl;
	//~ }


};

//~ void reindeer_index(uint k, string& fof,  string& color_dump_file, bool record_counts,  string& output, string& color_load_file, uint threads,  bool do_query_on_disk, bool quantize, bool do_log, uint m1, uint m3);
uint reindeer_query(string& output, string& output_query,   uint threshold,  string& query, uint threads, bool do_query_on_disk);
template class Reindeer_Index<uint8_t>;
template class Reindeer_Index<uint16_t>;
#endif
