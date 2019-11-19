#ifndef REN
#define REN



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
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../blight/blight.h"
#include "../blight/zstr.hpp"
#include "../blight/strict_fstream.hpp"
#include "utils.hpp"
#include "matrix_operation.hpp"
#include "query.hpp"
#include "build_index.hpp"


using namespace std;
using namespace chrono;


// MPHF options
uint m1(10);
uint m2(10);
uint m3(3);
uint c(1);
uint bit(0);
uint ex(0);

void reindeer_index(uint k, string& fof,  string& color_dump_file, bool record_counts, bool record_reads, string& output, string& color_load_file, uint threads, bool exact)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t color_number;
	kmer_Set_Light ksl(k,m1,m2,m3,threads,bit);
	int systemRet;
	// BUILD THE INDEX
	build_index(k, m1, m2, m3, c, bit, color_load_file, color_dump_file, fof, record_counts, record_reads, color_number, ksl, threads, exact, output);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index building and Coloration done total: "<< time_span12.count() << " seconds."<<endl;
}




void reindeer_query(uint k, string& output,string& output_query, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query, uint threads, bool exact)
{
	// QUERY //
	string fof( getRealPath("graphs.lst", output));
	string color_dump_file("");
	string color_load_file(getRealPath("reindeer_matrix", output));
	uint64_t color_number(get_color_number(fof));
	uint nb_threads(1);
	vector<unsigned char*> compr_minitig_color;
	vector<unsigned> compr_minitig_color_sizes;
	cout << "\nLoading index.." << endl;
	kmer_Set_Light* ksl = load_rle_index(k, color_load_file, color_dump_file, fof, record_counts, record_reads, color_number, nb_threads, exact, output, compr_minitig_color, compr_minitig_color_sizes);
	cout << "\nComputing query..." << endl;
	perform_query(*ksl, color_number, k, record_counts,  record_reads,  threshold, bgreat_paths_fof, query, output_query, threads, exact, compr_minitig_color, compr_minitig_color_sizes);
	//~ delete [] ksl; //todo fix blight destructor
}


#endif
