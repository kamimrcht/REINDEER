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
uint bit(6);
uint ex(0);

void reindeer_index(uint k, string& fof,  string& color_dump_file, bool record_counts, bool record_reads, string& output, string& color_load_file, uint threads, bool exact)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t color_number;
	vector<vector<uint8_t>> color_me_amaze;
	vector<vector<uint16_t>> color_me_amaze_counts;
	vector<vector<uint32_t>> color_me_amaze_reads;
	//~ kmer_Set_Light ksl(k,m1,m2,m3,c,bit,ex);
	kmer_Set_Light ksl(k,m1,m2,m3,c,bit);
	int systemRet;

	// BUILD THE INDEX
	//~ if (not graph.empty())
	//~ {
		high_resolution_clock::time_point t12 = high_resolution_clock::now();
		build_index(k, m1, m2, m3, c, bit, color_load_file, color_dump_file, fof, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, color_number, ksl, threads, exact, output);
		duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
		systemRet=system(("mv " + color_dump_file + " " + output ).c_str());
		cout<<"Index building and Coloration done: "<< time_span12.count() << " seconds."<<endl;
	//~ }
	//~ int systRet = system(cmd.c_str());
}


uint get_color_number(string& fof)
{	
	uint color(0);
	string line;
	ifstream fof_file(fof);
	while (not fof_file.eof())
	{
		getline(fof_file, line);
		if (not line.empty())
		{
			color++;
		}
	}
	return color;
}

void reindeer_query(uint k, string& output,string& output_query, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query, uint threads, bool exact)
{
	// QUERY //
	string graph = getRealPath( "bcalm_union_out/union_graph.unitigs.fa", output);
	string fof( getRealPath("graphs.lst", output));
	//~ string graph("/home/camillemarchet/dev/test/REINDEER/output_reindeer/bcalm_union_out/union_graph.unitigs.fa");
	//~ string fof(output + "/home/camillemarchet/dev/test/REINDEER/output_reindeer/graphs.lst");
	string color_dump_file("");
	string color_load_file(getRealPath("reindeer_matrix", output));
	//~ string color_load_file("/home/camillemarchet/dev/test/REINDEER/output_reindeer/reindeer_matrix");
	uint64_t color_number(get_color_number(fof));
	kmer_Set_Light ksl(k,m1,m2,m3,c,bit);
	vector<vector<uint8_t>> color_me_amaze;
	vector<vector<uint16_t>> color_me_amaze_counts;
	vector<vector<uint32_t>> color_me_amaze_reads;

	build_index(k, m1, m2, m3, c, bit, color_load_file, color_dump_file, fof, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, color_number, ksl, threads, exact, output);

	//~ color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.total_nb_minitigs,0));
	//~ color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.number_super_kmer,0));
	//~ color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl.total_nb_minitigs,0));
	//~ ksl.construct_index_fof(fof);
	cout << "\nComputing query..." << endl;
	perform_query(ksl, color_number, color_me_amaze,  color_me_amaze_counts,color_me_amaze_reads, k, record_counts,  record_reads,  threshold, bgreat_paths_fof, query, output_query, threads, exact);
}


#endif
