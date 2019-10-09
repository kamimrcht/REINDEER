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
//~ int main(int argc, char ** argv){
	//~ omp_set_nested(1);

	//~ char ch;
	//~ string input,query,fof, color_dump_file(""), color_load_file(""), bgreat_paths_fof("");
	//~ uint k(31);
	//~ uint m1(10);
	//~ uint m2(10);
	//~ uint m3(3);
	//~ uint c(1);
	//~ uint bit(6);
	//~ uint ex(0);
	//~ bool record_counts(false);
	//~ bool record_reads(false);
	//~ uint threshold(30);
	//~ while ((ch = getopt (argc, argv, "g:q:k:m:n:s:t:b:e:o:w:l:S:c:r:p:")) != -1){
		//~ switch(ch){
			//~ case 'q':
				//~ query=optarg;
				//~ break;
			//~ case 'w':
				//~ color_dump_file=optarg;
				//~ break;
			//~ case 'l':
				//~ color_load_file=optarg;
				//~ break;
			//~ case 'g':
				//~ input=optarg;
				//~ break;
			//~ case 'o':
				//~ fof=optarg;
				//~ break;
			//~ case 'k':
				//~ k=stoi(optarg);
				//~ break;
			//~ case 'm':
				//~ m1=stoi(optarg);
				//~ break;
			//~ case 'n':
				//~ m2=stoi(optarg);
				//~ break;
			//~ case 's':
				//~ m3=stoi(optarg);
				//~ break;
			//~ case 't':
				//~ c=stoi(optarg);
				//~ break;
			//~ case 'e':
				//~ ex=stoi(optarg);
				//~ break;
			//~ case 'b':
				//~ bit=stoi(optarg);
				//~ break;
			//~ case 'S':
				//~ threshold=stoi(optarg);
				//~ break;
			//~ case 'c':
				//~ record_counts=true;
				//~ break;
			//~ case 'r':
				//~ record_reads=true;
				//~ break;
			//~ case 'p':
				//~ bgreat_paths_fof=optarg;
				//~ break;
		//~ }
	//~ }
	//~ if(input=="" or (fof=="" and color_load_file=="") or k==0){
		//~ cout
		//~ <<"Mandatory arguments"<<endl
		//~ <<"-g graph file constructed fom all your file"<<endl
		//~ <<"-o your original files in a file of file OR -l a binary color matrix"<<endl
		//~ <<"-k k value used for graph "<<endl<<endl

		//~ <<"Performances arguments"<<endl
		//~ <<"-m minimizer size (9)"<<endl
		//~ <<"-n to create 4^n mphf (7). More mean slower construction but better index, must be <=m"<<endl
		//~ <<"-s to use 4^s files (3). More reduce memory usage and use more files, must be <=n"<<endl
		//~ <<"-t core used (1)"<<endl
		//~ <<"-b bit saved to encode positions (6). Will reduce the memory usage of b bit per kmer but query have to check 2^b kmers"<<endl << endl

		//~ <<"Options"<<endl
		//~ <<"-S minimum percent of query k-mers present in a dataset (default 30)" << endl
		//~ <<"-c record the k-mers counts"<<endl
		//~ <<"-r output the reads containing the query k-mers"<<endl
		//~ <<"-p bgreat paths file of file"<<endl
		//~ <<"-q query file (a fasta file). If no query file is provided, Reindeer will interactively ask for queries in the terminal."<< endl

		//~ <<"Serialization arguments"<<endl
		//~ <<"-w file where to write the index colors (-o is mandatory)"<<endl
		//~ <<"-l file from which index colors are loaded (do not use with -o)"<<endl;
		//~ return 0;
	//~ }
	//vector<mutex> MUTEXES(1000);
	//~ high_resolution_clock::time_point t1 = high_resolution_clock::now();
	//~ uint64_t color_number;
	//~ vector<vector<uint8_t>> color_me_amaze;
	//~ vector<vector<uint16_t>> color_me_amaze_counts;
	//~ vector<vector<uint32_t>> color_me_amaze_reads;
	//~ kmer_Set_Light ksl(k,m1,m2,m3,c,bit,ex);
	
	//~ // BUILD THE INDEX
	//~ build_index(k, m1, m2, m3, c, bit, ex, input, color_load_file, color_dump_file, fof, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, color_number, ksl);
	//~ high_resolution_clock::time_point t12 = high_resolution_clock::now();
	//~ duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	//~ cout<<"Index building and Coloration done: "<< time_span12.count() << " seconds."<<endl;
	//~ // QUERY //
	//~ perform_query(ksl, color_number, color_me_amaze,  color_me_amaze_counts,color_me_amaze_reads, k, record_counts,  record_reads,  threshold, bgreat_paths_fof, query);
	//~ return 0;
//~ }


void reindeer_index(uint k, string& graph, string& fof,  string& color_dump_file, bool record_counts, bool record_reads, string& output, string& color_load_file)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t color_number;
	vector<vector<uint8_t>> color_me_amaze;
	vector<vector<uint16_t>> color_me_amaze_counts;
	vector<vector<uint32_t>> color_me_amaze_reads;
	kmer_Set_Light ksl(k,m1,m2,m3,c,bit,ex);
	int systemRet;

	// BUILD THE INDEX
	if (not graph.empty())
	{
		high_resolution_clock::time_point t12 = high_resolution_clock::now();
		build_index(k, m1, m2, m3, c, bit, ex, graph, color_load_file, color_dump_file, fof, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, color_number, ksl);
		duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
		systemRet=system(("mv " + color_dump_file + " " + output ).c_str());
		cout<<"Index building and Coloration done: "<< time_span12.count() << " seconds."<<endl;
	}
}

//~ void reindeer_query(uint k, string& graph, string& fof,  string& color_dump_file, string& color_load_file, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query)
void reindeer_query(uint k, string& output,string& output_query, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query)
{
	// QUERY //
	string graph = getRealPath( "bcalm_union_out/union_graph.unitigs.fa", output);
	string fof( getRealPath("graphs.lst", output));
	//~ string graph("/home/camillemarchet/dev/test/REINDEER/output_reindeer/bcalm_union_out/union_graph.unitigs.fa");
	//~ string fof(output + "/home/camillemarchet/dev/test/REINDEER/output_reindeer/graphs.lst");
	string color_dump_file("");
	string color_load_file(getRealPath("reindeer_matrix", output));
	//~ string color_load_file("/home/camillemarchet/dev/test/REINDEER/output_reindeer/reindeer_matrix");
	uint64_t color_number;
	vector<vector<uint8_t>> color_me_amaze;
	vector<vector<uint16_t>> color_me_amaze_counts;
	vector<vector<uint32_t>> color_me_amaze_reads;
	kmer_Set_Light ksl(k,m1,m2,m3,c,bit,ex);
	build_index(k, m1, m2, m3, c, bit, ex, graph, color_load_file, color_dump_file, fof, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, color_number, ksl);

	perform_query(ksl, color_number, color_me_amaze,  color_me_amaze_counts,color_me_amaze_reads, k, record_counts,  record_reads,  threshold, bgreat_paths_fof, query, output_query);
}


#endif
