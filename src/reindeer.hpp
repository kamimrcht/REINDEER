
#include "../blight/blight.h"
#include "../blight/strict_fstream.hpp"
#include "../blight/zstr.hpp"
#include "build_index.hpp"
#include "utils.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <ctype.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <mutex>
#include <omp.h>
#include <pthread.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#ifndef REN
#define REN

using namespace std;
using namespace chrono;

template <class T>
class Reindeer_Index {
public:
    // MPHF options
    uint m2;
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
    string output_format; // output format of values (none/average/normalized/sum)
    bool output_monotigs;

    //variables
    vector<pair<string,uint64_t>> kmers_by_file; // total of kmers for each file
    uint threshold;
    vector<long> position_in_file;
    kmer_Set_Light* ksl;

    // other counts modes
    bool quantize; // record quantization
    bool do_log; // record log

    // color variables
    uint64_t nb_colors; // number of samples
    uint64_t nb_monotig {0};
    long nb_eq_class;

    vector<unsigned char*> compressed_monotig_color;
    vector<unsigned> compressed_monotig_color_sizes;

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
    Reindeer_Index(uint pk, string& pfof, bool precord_counts, string& preindeer_index_files, uint pthreads, bool pdo_query_on_disk, bool pquantize, bool pdo_log, uint pm1, uint pm3, bool dele_tmp, bool poutput_monotigs);
    Reindeer_Index(string& output, string& output_query, uint threshold, string& query, uint threads, bool dele_tmp, string& output_format);
    // methods
    void print_Reindeer();
    void load_index();
    void querying(string query);

    void build_index();
    void read_info();
    void do_coloring();
    kmer_Set_Light* load_rle_index();
    void write_matrix_in_bucket_files();
    void write_eq_class_matrix(vector<ofstream*>& all_files, ofstream* out_info);

    //query
    void perform_query(string& query);
    void query_by_file(uint& counter, string& entry);
    vector<long> get_position_vector_query_disk();
    void doQuery(string& input, string& name, vector<vector<uint32_t>>& query_unitigID);
    void get_colors_counts_query_eq_classes(vector<int64_t>& kmer_ids, vector<vector<uint16_t>>& query_counts);
    // for all queried k-mers, get the colors/counts in vector<vector<uint16_t>>& query_counts
    void get_colors_counts(vector<int64_t>& kmer_ids, vector<int64_t>& kmers_colors, vector<vector<uint16_t>>& query_counts);
    void write_output(string& toWrite, string& header, vector<vector<uint16_t>>& query_counts);
    // compute a string that sums up the count(s) for each dataset
    vector<uint> write_count_output(vector<vector<uint16_t>>& query_counts, vector<string>& toW, vector<string>& color_counts);
    void write_results_above_threshold(string& toWrite, vector<string>& color_counts, string& header, vector<uint>& covered_positions);
};

template class Reindeer_Index<uint8_t>;
template class Reindeer_Index<uint16_t>;
#endif
