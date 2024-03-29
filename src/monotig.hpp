#include "../blight/bbhash.h"
#include "../blight/blight.h"
#include "../blight/common.h"
#include "../blight/robin_hood.h"
#include "../blight/utils.h"
#include "../blight/zstr.hpp"
#include "../trle/trle.h"
#include "utils.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <tmmintrin.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#ifndef MINIT
#define MINIT

using namespace std;
using namespace chrono;

static ofstream monotigs_registry;

void set_monotigs_output(string& filename);

uint16_t abundance_at(uint8_t index);

void init_discretization_scheme();

uint8_t return_count_bin(uint16_t abundance);

uint16_t parseCoverage_bin(const string& str);

bool equal_nonull(const vector<uint16_t>& V1, const vector<uint16_t>& V2);

bool similar_count(const vector<uint16_t>& V1, const vector<uint16_t>& V2);

string compress_vector(const vector<uint16_t>& V);

void construct_index_fof(const string& input_file, vector<pair<string,uint64_t>>& kmers_by_file, const string& tmp_dir, int colormode);

vector<uint16_t> get_count_vector(const vector<pair<uint16_t, uint16_t>>& V, uint64_t number_color);
vector<uint8_t> get_color_vector(const vector<pair<uint16_t, uint16_t>>& V, uint64_t number_color);
void merge_super_buckets_mem(const string& input_file, uint64_t number_color, string& out_name, uint64_t number_pass, int colormode);

kmer select_good_successor(const robin_hood::unordered_node_map<kmer, kmer_context>& kmer2context, const kmer& start);

void write_buffer_count(vector<string>& buffers, ofstream* out, vector<uint16_t>& headerV, string& seq2dump, int32_t minimi);
void write_buffer_color(vector<string>& buffers, ofstream* out, vector<uint8_t>& headerV, string& seq2dump, int32_t minimi);

void get_monocolor_minitigs_mem(vector<robin_hood::unordered_node_map<kmer, kmer_context>>& min2kmer2context, ofstream* out, const vector<int32_t>& mini, uint64_t number_color, int colormode);

uint16_t parseCoverage(const string& str);

void create_super_buckets_list(const vector<string>& input_files, vector<pair<string,uint64_t>>& kmers_by_file);

#endif
