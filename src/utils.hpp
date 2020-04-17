

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
#include "../blight/zstr.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <linux/limits.h>

#ifndef UTILS
#define UTILS


using namespace std;


void throw_character_issue() ;
void throw_size_issue() ;
bool check_character(string& s);

uint64_t xorshift ( uint64_t x );

string getLineFasta_buffer(ifstream* in);

vector<string> getLineFasta_buffer2(ifstream* in, uint stop, uint k);

bool is_empty_file(ifstream& file);

int dirExists(string& path);

inline bool exists_test(const string& name) {
	return ( access( name.c_str(), F_OK ) != -1 );
}

vector<string> split_utils(const string &s, char delim);

double parseCoverage_utils(const string& str);


uint32_t unitig_toui32(const string& s);


void parse_bgreat_output(string& input, vector<vector<uint64_t>>& unitigs_to_nodes);


vector<unsigned char> RLE16C(const vector<uint16_t>&V);


vector<uint16_t> RLE16D(const vector<uint8_t>&V);


vector<uint8_t> RLE8C(const vector<uint8_t>&V);

vector<uint8_t> RLE8D(const vector<uint8_t>&V);



void new_paired_end_file(string& input, string& input2, string& output_file, bool fastq);


void interleave_paired_end(string& fof, string& output);

uint64_t harmonic_mean(vector<uint64_t>& counts);

string getRealPath(string file, string& dir);

uint64_t get_color_number(string& fof);

#endif
