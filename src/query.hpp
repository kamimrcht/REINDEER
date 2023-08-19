

#include "../blight/blight.h"
#include "../trle/trle.h"
#include "reindeer.hpp"
#include "utils.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
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
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef QQ
#define QQ

using namespace std;
using namespace chrono;

// decode rl encoded vector (color/count of ONE monotig)
unsigned char* decode_vector(unsigned char* monotig_counts, unsigned vector_size, uint64_t color_number, bool record_counts);

// for ONE monotig, get a vector of its counts/colors (uint) from the rle encoded matrix
vector<uint16_t> get_count_monotig(unsigned char* monotig_counts, unsigned vector_size, uint64_t color_number, bool record_counts);

long get_matrix_line_query_disk(int64_t rank, unsigned char* color, unsigned& line_size, vector<long>& position_in_file, ifstream& in);

long get_matrix_line_query(int64_t rank, unsigned char* color, unsigned& line_size, vector<long>& position_in_file, vector<unsigned char*>& compr_monotig_color);

#endif
