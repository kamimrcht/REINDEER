
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



void reindeer_index(uint k, string& fof,  string& color_dump_file, bool record_counts,  string& output, bool output_monotigs, string& color_load_file, uint threads,  bool do_query_on_disk, bool quantize, bool do_log, uint m1, uint m3);
uint reindeer_query(string& output, string& output_query,   uint threshold,  string& query, uint threads, bool do_query_on_disk);

#endif
