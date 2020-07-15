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
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "../blight/blight.h"

using namespace std;
using namespace chrono;

int main(int argc, char** argv) {
	char ch;
	string input, inputfof, query;
	uint k(31);
	uint m1(10);
	uint m2(0);
	uint m3(4);
	uint c(1);
	uint bit(0);
	bool full(false), dump(false), hash(false);
	while ((ch = getopt(argc, argv, "hdag:q:k:m:n:s:t:b:e:f:")) != -1) {
		switch (ch) {
			case 'h': hash = true; break;
			case 'a': full = true; break;
			case 'd': dump = true; break;
			case 'q': query = optarg; break;
			case 'g': input = optarg; break;
			case 'f': inputfof = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'm': m1 = stoi(optarg); break;
			case 'n': m2 = stoi(optarg); break;
			case 's': m3 = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
			case 'b': bit = stoi(optarg); break;
		}
	}

	if (input == "" and query != "") {
		input = query;
	}
	if ((input == "" and inputfof == "") or k == 0) {
		cout << "Core arguments:" << endl
		     << "	-f file of  file" << endl
		     << "	-q query file" << endl
		     << "	-k k value used for graph (31) " << endl
		     << "By default only presence benchmark is performed \nUse -h for hash benchmark or -a for a complete test and verification" << endl
		     << "Performances arguments:" << endl
		     << "	-m minimizer size (10)" << endl
		     << "	-n to create 4^n mphf (m). More mean slower construction but better index, must be <=m" << endl
		     << "	-s to use 4^s files (4). More reduce memory usage and use more files, must be <=n" << endl
		     << "	-t core used (1)" << endl
		     << "	-b bit saved to encode positions (0). Will reduce the memory usage of b bit per kmer but query have to check 2^b kmers" << endl;
		return 0;
	}

	{
		cout << "I use -k " + to_string(k) + " -m  " + to_string(m1) + " -n  " + to_string(m2) + " -s  " + to_string(m3) + " -t " + to_string(c) + " -b "
		          + to_string(bit)
		     << endl;
		kmer_Set_Light ksl(k, m1, m2, m3, c, bit);
		if (inputfof != "") {
			cout << "Build index from list of file " << inputfof << " in folder wdir (if it does not exist please create it) " << endl;
			//~ ksl.construct_index_fof(inputfof);
			ksl.construct_index_fof(inputfof, "wdir", 1);
			// 0 for colors
			// 1 for exact counts
			// 2 for binned counts
			// 3 for log counts
		}

		if (not query.empty()) {
			cout << "NEW TESTS" << endl;
			if (hash) {
				ksl.file_query_hases(query, false);
			} else {
				ksl.file_query_all_test(query, full);
				vector<int64_t>monotig_id=ksl.get_rank_query("AGAAGGAGAAGGGAAAGGAAAGAG"); 
				if((not monotig_id.empty()) and monotig_id.back()>=0)
					cout << "OK" << endl;
				else
					cout << "BUG" << endl;
			}
		}

		if (dump) {
			cout << "DUMP" << endl;
			//~ ksl.dump_and_destroy("index.txt.gz");
			ksl.dump_disk("wdir/index.txt.gz");
			if (not query.empty()) {
				kmer_Set_Light ksl2("wdir/index.txt.gz");
				//~ for (uint i(0); i < 1; ++i)
					//~ ksl2.file_query_all_test(query, full);
				vector<int64_t>monotig_id=ksl2.get_rank_query("AGAAGGAGAAGGGAAAGGAAAGAG"); 
				if((not monotig_id.empty()) and monotig_id.back()>=0)
					cout << "OK" << endl;
				else
					cout << "BUG" << endl;
			}
		}

		cout << "I am glad you are here with me. Here at the end of all things." << endl;
	}
	return 0;
}
