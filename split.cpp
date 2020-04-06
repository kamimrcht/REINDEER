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
#include "blight.h"

using namespace std;
using namespace chrono;

int main(int argc, char** argv) {
	char ch;
	string input, query;
	uint k(31);
	uint m1(10);
	uint m2(0);
	uint m3(3);
	uint c(1);
	uint bit(0);
	while ((ch = getopt(argc, argv, "g:q:k:m:n:s:t:b:e:")) != -1) {
		switch (ch) {
			case 'g': input = optarg; break;
			case 'k': k = stoi(optarg); break;
			case 'm': m1 = stoi(optarg); break;
			case 's': m3 = stoi(optarg); break;
			case 't': c = stoi(optarg); break;
		}
	}

	if (input == "" or k == 0) {
		cout << "Mandatory arguments" << endl
		     << "-g fasta file" << endl
		     << "-k k value used for graph (31) " << endl
		     << endl
		     << "-m minimizer size (10)" << endl
		     << "-s to use 4^s files (3). More reduce memory usage and use more files, must be <=n" << endl
		     << "-t core used (1)" << endl;
		return 0;
	}
	{
		cout << "I use -g " + input + to_string(k) + " -m  " + to_string(m1) + " -s  " + to_string(m3) + " -t " + to_string(c) << endl;
		kmer_Set_Light ksl(k, m1, m1, m3, c, bit);
		ksl.create_super_buckets_regular(input, 1);
	}
	return 0;
}
