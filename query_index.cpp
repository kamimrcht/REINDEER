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
	string input_index_file, query;
	while ((ch = getopt(argc, argv, "g:q:")) != -1) {
		switch (ch) {
			case 'q': query = optarg; break;
			case 'g': input_index_file = optarg; break;
		}
	}

	if (input_index_file == "" or query == "") {
		cout << "Core arguments:" << endl
		     << "	-g blight index file" << endl
		     << "	-q query string" << endl;
		return 0;
	}

	{
		kmer_Set_Light* ksl= new kmer_Set_Light(input_index_file);

		vector<bool> res = ksl->get_presence_query(query);
		cout << "Resulting presence vector :" << endl << "    ";
		for (const auto& b : res) {
		        cout << b << " ";
		}
		cout << endl;
	}
	return 0;
}
