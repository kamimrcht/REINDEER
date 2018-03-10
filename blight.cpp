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

#include "ksl.h"



using namespace std;
using namespace chrono;




int main(int argc, char ** argv){
	if(argc<7){
		cout<<"[graph file] [kmer size] [minimizer size] [number superbuckets] [number core] [query file] "<<endl;
		exit(0);
	}
	string input(argv[1]);
	uint k(stoi(argv[2]));
	uint m(stoi(argv[3]));
	uint n(stoi(argv[4]));
	uint c(stoi(argv[5]));
	string query(argv[6]);
	kmer_Set_Light ksl(k,m,n,c);

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	ksl.abundance_minimizer_construct(input);
	ksl.create_super_buckets(input);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	//~ cout<<"Super bucket created: "<< time_span12.count() << " seconds."<<endl;

	ksl.read_super_buckets("_out");
	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);
	//~ cout<<"Bucket loaded: "<< time_span13.count() << " seconds."<<endl;

	//~ ksl.create_mphf();
	//~ high_resolution_clock::time_point t14 = high_resolution_clock::now();
	//~ duration<double> time_span14 = duration_cast<duration<double>>(t14 - t13);
	//~ cout<<"Mphf constructed: "<< time_span14.count() << " seconds."<<endl;

	//~ ksl.fill_positions();
	//~ high_resolution_clock::time_point t15 = high_resolution_clock::now();
	//~ duration<double> time_span15 = duration_cast<duration<double>>(t15- t14);
	//~ cout<<"Position  filled: "<< time_span15.count() << " seconds."<<endl;

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole indexing took me " << time_span.count() << " seconds."<< endl;

	ksl.multiple_query(query);
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
	cout << "The whole query took me " << time_span2.count() << " seconds."<<endl;
	cout<<"I am glad you are here with me. Here at the end of all things."<<endl;

	return 0;
}

