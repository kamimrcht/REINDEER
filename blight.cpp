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
	char ch;
	string input,query;
	uint k(0);
	uint m1(7);
	uint m2(5);
	uint m3(3);
	uint c(1);
	uint bit(6);
	while ((ch = getopt (argc, argv, "g:q:k:m:n:s:t:b:")) != -1){
		switch(ch){
			case 'q':
				query=optarg;
				break;
			case 'g':
				input=optarg;
				break;
			case 'k':
				k=stoi(optarg);
				break;
			case 'm':
				m1=stoi(optarg);
				break;
			case 'n':
				m2=stoi(optarg);
				break;
			case 's':
				m3=stoi(optarg);
				break;
			case 't':
				c=stoi(optarg);
				break;
			case 'b':
				bit=stoi(optarg);
				break;
		}
	}

	if(query=="" or input=="" or k==0){
		cout
		<<"Mandatory arguments"<<endl
		<<"-g graph file"<<endl
		<<"-q query file"<<endl
		<<"-k k value used for graph "<<endl<<endl
		<<"-m minimizer size (7)"<<endl
		<<"-n to create 4^n mphf (5). More mean slower construction but better index, must be <=m"<<endl
		<<"-s to use 4^s files (3). More reduce memory usage and use more files, must be <=n"<<endl
		<<"-t core used (1)"<<endl
		<<"-b bit saved to encode positions (6). Will reduce the memory usage of b bit per kmer but query have to check 2^b kmers"<<endl;
		return 0;
	}
	kmer_Set_Light ksl(k,m1,m2,m3,c,bit);

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	ksl.abundance_minimizer_construct(input);
	cout<<"Abundance minimizer computed"<<endl;

	ksl.create_super_buckets(input);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Super bucket created: "<< time_span12.count() << " seconds."<<endl;

	ksl.read_super_buckets("_out");
	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole indexing took me " << time_span.count() << " seconds."<< endl;
	cin.get();
	ksl.file_query(query,false);
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
	cout << "The serial query took me " << time_span2.count() << " seconds."<<endl;

	ksl.file_query(query,true);
	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	duration<double> time_span3 = duration_cast<duration<double>>(t4 - t3);
	cout << "The optimized query took me " << time_span3.count() << " seconds."<<endl;

	cout<<"I am glad you are here with me. Here at the end of all things."<<endl;

	return 0;
}

