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




int main(int argc, char ** argv){
	char ch;
	string input,inputfof,query;
	uint k(31);
	uint m1(9);
	uint m2(0);
	uint m3(4);
	uint c(1);
	uint bit(0);
	bool full(false),dump(false);
	while ((ch = getopt (argc, argv, "dag:q:k:m:n:s:t:b:e:f:")) != -1){
		switch(ch){
			case 'a':
				full=true;
				break;
			case 'd':
				dump=true;
				break;
			case 'q':
				query=optarg;
				break;
			case 'g':
				input=optarg;
				break;
			case 'f':
				inputfof=optarg;
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


	if(query=="" and input!=""){
		query=input;
	}
	if(input=="" and query!=""){
		input=query;
	}
	cout<<inputfof<<endl;
	if((input=="" and inputfof=="") or k==0){
		cout
		<<"Mandatory arguments"<<endl
		<<"-g graph file"<<endl
		<<"-q query file"<<endl
		<<"-k k value used for graph (31) "<<endl<<endl
		<<"-m minimizer size (10)"<<endl
		<<"-n to create 4^n mphf (m). More mean slower construction but better index, must be <=m"<<endl
		<<"-s to use 4^s files (3). More reduce memory usage and use more files, must be <=n"<<endl
		<<"-t core used (1)"<<endl
		<<"-b bit saved to encode positions (6). Will reduce the memory usage of b bit per kmer but query have to check 2^b kmers"<<endl;
		return 0;
	}

	{
		cout<<"I use -k "+to_string(k)+" -m  "+to_string(m1)+" -n  "+to_string(m2)+" -s  "+to_string(m3)+" -t "+to_string(c)+" -b "+to_string(bit)<<endl;
		kmer_Set_Light ksl(k,m1,m2,m3,c,bit);
		if(inputfof!=""){
			cout<<"Build index from list of file "<<inputfof<<endl;
			ksl.construct_index_fof(inputfof);
		}else{
			cout<<"Build index from file "<<input<<endl;
			ksl.construct_index(input);
		}

		cout<<"NEW TESTS"<<endl;

		ksl.file_query_all_test(query,full);

		if(dump){
			cout<<"DUMP"<<endl;

			ksl.dump_disk("index.txt");
			kmer_Set_Light ksl2("index.txt");
			for(uint i(0);i<1;++i)
				ksl2.file_query_all_test(query,full);
		}

		cout<<"I am glad you are here with me. Here at the end of all things."<<endl;
	}
	return 0;
}

