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
	string input,query,fof;
	uint k(0);
	uint m1(9);
	uint m2(8);
	uint m3(3);
	uint c(1);
	uint bit(6);
	uint ex(0);
	while ((ch = getopt (argc, argv, "g:q:k:m:n:s:t:b:e:o:")) != -1){
		switch(ch){
			case 'q':
				query=optarg;
				break;
			case 'g':
				input=optarg;
				break;
			case 'o':
				fof=optarg;
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
			case 'e':
				ex=stoi(optarg);
				break;
			case 'b':
				bit=stoi(optarg);
				break;
		}
	}

	if(query=="" or input=="" or fof=="" or k==0){
		cout
		<<"Mandatory arguments"<<endl
		<<"-g graph file constructed fom all your file"<<endl
		<<"-o your original files in a file of file"<<endl
		<<"-q query file"<<endl
		<<"-k k value used for graph "<<endl<<endl

		<<"Performances arguments"<<endl
		<<"-m minimizer size (9)"<<endl
		<<"-n to create 4^n mphf (7). More mean slower construction but better index, must be <=m"<<endl
		<<"-s to use 4^s files (3). More reduce memory usage and use more files, must be <=n"<<endl
		<<"-t core used (1)"<<endl
		<<"-b bit saved to encode positions (6). Will reduce the memory usage of b bit per kmer but query have to check 2^b kmers"<<endl;
		return 0;
	}
	{
		// I BUILD THE INDEX
		kmer_Set_Light ksl(k,m1,m2,m3,c,bit,ex);
		// IF YOU DONT KNOW WHAT TO DO THIS SHOULD WORKS GOOD -> kmer_Set_Light ksl(KMERSIZE,9,9,3,CORE_NUMBER,6,0);
		ksl.construct_index(input);

		// I PARSE THE FILE OF FILE
		ifstream fofin(fof);
		vector <string> file_names;
		string file_name;
		while(not fofin.eof()){
			getline(fofin,file_name);
			if(not file_name.empty()){
				file_names.push_back(file_name);
			}
		}
		uint64_t color_number(file_names.size());

		// I ALLOCATE THE COLOR VECTOR
		vector<bool> color_me_amaze((ksl.number_kmer)*color_number,false);
		//NOT VERY SMART I KNOW...

		// FOR EACH LINE OF EACH INDEXED FILE
		for(uint i_file(0);i_file<file_names.size();++i_file){
			ifstream in(file_names[i_file]);
			string line;
			vector<int64_t> kmer_ids;
			while(not in.eof()){
				getline(in,line);
				if(line[0]!= 'A' and line[0]!= 'C' and line[0]!= 'G' and line[0]!= 'T' ){
					continue;
				}
				//~ cout<<line<<endl;
				// I GOT THE IDENTIFIER OF EACH KMER
				kmer_ids=ksl.query_sequence_hash(line);
				//~ cout<<"HASH"<<endl;
				for(uint64_t i(0);i<kmer_ids.size();++i){
					//I COLOR THEM
					if(kmer_ids[i]>=0){
						//~ cout<<kmer_ids[i]*color_number+i_file<<"\n"<<endl;
						color_me_amaze[kmer_ids[i]*color_number+i_file]=true;
						//~ cout<<color_me_amaze.size()<<"cgood"<<endl;;
					}
				}
			}
		}


		cout<<"GO QUERY"<<endl;
		ifstream query_file(query);

		string line;
		vector<int64_t> kmer_ids;
		// FOR EACH LINE OF THE QUERY FILE
		while(not query_file.eof()){
			getline(query_file,line);
			if(line[0]!= 'A' and line[0]!= 'C' and line[0]!= 'G' and line[0]!= 'T' ){
				continue;
			}
			// I GOT THEIR INDICES
			kmer_ids=ksl.query_sequence_hash(line);
			for(uint64_t i(0);i<kmer_ids.size();++i){
				// KMERS WITH NEGATIVE INDICE ARE ALIEN/STRANGER/COLORBLIND KMERS
				if(kmer_ids[i]<0){
					cout<<"I WAS NOT INDEXED\n";
					continue;
				}
				// I KNOW THE COLORS OF THIS KMER !... I'M BLUE DABEDI DABEDA...
				cout<<"I HAVE BEEN SEEN IN THE FILE ";
				for(uint64_t i_color(0);i_color<color_number;++i_color){
					if(color_me_amaze[kmer_ids[i]*color_number+i_color]){
						cout<<i_color<<" ";
					}
				}
				cout<<"\n";
			}
		}
	}
	return 0;
}

