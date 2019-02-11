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
#include <omp.h>
#include "blight.h"



using namespace std;
using namespace chrono;



int main(int argc, char ** argv){
	omp_set_nested(1);

	char ch;
	string input,query,fof;
	uint k(0);
	uint m1(10);
	uint m2(10);
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
		// IF YOU DONT KNOW WHAT TO DO THIS SHOULD WORKS GOOD -> kmer_Set_Light ksl(KMERSIZE,10,10,3,CORE_NUMBER,6,0);
		ksl.construct_index(input);


		high_resolution_clock::time_point t1 = high_resolution_clock::now();


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

		// I ALLOCATE THE ABUNDANCE VECTOR
		vector<uint8_t> abundance(ksl.number_kmer,0);

		// FOR EACH LINE OF EACH INDEXED FILE
		uint i_file;
		#pragma omp parallel for num_threads(file_names.size())
		for(i_file=0;i_file<file_names.size();++i_file){
			ifstream in(file_names[i_file]);
			#pragma omp parallel num_threads(10)
			{
				vector<string> lines;
				string line;
				vector<int64_t> kmer_ids;
				while(not in.eof()){
					#pragma omp critical(i_file)
					{
						for(uint i(0);i<100;++i){
							getline(in,line);
							lines.push_back(line);
						}
					}
					for(uint i(0);i<100;++i){
						line=lines[i];
						if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T'){
							// I GOT THE IDENTIFIER OF EACH KMER
							kmer_ids=ksl.query_sequence_hash(line);
							for(uint64_t i(0);i<kmer_ids.size();++i){
								//I INCREMENT  THEIR ABUNDANCE
								if(kmer_ids[i]>=0){
									#pragma omp critical(color)
									{
										abundance[kmer_ids[i]]++;
									}
								}
							}
						}
					}
					lines={};
				}
			}
		}

		high_resolution_clock::time_point t12 = high_resolution_clock::now();
		duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
		cout<<"Abundancy computed: "<< time_span12.count() << " seconds."<<endl;

		ofstream out("outabundance");

		ifstream query_file(query);
		#pragma omp parallel
		{
			string qline;
			vector<string> lines;
			// FOR EACH LINE OF THE QUERY FILE
			while(not query_file.eof()){
				#pragma omp critical(i_file)
				{
					for(uint i(0);i<1000;++i){
						getline(query_file,qline);
						if(qline.empty()){break;}
						lines.push_back(qline);
					}
				}
				uint i;
				#pragma omp for ordered
				for(i=(0);i<lines.size();++i){
					string toWrite;
					string line=lines[i];
					if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T'){
						// I GOT THEIR INDICES
						vector<int64_t> kmer_ids=ksl.query_sequence_hash(line);
						for(uint64_t i(0);i<kmer_ids.size();++i){
							// KMERS WITH NEGATIVE INDICE ARE ALIEN/STRANGER/COLORBLIND KMERS
							if(kmer_ids[i]>=0){
							// I KNOW THE ABUNDANCY OF THIS KMER !
								toWrite+=to_string(abundance[kmer_ids[i]])+"	";
							}
							toWrite+="\n";
						}
					}
					#pragma omp ordered
					out<<toWrite;
				}
				lines={};
			}
		}

		high_resolution_clock::time_point t13 = high_resolution_clock::now();
		duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);
		cout<<"Query done: "<< time_span13.count() << " seconds."<<endl;
	}
	return 0;
}

