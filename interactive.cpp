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
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>
#include "blight.h"
#include "zstr.hpp"
#include "strict_fstream.hpp"


using namespace std;
using namespace chrono;



inline bool exists_test(const string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}


void write_color_matrix(const string& output_file, vector<vector<uint8_t>>& color_matrix){
	auto out=new zstr::ofstream(output_file);
	uint64_t color_number(color_matrix.size());
	uint64_t line_size(color_matrix[0].size());
	uint i(0);

	out->write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t));
	out->write(reinterpret_cast<char*>(&line_size),sizeof(uint64_t));
	//~ auto out=new zstr::ofstream(output_file);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		auto point =&(color_matrix[vec][0]);
		out->write((char*)point,(line_size));
	}
	delete out;
}



vector<vector<uint8_t>> load_written_matrix(const string& input_file){
	if(not exists_test(input_file))
	{
		cout<<"File problem"<<endl;
		exit(0);
	}
	uint64_t color_number;
	uint64_t line_size;
	auto in=new zstr::ifstream(input_file);
	in-> read(reinterpret_cast<char *>(&color_number), sizeof(uint64_t));
	in-> read(reinterpret_cast<char *>(&line_size), sizeof(uint64_t));
	cout << "line size " << line_size << " color number " << color_number << endl;
	vector<uint8_t> colorV(line_size,0);
	vector<vector<uint8_t>> color_matrix(color_number,colorV);
	//~ auto in=new zstr::ifstream(input);
	//~ assign(bloom_size/8,0);
	for (uint vec(0); vec < color_matrix.size(); ++vec)
	{
		in->read( (char*)((color_matrix[vec].data())) , line_size);
	}
	uint i(0);
	delete(in);
	return color_matrix;
}



void doQuery(string input, string name, kmer_Set_Light& ksl, uint64_t color_number, vector<vector<uint8_t>>& color_me_amaze, uint k){
	ifstream query_file(input);
	ofstream out(name);

	//~ #pragma omp parallel
	uint64_t num_seq(0);
	mutex mm;
	{
		string qline;
		vector<string> lines;
		// FOR EACH LINE OF THE QUERY FILE
		while(not query_file.eof()){
			#pragma omp critical(i_file)
			{
				for(uint i(0);i<4000;++i){
					getline(query_file,qline);
					if(qline.empty()){break;}
					lines.push_back(qline);
				}
			}
			uint i;
			string header;
			#pragma omp for ordered
			for(i=(0);i<lines.size();++i){
				string toWrite;
				string line=lines[i];
				if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T'){
					vector<int64_t> kmers_colors;
					// I GOT THEIR INDICES
					vector<int64_t> kmer_ids=ksl.query_sequence_minitig(line);
					for(uint64_t i(0);i<kmer_ids.size();++i){
						// KMERS WITH NEGATIVE INDICE ARE ALIEN/STRANGER/COLORBLIND KMERS
						if(kmer_ids[i]>=0){

						// I KNOW THE COLORS OF THIS KMER !... I'M BLUE DABEDI DABEDA...
							for(uint64_t i_color(0);i_color<color_number;++i_color){
								if(color_me_amaze[i_color][kmer_ids[i]]){
									kmers_colors.push_back(i_color);
								}

							}
						}
					}
					if (not  kmers_colors.empty())
					{
						sort(kmers_colors.begin(), kmers_colors.end());
						vector<pair<uint64_t, double_t>> percents;
						int64_t val(-1);
						for (uint64_t i_col(0); i_col < kmers_colors.size(); ++i_col)
						{
							if (kmers_colors[i_col] != val){
								percents.push_back({kmers_colors[i_col],1});
								val = kmers_colors[i_col];
							} else {
								percents.back().second++;
							}
						}
						toWrite += header ;
						for (uint per(0); per < percents.size(); ++per)
						{
							percents[per].second = percents[per].second *100 /(line.size() -k);
							if (percents[per].second >= 30 and percents[per].second < 101){
								toWrite += " dataset" + to_string(percents[per].first+1) + ":" + to_string(percents[per].second) + "%";
							}
						}
						toWrite += "\n";
					}
				}
				else if (line[0]=='@' or line[0]=='>')
				{
					//~ mm.lock();
					header = line;
					//~ mm.unlock();
				}
				#pragma omp ordered
				if (toWrite != header +"\n")
				{
					out<<toWrite;
				}
			}
			lines={};
		}
	}
	query_file.close();
	out.close();
}



int main(int argc, char ** argv){
	omp_set_nested(1);

	char ch;
	string input,query,fof, color_dump_file(""), color_load_file("");
	uint k(31);
	uint m1(10);
	uint m2(10);
	uint m3(3);
	uint c(1);
	uint bit(6);
	uint ex(0);
	while ((ch = getopt (argc, argv, "g:q:k:m:n:s:t:b:e:o:w:l:")) != -1){
		switch(ch){
			case 'q':
				query=optarg;
				break;
			case 'w':
				color_dump_file=optarg;
				break;
			case 'l':
				color_load_file=optarg;
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
	//~ cout << input << " " << fof << " " << color_load_file << k << endl;
	if(input=="" or (fof=="" and color_load_file=="") or k==0){
		cout
		<<"Mandatory arguments"<<endl
		<<"-g graph file constructed fom all your file"<<endl
		<<"-o your original files in a file of file OR -l a binary color matrix"<<endl
		<<"-k k value used for graph "<<endl<<endl

		<<"Performances arguments"<<endl
		<<"-m minimizer size (9)"<<endl
		<<"-n to create 4^n mphf (7). More mean slower construction but better index, must be <=m"<<endl
		<<"-s to use 4^s files (3). More reduce memory usage and use more files, must be <=n"<<endl
		<<"-t core used (1)"<<endl
		<<"-b bit saved to encode positions (6). Will reduce the memory usage of b bit per kmer but query have to check 2^b kmers"<<endl << endl

		<<"Serialization arguments"<<endl
		<<"-w file where to write the index colors (-o is mandatory)"<<endl
		<<"-l file from which index colors are loaded (do not use with -o)"<<endl;
		return 0;
	}
	{
		//~ vector<mutex> MUTEXES(1000);
		// I BUILD THE INDEX
		kmer_Set_Light ksl(k,m1,m2,m3,c,bit,ex);
		// IF YOU DONT KNOW WHAT TO DO THIS SHOULD WORKS GOOD -> kmer_Set_Light ksl(KMERSIZE,10,10,3,CORE_NUMBER,6,0);
		ksl.construct_index(input);



		high_resolution_clock::time_point t1 = high_resolution_clock::now();

		// I PARSE THE FILE OF FILE


		// I ALLOCATE THE COLOR VECTOR
		uint64_t color_number;
		vector<uint8_t> colorV(ksl.number_kmer,0);
		vector<vector<uint8_t>> color_me_amaze;
		//~ vector<vector<uint8_t>> color_me_amaze(color_number,colorV);
		if (color_load_file.empty()) // use colors from the file of file
		{
			ifstream fofin(fof);
			vector <string> file_names;
			string file_name;
			while(not fofin.eof()){
				getline(fofin,file_name);
				if(not file_name.empty()){
					file_names.push_back(file_name);
				}
			}

			//~ cout << "here" << endl;
			color_number = file_names.size();
			for (uint c(0); c < color_number; ++c){
				color_me_amaze.push_back(colorV);
			}
			//NOT VERY SMART I KNOW...

			// FOR EACH LINE OF EACH INDEXED FILE
			uint i_file;
			#pragma omp parallel for
			for(i_file=0;i_file<file_names.size();++i_file){
				ifstream in(file_names[i_file]);
				string read;
				//~ #pragma omp parallel num_threads(c)
				{
					vector<string> lines;
					while(not in.eof()){
						//~ #pragma omp critical(i_file)
						{
							for(uint i(0);i<1000;++i){
								getline(in,read);
								lines.push_back(read);
							}
						}
						uint i_buffer;
						//~ #pragma omp parallel for
						for(i_buffer=0;i_buffer<1000;++i_buffer){
							string line=lines[i_buffer];
							if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T'){
								// I GOT THE IDENTIFIER OF EACH KMER
								auto kmer_ids=ksl.query_sequence_hash(line);
								for(uint64_t i(0);i<kmer_ids.size();++i){
									//I COLOR THEM
									if(kmer_ids[i]>=0){
										//~ #pragma omp critical(color)
										//~ MUTEXES[(kmer_ids[i]*color_number+i_file)%1000].lock();
										{
											color_me_amaze[i_file][kmer_ids[i]]=1;
										}
										//~ MUTEXES[(kmer_ids[i]*color_number+i_file)%1000].unlock();
									}
								}
							}
						}
						lines={};
					}
				}
			}
			if (not color_dump_file.empty())
			{
				write_color_matrix(color_dump_file, color_me_amaze);
			}
			//~ cout << "here2" << endl;
		} else
		{ // use color from file on disk
			color_me_amaze = load_written_matrix(color_load_file);
			color_number = color_me_amaze.size();
		}
		//~ write_color_matrix("test_serial.txt", color_me_amaze);
		//~ color_me_amaze = load_written_matrix("test_serial.txt");
		high_resolution_clock::time_point t12 = high_resolution_clock::now();
		duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
		cout<<"Coloration done: "<< time_span12.count() << " seconds."<<endl;

		// query //
		uint counter(0),patience(0);
		while(true){
			char str[256];
			cout << "Enter the name of the query file: ";

			cin.get (str,256);
			cin.clear();
			cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			//~ cout << str << endl;
			string entry(str);
			if (entry == ""){
				if(patience>0){
					cout<<"See you soon !"<<endl;
					exit(0);
				}else{
					cout<<"Empty query file !  Type return again if you  want to quit"<<endl;
					patience++;
				}
			}else{
				patience=0;
				if(exists_test(entry)){
					string outName("out_query_BLight" + to_string(counter) + ".out");
					doQuery(entry, outName, ksl, color_number, color_me_amaze, k);
					memset(str, 0, 255);
					counter++;

					//~ high_resolution_clock::time_point t13 = high_resolution_clock::now();
					//~ duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);
					//~ cout<<"Query done: "<< time_span13.count() << " seconds."<<endl;
				}else{
					cout << "The entry is not a file" << endl;
				}
			}
		}
	}
	return 0;
}

