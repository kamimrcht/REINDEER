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


//todo typedef for vector 16 bits (counts) / 8 bits (pres/abs)

inline bool exists_test(const string& name) {
    return ( access( name.c_str(), F_OK ) != -1 );
}


double parseCoverage(const string& str){
	size_t pos(str.find("km:f:"));
	if(pos==string::npos)
	{
		pos=(str.find("KM:f:"));
	}
	if(pos==string::npos){
		return 1;
	}
	uint i(1);
	while(str[i+pos+5]!=' '){
		++i;
	}
	// WARNING WE RETURN COUNTS BETWEEN 0 AND 255
	float result(stof(str.substr(pos+5,i)));
	if (result > 255){
		return 255;
	}
	else
	{
		return result;
	}
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
void write_color_matrix_counts(const string& output_file, vector<vector<uint16_t>>& color_matrix){
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
vector<vector<uint16_t>> load_written_matrix_counts(const string& input_file){
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
	vector<uint16_t> colorV(line_size,0);
	vector<vector<uint16_t>> color_matrix(color_number,colorV);
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



void doQuery(string input, string name, kmer_Set_Light& ksl, uint64_t color_number, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, uint k, bool record_counts, uint threshold){
	ifstream query_file(input);
	ofstream out(name);
	// #pragma omp parallel
	uint64_t num_seq(0);
	// mutex mm;
	// {
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
					vector<int64_t> color_counts(color_number,0);
					// I GOT THEIR INDICES
					vector<int64_t> kmer_ids=ksl.query_sequence_minitig(line);
					for(uint64_t i(0);i<kmer_ids.size();++i){
						// KMERS WITH NEGATIVE INDICE ARE ALIEN/STRANGER/COLORBLIND KMERS
						if(kmer_ids[i]>=0){
						// I KNOW THE COLORS OF THIS KMER !... I'M BLUE DABEDI DABEDA...
							if (not record_counts)
							{
								for(uint64_t i_color(0);i_color<color_number;++i_color){
									if(color_me_amaze[i_color][kmer_ids[i]]){
										kmers_colors.push_back(i_color);
										//~ if (record_counts)
										//~ {
											//~ color_counts[i_color] = color_me_amaze[i_color][kmer_ids[i]];
										//~ }
									}

								}
							} else {
								//~ cout << "here 1" << endl;
								//~ cout << 
								for(uint64_t i_color(0);i_color<color_number;++i_color)
								{
									if(color_me_amaze_counts[i_color][kmer_ids[i]])
									{
										kmers_colors.push_back(i_color);
										color_counts[i_color] = color_me_amaze_counts[i_color][kmer_ids[i]];
									}
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
							if (not record_counts)
							{
								percents[per].second = percents[per].second *100 /(line.size() -k +1);
								if (percents[per].second >= (double_t) threshold )
								{
							
									toWrite += " dataset" + to_string(percents[per].first+1) + ":" + to_string((int)(percents[per].second*10)/10) + "%";
								}
							}
							else
							{
								toWrite += " dataset" + to_string(percents[per].first+1) + ":" +  to_string(color_counts[percents[per].first]);
							}
						}
						toWrite += "\n";
					}
				}
				else if (line[0]=='@' or line[0]=='>')
				{
					header = line;
				}
				#pragma omp ordered
				if (toWrite != header +"\n")
				{
					out<<toWrite;
				}
			}
			lines={};
		}
	// }
	query_file.close();
	out.close();
}


void doColoring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light& ksl, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts,bool record_counts, uint k, uint64_t& color_number){
	if (color_load_file.empty()){ // use colors from the file of file
		//~ cout << "here" << endl;
		vector <string> file_names;
		if(exists_test(fof)){
			ifstream fofin(fof);
			string file_name;
			while(not fofin.eof()){
				getline(fofin,file_name);
				if(not file_name.empty()){
					if(exists_test(file_name)){
						file_names.push_back(file_name);
					}else{
						cout<<file_name<<" is not here"<<endl;
					}
				}
			}
		}else{
			cout<<"File of file problem"<<endl;
		}
		//~ cout << "here 111 " << endl;
		color_number = file_names.size();
		if (not record_counts)
		{
			color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.total_nb_minitigs,0));
		}
		else 
		{
			color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl.total_nb_minitigs,0));
			//~ cout << "here 3" << endl;
		}
		// FOR EACH LINE OF EACH INDEXED FILE
		uint i_file;
		#pragma omp parallel for
		for(i_file=0;i_file<file_names.size();++i_file){
			ifstream in(file_names[i_file]);
			string line;
			uint64_t count;

			// #pragma omp parallel num_threads(c)
			//{
				while(not in.eof()){
					// #pragma omp critical(i_file)
					// uint i_buffer;
					// #pragma omp parallel for
					// for(i_buffer=0;i_buffer<lines.size();++i_buffer){
						getline(in,line);
						if(line.empty()){continue;}
						if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T'){
							if (line.size() >= k){
							// I GOT THE IDENTIFIER OF EACH KMER
								auto minitig_ids=ksl.query_sequence_minitig(line);
								//~ cout << "here 31 " <<  minitig_ids.size() << " " << color_me_amaze.size() << endl;
								for(uint64_t i(0);i<minitig_ids.size();++i){
									//I COLOR THEM
									//~ cout << "check " << minitig_ids[i] << " " << color_me_amaze[0].size() <<  endl;
									
									if (not record_counts){
										if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze[0].size())
										{
											color_me_amaze[i_file][minitig_ids[i]]=1;
										}
									}else{
										//~ cout << "here 4" << endl;
										if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze_counts[0].size())
										{
											if (color_me_amaze_counts[i_file][minitig_ids[i]] == 0)
											{
												color_me_amaze_counts[i_file][minitig_ids[i]]=count;
											}
										}
											//~ if (line=="AAAAAAAAAACAAAAAATATAAAAAAAAAAAAAAAAAAA"){
												//~ cout << i_file << " " << minitig_ids[i] << " count " << count << endl;
											//~ }
										
									}
								}
								//~ if (line=="AAAAAAAAAACAAAAAATATAAAAAAAAAAAAAAAAAAA"){
												//~ cout << "end for" << endl;
								//~ }
							}
							//~ cout << "here 11" << endl;
						}
						else{
							if (record_counts){
								
								count = (uint16_t) parseCoverage(line);
								//~ if (line == ">1 LN:i:33 KC:i:6 km:f:2.0  L:-:808:+ L:-:814:+ L:-:160162:+"){
									//~ cout << count << endl;
									//~ cin.get();
								//~ }
							}
						}
					}
				//}
		}
		if (not color_dump_file.empty())
		{	
			if (not record_counts)
			{
				write_color_matrix(color_dump_file, color_me_amaze);
			} 
			else 
			{
				write_color_matrix_counts(color_dump_file, color_me_amaze_counts);
			}
		}
	}else{ // use color from file on disk
		if (not record_counts)
		{
			color_me_amaze = load_written_matrix(color_load_file);
			color_number = color_me_amaze.size();
		}
		else 
		{
			color_me_amaze_counts = load_written_matrix_counts(color_load_file);
			color_number = color_me_amaze_counts.size();
		}
	}
	
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
	bool record_counts(false);
	uint threshold(30);
	while ((ch = getopt (argc, argv, "g:q:k:m:n:s:t:b:e:o:w:l:S:c")) != -1){
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
			case 'S':
				threshold=stoi(optarg);
				break;
			case 'c':
				record_counts=true;
				break;
		}
	}
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

		<<"Options"<<endl
		<<"-c record the k-mers counts"<<endl
		<<"-S minimum percent of query k-mers present in a dataset (default 30)" << endl

		<<"Serialization arguments"<<endl
		<<"-w file where to write the index colors (-o is mandatory)"<<endl
		<<"-l file from which index colors are loaded (do not use with -o)"<<endl;
		return 0;
	}
	
	//~ vector<mutex> MUTEXES(1000);
	// I BUILD THE INDEX
	kmer_Set_Light ksl(k,m1,m2,m3,c,bit,ex);
	// IF YOU DONT KNOW WHAT TO DO THIS SHOULD WORKS GOOD -> kmer_Set_Light ksl(KMERSIZE,10,10,3,CORE_NUMBER,6,0);
	ksl.construct_index(input);
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	// I PARSE THE FILE OF FILE
	// I ALLOCATE THE COLOR VECTOR
	uint64_t color_number;
	vector<vector<uint8_t>> color_me_amaze;
	vector<vector<uint16_t>> color_me_amaze_counts;
	doColoring(color_load_file, color_dump_file, fof, ksl, color_me_amaze, color_me_amaze_counts, record_counts, k, color_number );
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
				high_resolution_clock::time_point t121 = high_resolution_clock::now();

				doQuery(entry, outName, ksl, color_number, color_me_amaze, color_me_amaze_counts, k, record_counts, threshold);
				memset(str, 0, 255);
				counter++;

				high_resolution_clock::time_point t13 = high_resolution_clock::now();
				duration<double> time_span13 = duration_cast<duration<double>>(t13 - t121);
				cout<<"Query done: "<< time_span13.count() << " seconds."<<endl;
			}else{
				cout << "The entry is not a file" << endl;
			}
		}
	}
	return 0;
}

