#ifndef QQ
#define QQ


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
#include <unordered_set>
#include <set>
#include <mutex>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include "matrix_operation.hpp"
#include "utils.hpp"
#include "../trle/trle.h"
using namespace std;
using namespace chrono;





unsigned char* decode_vector(unsigned char* minitig_counts, unsigned vector_size, uint64_t color_number)
{
	auto decoded_vector= new unsigned char[color_number*2];
	trled(minitig_counts, vector_size, decoded_vector, color_number*2);
	return decoded_vector;
}

vector<uint16_t> count_string_to_count_vector(unsigned char* count_char_minitig, unsigned size)
{
	vector<uint16_t> counts;
	for (uint i(0); i < size -1 ; i+=2)
	{
		counts.push_back((uint16_t)count_char_minitig[i] + (uint16_t)count_char_minitig[i+1]*256);
	}
	return counts;
}


vector<uint16_t> get_count_minitig(unsigned char* minitig_counts, unsigned vector_size, uint64_t color_number)
{
	unsigned char* decoded(decode_vector(minitig_counts, vector_size,color_number));
	vector<uint16_t> counts (count_string_to_count_vector(decoded, color_number*2));
	return counts;
}


void get_colors_counts(vector<int64_t>& kmer_ids, bool record_counts, uint64_t color_number, vector<int64_t>& kmers_colors, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts,vector<vector<uint32_t>>& color_me_amaze_reads, vector<vector<uint16_t>>& query_counts, vector<unsigned char*>& compr_minitig_color, vector<unsigned>&compr_minitig_color_size)
{
	vector<uint16_t> counts;
	for(uint64_t i(0);i<kmer_ids.size();++i)
	{
		if(kmer_ids[i]>=0)
		{
			vector<uint16_t> qcounts(get_count_minitig(compr_minitig_color[kmer_ids[i]], compr_minitig_color_size[kmer_ids[i]], color_number));
			if (not qcounts.empty()){
				query_counts.push_back( qcounts );
			}
		}
		
		
	}
	
}


void write_count_output(bool record_counts, vector<vector<uint16_t>>& query_counts, vector<string>& color_counts, uint64_t color_number )
{
		// compute a string that sums up the count(s) for each dataset
		vector<string> toW(color_number,"");
		for (auto&& c: query_counts)
		{
			for (uint color(0); color < c.size(); ++color)
			{
				toW[color]+=to_string(c[color]) + ":";
			}
		}
		for (uint str(0); str<toW.size(); ++str)
		{
			if (toW[str].empty())
			{
				color_counts.push_back("0");
			}
			else 
			{
				toW[str].pop_back();
				color_counts.push_back(toW[str]);
			}
		}
}


//~ void write_output(vector<int64_t>& kmers_colors, string& toWrite, bool record_reads, bool record_counts, vector<vector<uint32_t>>& query_unitigID,vector<vector<uint32_t>>& query_unitigID_tmp,  uint64_t& color_number, string& header, string& line, uint k, uint threshold, vector<string>& color_counts)
void write_output(vector<int64_t>& kmers_colors, string& toWrite, bool record_reads, bool record_counts, vector<vector<uint32_t>>& query_unitigID,vector<vector<uint32_t>>& query_unitigID_tmp,  uint64_t& color_number, string& header, string& line, uint k, uint threshold, vector<vector<uint16_t>>& query_counts)
{
	vector<string> color_counts;
	vector<string> toW(color_number,"");
	
	string nc;
	vector<string> last(color_number,"0");

	for (auto&& c: query_counts)
	{
		for (uint color(0); color < c.size(); ++color)
		{
			nc = to_string(c[color]);
			if (last[color] != nc )
			{
				toW[color]+=nc + ":";
				last[color] = nc;
			}
		}
		
	}
	for (uint str(0); str<toW.size(); ++str)
	{
		if (toW[str].empty())
		{
			color_counts.push_back("0");
		}
		else 
		{
			toW[str].pop_back();
			color_counts.push_back(toW[str]);
		}
	}
	vector<double_t> percent(color_number, 0);
	//~ for (auto&& c : kmers_colors)
	//~ {
		//~ percent[c]++;
	//~ }
	
	for (auto&& vec_c : query_counts)
	{
		for (uint c(0); c< vec_c.size(); ++c)
		{
			if (vec_c[c] > 0)
				percent[c]++;
		}
	}
	toWrite += header.substr(0,50) ;
	
	
	
	
	for (uint per(0); per < percent.size(); ++per)
	{
		percent[per] = percent[per] *100 /(line.size() -k +1);
		if (percent[per] >= (double_t) threshold )
		{
			if (not record_reads)
			{
				if (record_counts)
				{
					toWrite += "\t" + color_counts[per];
				} else {
					toWrite += "\t" + to_string((int)(percent[per]*10)/10) ;
				}
			}
			else
			{	//  appliquer aussi ici le threshold pour le cas on où query des reads
				query_unitigID.push_back(query_unitigID_tmp[percent[per]]);
			}
		} else {
			if (not record_reads)
				toWrite += "\t*";
		}
	}
	toWrite += "\n";
	//~ cout << "done" << endl;
}


void doQuery(string& input, string& name, kmer_Set_Light& ksl, uint64_t color_number, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts,vector<vector<uint32_t>>& color_me_amaze_reads, uint k, bool record_counts, bool record_reads, uint threshold,vector<vector<uint32_t>>& query_unitigID, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color, vector<unsigned >& compr_minitig_color_size){
	ifstream query_file(input);
	ofstream out(name);
	// #pragma omp parallel
	uint64_t num_seq(0);
	string qline;
	vector<string> lines;
	vector<vector<uint32_t>> query_unitigID_tmp;
	// FOR EACH LINE OF THE QUERY FILE
	while(not query_file.eof()){
		#pragma omp parallel num_threads(nb_threads)
		{
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
				if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T')
				{
					//~ cout << "new line" << endl;
					vector<int64_t> kmers_colors;
					vector<string> color_counts;
					vector<int64_t> kmer_ids;
					kmer_ids=ksl.get_rank_query(line);
					vector<vector<uint16_t>> query_counts;
					get_colors_counts(kmer_ids,  record_counts,  color_number, kmers_colors, color_me_amaze,  color_me_amaze_counts, color_me_amaze_reads, query_counts, compr_minitig_color, compr_minitig_color_size);
						write_output( kmers_colors, toWrite,  record_reads,  record_counts,  query_unitigID,query_unitigID_tmp,  color_number, header, line, k, threshold, query_counts);
				} else {
					if (line[0]=='@' or line[0]=='>')
						header = line;
				}
				#pragma omp ordered
				if (toWrite != header +"\n")
				{
					out<<toWrite;
				}
			}
		}
		lines={};
	}
	query_file.close();
	out.close();
}



void getReadsOfUnitig(string& bgreat_output_file, uint32_t unitigID, vector<uint64_t>& reads_u){
	string input(bgreat_output_file);
	vector<vector<uint64_t>> unitigs_to_nodes ;
	unitigs_to_nodes.reserve(10000);
	parse_bgreat_output(input, unitigs_to_nodes);
	dump_map("test_dump_vec", "test_dump_position", "test_dump_size", unitigs_to_nodes);
	load_unitig(unitigID,"test_dump_vec", "test_dump_position", "test_dump_size",reads_u);
	reads_u.clear();
	load_unitig(0,"test_dump_vec", "test_dump_position", "test_dump_size",reads_u);
}

//// todo: choisir le bon path file
//// todo: un threshold avant de renvoyer un id d'unitig dans les fonctions d'avant
void queryReadsID(vector<string>& bgreat_files, vector<vector<uint32_t>>& query_unitigID, string& outName)
{
	ofstream out(outName);
	vector<uint32_t> unitigID_vec;
	uint32_t real_ID;
	vector<uint64_t> reads_u;
	string toWriteTmp("");
	for (uint queryID(0); queryID < query_unitigID.size(); ++queryID)
	{
		unitigID_vec = query_unitigID[queryID];
		if (not unitigID_vec.empty())
		{
			for (uint graph(0); graph < unitigID_vec.size(); ++graph)
			{
				if (unitigID_vec[graph] > 0)
				{ // unitig IDs start at 1, a 0 means that no unitig matches in that graph
					real_ID = unitigID_vec[graph] -1;
					//todo ici on devrait passer un seul fichier pas le fof
					getReadsOfUnitig(bgreat_files[graph], real_ID, reads_u);
				}
				if (not reads_u.empty())
				{
					toWriteTmp +=  to_string(graph)+"-";
					sort(reads_u.begin(), reads_u.end());
					int prev_read(-1);
					for (auto&& read: reads_u)
					{
						if (((int)read) != prev_read) // remove duplicate ids
						{
							toWriteTmp += to_string(read) + ",";
						}
					} 
					toWriteTmp += " ";
				}
				
			}
		}
		if (not toWriteTmp.empty())
		{
			out << "dataset" << queryID <<":" <<toWriteTmp << endl;
		}
	}

}



string get_Bgreat_File(uint graph, string& bgread_output_file)
{
	string bg_file(""), line;
	ifstream fof(bgread_output_file);
	uint i(0);
	while(not fof.eof())
	{
		getline(fof,line);
		if (graph == i){break;}
		++i;
	}
	return bg_file;
}



void query_by_file(uint& counter, string& entry, kmer_Set_Light& ksl, uint64_t color_number,  vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, uint k, bool record_counts, bool record_reads, uint threshold, vector<string>& bgreat_files, string& output, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color, vector<unsigned >& compr_minitig_color_size)
{
	string outName(output + "/out_query_Reindeer" + to_string(counter) + ".out");
	vector<vector<uint32_t>> query_unitigID(color_number,{0});
	doQuery(entry, outName, ksl, color_number, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, k, record_counts, record_reads, threshold, query_unitigID, nb_threads, exact,  compr_minitig_color, compr_minitig_color_size);
	if (record_reads)
	{
		queryReadsID(bgreat_files, query_unitigID, outName);
	}
	counter++;
}


void perform_query(kmer_Set_Light& ksl, uint64_t color_number, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, uint k, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query, string& output, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color,  vector<unsigned>& compr_minitig_color_size)
{
	
	uint counter(0),patience(0);
	vector<string> bgreat_files;
	string entry;
	if (record_reads)
	{
		for (uint graph(0); graph < color_number; ++ graph)
		{
			bgreat_files.push_back(get_Bgreat_File(graph, bgreat_paths_fof));
		}

	}
	if (not query.empty())
	{
		if(exists_test(query))
		{
			high_resolution_clock::time_point t121 = high_resolution_clock::now();
			query_by_file(counter, query, ksl, color_number,  color_me_amaze,  color_me_amaze_counts, color_me_amaze_reads, k, record_counts, record_reads,  threshold, bgreat_files, output, nb_threads, exact, compr_minitig_color, compr_minitig_color_size);
			high_resolution_clock::time_point t13 = high_resolution_clock::now();
			duration<double> time_span13 = duration_cast<duration<double>>(t13 - t121);
			cout<<"Query done: "<< time_span13.count() << " seconds."<<endl;
		} else {
			cout<<"Empty query file !"<<endl;
		}
	} else {
		
		while(true){
			char str[256];
			cout << "Enter the name of the query file: ";

			cin.get (str,256);
			cin.clear();
			cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			entry = str;
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
					high_resolution_clock::time_point t121 = high_resolution_clock::now();
					query_by_file(counter, entry, ksl, color_number,  color_me_amaze,  color_me_amaze_counts, color_me_amaze_reads, k, record_counts, record_reads,  threshold, bgreat_files,output, nb_threads, exact,  compr_minitig_color, compr_minitig_color_size);
					memset(str, 0, 255);
					high_resolution_clock::time_point t13 = high_resolution_clock::now();
					duration<double> time_span13 = duration_cast<duration<double>>(t13 - t121);
					cout<<"Query done: "<< time_span13.count() << " seconds."<<endl;
				}else{
					cout << "The entry is not a file" << endl;
				}
			}
		}
	}
	delete(&ksl);
}

#endif
