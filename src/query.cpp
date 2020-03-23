#include "query.hpp"
using namespace std;
using namespace chrono;




// decode rl encoded vector (color/count of ONE minitig)
unsigned char* decode_vector(unsigned char* minitig_counts, unsigned vector_size, uint16_t color_number, bool record_counts)
{
	unsigned char* decoded_vector;
	unsigned sz;
	if (record_counts)
	{
		decoded_vector= new unsigned char[color_number*2 + 4096];
		 sz = trled(minitig_counts, vector_size, decoded_vector, color_number*2);
	}
	else
	{
		decoded_vector= new unsigned char[color_number + 1024];
		 sz = trled(minitig_counts, vector_size, decoded_vector, color_number);
	}
	return decoded_vector;
}

// convert char [] counts/colors to uint
vector<uint16_t> count_string_to_count_vector(unsigned char* count_char_minitig, unsigned size)
{
	vector<uint16_t> counts;
	for (uint i(0); i < size -1 ; i+=2)
	{
		counts.push_back((uint16_t)count_char_minitig[i] + (uint16_t)count_char_minitig[i+1]*256);
	}
return counts;
}

// convert char [] counts/colors to uint
vector<uint8_t> count_string_to_count_vector8(unsigned char* count_char_minitig, unsigned size)
{
	vector<uint8_t> counts;
	for (uint i(0); i < size  ; i++)
	{
		counts.push_back((uint8_t)count_char_minitig[i]);
	}
return counts;
}

// for ONE minitig, get a vector of its counts/colors (uint) from the rle encoded matrix
vector<uint16_t> get_count_minitig(unsigned char* minitig_counts, unsigned vector_size, uint16_t color_number, bool record_counts)
{
	unsigned char* decoded(decode_vector(minitig_counts, vector_size,color_number, record_counts));
	vector<uint16_t> counts (count_string_to_count_vector(decoded, color_number*2));
	delete [] decoded;
	return counts;
}



void get_position_vector_query_disk(vector<long>& position_in_file, string& position_file_name, uint64_t nb_minitig)
{
	zstr::ifstream * in = new zstr::ifstream(position_file_name); //positions file
	long position;
	uint nb(0);
	while (nb < nb_minitig)
	{
		in->read(reinterpret_cast<char*>(&position), sizeof(long));
		position_in_file.push_back(position);
		++nb;
	}
	delete in;
}

void get_matrix_line_query_disk(int64_t rank, unsigned char* color, unsigned& line_size, long position_in_matrix, ifstream& in )
{
	in.seekg(position_in_matrix);
	in.read(reinterpret_cast<char *>(&line_size), sizeof(unsigned));
	in.read(reinterpret_cast<char *>(&rank), sizeof(int64_t));
	in.read((char*)(color) , line_size);

}

long get_matrix_line_query(int64_t	rank, unsigned char* color, unsigned& line_size, vector<long>& position_in_file, vector<unsigned char*>& compr_minitig_color )
{
	
	long pos =  position_in_file[rank];
	color = compr_minitig_color[pos];
	return pos;
}

void get_colors_counts_query_eq_classes(vector<int64_t>& kmer_ids,   uint16_t color_number, vector<vector<uint16_t>>& query_counts,  vector<unsigned char*>& compr_minitig_color, vector<unsigned>&compr_minitig_color_size,vector<long>& position_in_file, bool record_counts, vector<vector<uint8_t>>& query_colors, bool do_query_on_disk, string& rd_file)
{

	ifstream in(rd_file);
	unsigned char* color;
	color = new unsigned char[2*color_number+1204];
	unsigned size;
	//uint count(0);
	int64_t lastId(-1);
	vector<uint16_t> qcounts, lastV;
	vector<uint8_t> qcolors;
	for(uint64_t i(0);i<kmer_ids.size();++i)
	{
		if(kmer_ids[i]>=0)
		{
			// no need to compute for a kmer if it comes from the same minitig than the previous, just copy the result
			if ( kmer_ids[i] == lastId)
			{
				lastId = kmer_ids[i];
				query_counts.push_back(lastV);
			}
			else
			{
				long pos ;
				unsigned char* lo;
				if (do_query_on_disk)
				{
					pos = position_in_file[kmer_ids[i]];
					int64_t rank;
					get_matrix_line_query_disk(rank, color, size, pos, in);
					vector<uint8_t> count_int = count_string_to_count_vector8(color, size);
					if (record_counts)
					{
						qcounts = RLE16D(count_int);
					}
					else
					{
						qcounts = {};
						vector<uint8_t> qcounts8 = RLE8D(count_int);
						copy(qcounts8.begin(), qcounts8.end(), back_inserter(qcounts));
					}
				}
				else
				{
					
					pos = get_matrix_line_query(kmer_ids[i], color, size, position_in_file, compr_minitig_color);
					lo = decode_vector((unsigned char*)&compr_minitig_color[pos][0], compr_minitig_color_size[pos], color_number, record_counts);
					if (record_counts)
					{
						qcounts = count_string_to_count_vector(lo, color_number*2);
					}
					else
					{
						qcounts = {};
						vector<uint8_t> qcounts8 = count_string_to_count_vector8(lo, compr_minitig_color_size[pos]);
						copy(qcounts8.begin(), qcounts8.end(), back_inserter(qcounts));
					}
					delete lo;	
				}
				lastV = qcounts;
				if (not qcounts.empty())
				{
					query_counts.push_back( qcounts );
				}
				
			}
		}
	}
	delete color;
	in.close();
}


// for all queried k-mers, get the colors/counts in vector<vector<uint16_t>>& query_counts
void get_colors_counts(vector<int64_t>& kmer_ids, bool record_counts, uint16_t color_number, vector<int64_t>& kmers_colors, vector<vector<uint16_t>>& query_counts, vector<unsigned char*>& compr_minitig_color, vector<unsigned>&compr_minitig_color_size, vector<long>& position_in_file)
{
	vector<uint16_t> counts;
	int64_t lastId(-1);
	vector<uint16_t> qcounts, lastV;
	for(int64_t i(0);i<kmer_ids.size();++i)
	{
		if(kmer_ids[i]>=0)
		{
			// no need to compute for a kmer if it comes from the same minitig than the previous, just copy the result
			if ( kmer_ids[i] == lastId)
			{
				lastId = kmer_ids[i];
				query_counts.push_back(lastV);
			}
			else
			{
				qcounts = get_count_minitig(compr_minitig_color[kmer_ids[i]], compr_minitig_color_size[kmer_ids[i]], color_number, record_counts);
				lastV = qcounts;
				if (not qcounts.empty())
					query_counts.push_back( qcounts );
		
			}
		}

	}
}


// compute a string that sums up the count(s) for each dataset
void write_count_output(bool record_counts, vector<vector<uint16_t>>& query_counts, uint16_t color_number , vector<string>& toW, vector<string>& color_counts)
{
	string nc("*");
	vector<string> last(color_number,"*");
	for (auto&& c: query_counts)
	{
		for (uint color(0); color < c.size(); ++color)
		{
			nc = to_string(c[color]);
			if (nc == "0")
					nc = "*";
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
			color_counts.push_back("*");
		}
		else 
		{
			toW[str].pop_back();
			while (toW[str].back() == '*' and (toW[str].size() > 1))
			{
				toW[str].pop_back();
				toW[str].pop_back();
			}
			color_counts.push_back(toW[str]);
		}
	}
}

void write_results_above_threshold(string& toWrite, vector<vector<uint16_t>>& query_counts, uint16_t color_number , vector<string>& toW, vector<string>& color_counts,  string& header, bool record_counts, uint threshold, string& line, uint k)
{
	vector<double_t> percent(color_number, 0);
	
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
			if (record_counts)
				toWrite += "\t" + color_counts[per];
			else 
				toWrite += "\t" + to_string((int)(percent[per]*10)/10) ;
		} else {
				toWrite += "\t*";
		}
	}
	toWrite += "\n";
}

void write_output(vector<int64_t>& kmers_colors, string& toWrite, bool record_reads, bool record_counts, vector<vector<uint32_t>>& query_unitigID,vector<vector<uint32_t>>& query_unitigID_tmp,  uint16_t& color_number, string& header, string& line, uint k, uint threshold,  vector<vector<uint16_t>>& query_counts, vector<vector<uint8_t>>& query_colors)
{
	vector<string> color_counts;
	vector<string> toW(color_number,"");
	write_count_output(record_counts, query_counts, color_number, toW, color_counts );
	write_results_above_threshold( toWrite,query_counts, color_number , toW,  color_counts,   header,  record_counts, threshold, line, k);
}






void doQuery(string& input, string& name, kmer_Set_Light& ksl, uint16_t& color_number, uint k, bool record_counts, bool record_reads, uint threshold,vector<vector<uint32_t>>& query_unitigID, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color, vector<unsigned >& compr_minitig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb, uint64_t nb_minitig){

	ifstream query_file(input);
	ofstream out(name);
	uint64_t num_seq(0);
	string qline;
	mutex mm;
	vector<string> lines;
	vector<vector<uint32_t>> query_unitigID_tmp;
	vector<long> position_in_file;
	// FOR EACH LINE OF THE QUERY FILE
	string position_file_name(rd_file+"_position.gz");
	get_position_vector_query_disk(position_in_file,  position_file_name,nb_minitig);
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
			
			#pragma omp for ordered
			//~ for(i=(0);i<lines.size();++i){
			for(i=(0);i<lines.size();i+=2){
				uint j(i);
				
				if (i%1000 == 0)
					cout << "-";
				string toWrite;
				string header;
				while (j < i+2)
				{
					string line=lines[j];
					if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T')
					{
						vector<int64_t> kmers_colors;
						vector<string> color_counts;
						vector<int64_t> kmer_ids;
						kmer_ids=ksl.get_rank_query(line);
						vector<vector<uint16_t>> query_counts;
						vector<vector<uint8_t>> query_colors;

						get_colors_counts_query_eq_classes( kmer_ids, color_number, query_counts,  compr_minitig_color, compr_minitig_color_size, position_in_file, record_counts, query_colors, do_query_on_disk, rd_file);
						mm.lock();
						write_output( kmers_colors, toWrite,  record_reads,  record_counts,  query_unitigID,query_unitigID_tmp,  color_number, header, line, k, threshold, query_counts, query_colors);
						mm.unlock();
					} else {
						if (line[0]=='>')
							header = line;
					}
					j++;
				}
				#pragma omp ordered
				mm.lock();
				if (toWrite != header +"\n")
				{
					out<<toWrite;
				}
				mm.unlock();
			}
		}
		lines={};
	}
	cout << endl;
	query_file.close();
	out.close();
}






void query_by_file(uint& counter, string& entry, kmer_Set_Light& ksl, uint16_t& color_number,   uint k, bool record_counts, bool record_reads, uint threshold, vector<string>& bgreat_files, string& output, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color, vector<unsigned >& compr_minitig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb, uint64_t	 nb_minitig)
{
	string outName(output + "/out_query_Reindeer" + to_string(counter) + ".out");
	vector<vector<uint32_t>> query_unitigID(color_number,{0});
	doQuery(entry, outName, ksl, color_number,  k, record_counts, record_reads, threshold, query_unitigID, nb_threads, exact,  compr_minitig_color, compr_minitig_color_size, do_query_on_disk, rd_file, eq_class_nb,  nb_minitig);
	counter++;
}


void perform_query(kmer_Set_Light& ksl, uint16_t& color_number,  uint k, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query, string& output, uint nb_threads, bool exact, vector<unsigned char*>& compr_minitig_color,  vector<unsigned>& compr_minitig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb)
{
	uint64_t nb_minitig;
	ifstream nb_minit_f(rd_file + "_minitig_nb");
	nb_minit_f.read(reinterpret_cast<char *>(&nb_minitig), sizeof(uint64_t));
	uint counter(0),patience(0);
	vector<string> bgreat_files;
	string entry;
	if (not query.empty())
	{
		if(exists_test(query))
		{
			high_resolution_clock::time_point t121 = high_resolution_clock::now();
			query_by_file(counter, query, ksl, color_number, k, record_counts, record_reads,  threshold, bgreat_files, output, nb_threads, exact, compr_minitig_color, compr_minitig_color_size, do_query_on_disk, rd_file, eq_class_nb, nb_minitig);
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
					query_by_file(counter, entry, ksl, color_number,   k, record_counts, record_reads,  threshold, bgreat_files,output, nb_threads, exact,  compr_minitig_color, compr_minitig_color_size, do_query_on_disk, rd_file, eq_class_nb, nb_minitig);
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
}

