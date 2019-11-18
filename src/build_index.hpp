#ifndef BUI
#define BUI



#include "../trle/trle.h"

using namespace std;


vector<bool>  get_colors_minitigs(string& line)
{
	vector<bool> colors;
	vector<string> colors_minitig = split_utils(line,':');
	for (uint c(1); c < colors_minitig.size(); ++c) // convert the bit string to a bit vector
	{
		colors.push_back(stoi(colors_minitig[c]) > 0);
	}
	
	return colors;
}

vector<uint16_t> get_counts_minitigs(string& line)
{
	vector<uint16_t> counts;
	vector<string> colors_minitig = split_utils(line,':');
	for (uint c(1); c < colors_minitig.size(); ++c) // convert the bit string to a bit vector
	{
		counts.push_back((uint16_t) stoi(colors_minitig[c]));
	}
	return counts;
}

string get_counts_minitigs_char(string& line)
{
	//~ vector<uint16_t> counts;
	vector<string> colors_minitig = split_utils(line,':');
	string counts;
	for (uint c(1); c < colors_minitig.size(); ++c) // convert the bit string to a bit vector
	{
		//~ counts.push_back(colors_minitig[c]);
		//~ counts.push_back(':');
		counts+=(colors_minitig[c]);
		counts+=(':');
	}
	return counts;
}



void load_matrix(bool record_counts, bool record_reads, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, uint64_t& color_number, string& color_load_file)
{
	if (not record_counts)
		{
			if (not record_reads)
			{
				color_me_amaze = load_written_matrix(color_load_file);
				color_number = color_me_amaze.size();
			}
			else
			{
				color_me_amaze_reads = load_written_matrix_reads(color_load_file);
				color_number = color_me_amaze.size();
			}
	} else {
		color_me_amaze_counts = load_written_matrix_counts(color_load_file);
		color_number = color_me_amaze_counts.size();
	}
}

//~ void load_rle_matrix(bool record_counts, bool record_reads, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, uint64_t& color_number, string& color_load_file)
//~ {
	//~ vector<unsigned> vector_sizes;
	//~ uint64_t color_number;
	//~ uint64_t minitig_number;
	//~ load_compressed_vectors(color_load_file, vector_sizes, color_number,minitig_number);
//~ }


void build_matrix(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file)
{
	auto minitigs_file = new zstr::ifstream("_blmonocolor.fa.gz");
	string minitig;
	uint16_t count;
	uint32_t unitigID;
	vector<bool> colors;
	uint64_t nb_minitigs(0);

	//~ vector <unsigned char*> compressed_colors;
	//~ vector <unsigned> compressed_colors_size;
	vector<uint16_t> counts;
	//~ string counts;
	ofstream out(output_file);
	//+/ number of colors

	// wait before the real value
	out.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t)); // number of minitigs
	out.write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	unsigned line_size;
	string header;
	while(not minitigs_file->eof())
	{
		
		getline(*minitigs_file, minitig);
		if(minitig.empty()){continue;}
		string minitig_colorset(color_number, '0');
		string minitig_countlist;
		vector<int64_t> minitig_id;
		if (minitig[0] == 'A' or minitig[0] == 'C' or minitig[0] == 'G' or minitig[0] == 'T')
		{
			++nb_minitigs;
			minitig_id=ksl->get_rank_query(minitig); //normalement un seul id car on query tous les kmers d'un minitig

			uint64_t i(0);
					if(minitig_id[i]>=0)
					{
						// one vector corresponding to the minitig count/colors
						// convert to string for the compression
						uint n(counts.size()*2);
						uint nn(counts.size()*2 + 1024);
						unsigned char comp[nn];
						unsigned char *in = (unsigned char*)&counts[0];
						unsigned compr_vector_size = trlec(in, n, comp) ; //todo can I write it now on the disk
						out.write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned));
						int64_t tw(minitig_id[i]);
						out.write(reinterpret_cast<char*>(&tw),sizeof(int64_t));
						out.write((const char*)comp,(compr_vector_size));
					}
		} else { //header
			if (record_counts)
			{
				header = minitig;
				counts = get_counts_minitigs(minitig); // vector of ints
			} else {
				colors = get_colors_minitigs(minitig);
			}
		}
	}
	out.seekp(0, ios::beg);
	out.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t));
	delete(minitigs_file);
	out.close();
}


void dump_compressed_vectors(vector<unsigned>& minitigs_colors_size,vector<unsigned char*>& compr_minitigs_colors, const string& output_file, uint64_t color_number)
{
	ofstream out(output_file);
	out.write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	uint64_t minit_nb(compr_minitigs_colors.size());
	out.write(reinterpret_cast<char*>(&minit_nb),sizeof(uint64_t)); // number of minitigs
	unsigned line_size;
	for (uint i(0); i < compr_minitigs_colors.size(); ++i)
	{
		line_size = minitigs_colors_size[i];
		out.write(reinterpret_cast<char*>(&line_size),sizeof(unsigned));
		auto point =&(compr_minitigs_colors[i][0]);
		out.write((char*)point,(line_size));
	
	}
	
}



void dump_matrix(bool record_counts, bool record_reads, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, uint64_t& color_number, string& color_dump_file)
{
	if (not record_counts)
	{
		if (not record_reads)
		{
			write_color_matrix(color_dump_file, color_me_amaze);

		}
		else
		{
			write_color_matrix_reads(color_dump_file, color_me_amaze_reads);
		}
	} 
	else 
	{
		write_color_matrix_counts(color_dump_file, color_me_amaze_counts);
	}
}

void do_coloring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*>& compr_minitig_color,vector<unsigned>& compr_minitig_color_size)
{
	//~ vector <unsigned char*> compressed_colors;
	//~ vector <unsigned> compressed_colors_size;
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
	color_number = file_names.size();
	if (not record_counts)
	{
		color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl->number_super_kmer,0));
	}
	else 
	{
		color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl->number_super_kmer,0));
	}
	if (not color_load_file.empty())//todo: faire une fonction qui load tout l'index (blight se dump aussi)
	{
		//~ load_matrix( record_counts, record_reads,  color_me_amaze, color_me_amaze_counts,  color_me_amaze_reads,  color_number, color_load_file);
		//~ vector<unsigned> vector_sizes;
		uint64_t color_number;
		uint64_t minitig_number;
		compr_minitig_color = load_compressed_vectors(color_load_file, compr_minitig_color_size, color_number,minitig_number);
	} else {
		build_matrix( color_load_file,  color_dump_file,  fof,  ksl,  color_me_amaze,  color_me_amaze_counts, color_me_amaze_reads,  record_counts, record_reads,  k, color_number, nb_threads,  exact, output, compr_minitig_color, compr_minitig_color_size, color_dump_file);
	}
	//~ if (not color_dump_file.empty())
	//~ {	
		//dump_matrix( record_counts, record_reads,  color_me_amaze, color_me_amaze_counts,  color_me_amaze_reads,  color_number, color_dump_file);
		//~ dump_compressed_vectors(compr_minitig_color_size, compr_minitig_color, color_dump_file, color_number);
	//~ }
}



//~ kmer_Set_Light* load_index(uint k, string& color_load_file, string& color_dump_file, string& fof, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint64_t& color_number, uint nb_threads, bool exact, string& output)
//~ {
	//~ kmer_Set_Light* ksl= new kmer_Set_Light(output + "/reindeer_index.gz");
	//~ do_coloring(color_load_file, color_dump_file, fof, ksl, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, k, color_number, nb_threads, exact, output);
	//~ string cmd("rm -f _blmonocolor.fa.gz");
	//~ int sysRet(system(cmd.c_str()));
	//~ return ksl;
//~ }



kmer_Set_Light* load_rle_index(uint k, string& color_load_file, string& color_dump_file, string& fof, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*> &compr_minitig_color,  vector<unsigned>& compr_minitig_color_sizes)
{
	//~ kmer_Set_Light* ksl= new kmer_Set_Light(output + "/reindeer_index.gz");
	kmer_Set_Light* ksl= new kmer_Set_Light(output + "/reindeer_index.gz");
	do_coloring(color_load_file, color_dump_file, fof, ksl, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, k, color_number, nb_threads, exact, output, compr_minitig_color, compr_minitig_color_sizes);
	string cmd("rm -f _blmonocolor.fa.gz");
	int sysRet(system(cmd.c_str()));
	return ksl;
}




void build_index(uint k, uint m1,uint m2,uint m3, uint c, uint bit, string& color_load_file, string& color_dump_file, string& fof, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint64_t& color_number, kmer_Set_Light& ksl, uint nb_threads, bool exact, string& output)
{
	ksl.construct_index_fof(fof, true);
	vector<unsigned char*>compr_minitig_color;
	vector<unsigned> compr_minitig_color_size;
	do_coloring(color_load_file, color_dump_file, fof, &ksl, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, k, color_number, nb_threads, exact, output, compr_minitig_color, compr_minitig_color_size);
	ksl.dump_disk(output + "/reindeer_index.gz");
	string cmd("rm -f _blmonocolor.fa.gz");
	int sysRet(system(cmd.c_str()));
}


#endif
