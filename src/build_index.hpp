#ifndef BUI
#define BUI




using namespace std;

// read colors from bcalm headers
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

// read counts from bcalm headers
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




// build color/count matrix and dump it
void build_matrix(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file)
{
	auto minitigs_file = new zstr::ifstream("_blmonocolor.fa.gz");
	string minitig;
	uint16_t count;
	uint32_t unitigID;
	uint64_t nb_minitigs(0);
	vector<uint16_t> counts;
	ofstream out(output_file);
	// wait before the real value
	out.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t)); // number of minitigs
	out.write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	unsigned line_size;
	string header;
	unsigned char *in;
	while(not minitigs_file->eof())
	{
		
		getline(*minitigs_file, minitig);
		if(minitig.empty()){continue;}
		vector<int64_t> minitig_id;
		if (minitig[0] == 'A' or minitig[0] == 'C' or minitig[0] == 'G' or minitig[0] == 'T')
		{
			++nb_minitigs;
			minitig_id=ksl->get_rank_query(minitig); 
			uint64_t i(0);
			if(minitig_id[i]>=0) // even if minitig id is longer than 1, all ids should be the same
				dump_compressed_vector(counts, minitig_id[i], out, in);
		} else { //header
			header = minitig;
			counts = get_counts_minitigs(minitig); // vector of uints, colors and counts have the same encoding
			
		}
	}
	out.seekp(0, ios::beg);
	out.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t));
	delete(minitigs_file);
	out.close();
	memset(&in, 0, counts.size()*2);

}






// color using minitig file: either build and dump the color matrix during the index construction, or load it during the query
void do_coloring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*>& compr_minitig_color,vector<unsigned>& compr_minitig_color_size)
{
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
	if (not color_load_file.empty())
	{
		uint64_t color_number;
		uint64_t minitig_number;
		compr_minitig_color = load_compressed_vectors(color_load_file, compr_minitig_color_size, color_number,minitig_number);
	} else {
		build_matrix( color_load_file,  color_dump_file,  fof,  ksl,  record_counts, record_reads,  k, color_number, nb_threads,  exact, output, compr_minitig_color, compr_minitig_color_size, color_dump_file);
	}
}


// load dumped index(+colors)
kmer_Set_Light* load_rle_index(uint k, string& color_load_file, string& color_dump_file, string& fof, bool record_counts, bool record_reads, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*> &compr_minitig_color,  vector<unsigned>& compr_minitig_color_sizes)
{
	kmer_Set_Light* ksl= new kmer_Set_Light(output + "/reindeer_index.gz");
	do_coloring(color_load_file, color_dump_file, fof, ksl, record_counts, record_reads, k, color_number, nb_threads, exact, output, compr_minitig_color, compr_minitig_color_sizes);
	string cmd("rm -f _blmonocolor.fa.gz");
	int sysRet(system(cmd.c_str()));
	return ksl;
}



// build index from new file
void build_index(uint k, uint m1,uint m2,uint m3, uint c, uint bit, string& color_load_file, string& color_dump_file, string& fof, bool record_counts, bool record_reads, uint64_t& color_number, kmer_Set_Light& ksl, uint nb_threads, bool exact, string& output)
{
	cout << "Minitig coloring..."<< endl;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	// apply minitig merge (-> MMM) with rule regarding colors or counts
	if (record_counts)
		ksl.construct_index_fof(fof, true);
	else
		ksl.construct_index_fof(fof, false);
	vector<unsigned char*>compr_minitig_color;
	vector<unsigned> compr_minitig_color_size;
	do_coloring(color_load_file, color_dump_file, fof, &ksl, record_counts, record_reads, k, color_number, nb_threads, exact, output, compr_minitig_color, compr_minitig_color_size);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Coloration done: "<< time_span12.count() << " seconds."<<endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	cout << "Dump index..."<< endl;
	ksl.dump_disk(output + "/reindeer_index.gz");
	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	string cmd("rm -f _blmonocolor.fa.gz");
	int sysRet(system(cmd.c_str()));
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t2);
	cout<<"Index written on disk: "<< time_span13.count() << " seconds."<<endl;

}


#endif
