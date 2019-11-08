#ifndef BUI
#define BUI





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

void do_coloring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light& ksl, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output)
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
	if (not record_counts)
	{
		color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.number_super_kmer,0));
		//~ color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.total_nb_minitigs,0));
	}
	else 
	{
		//~ color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl.total_nb_minitigs,0));
		color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl.number_super_kmer,0));
	}
	if (not color_load_file.empty())//todo: faire une fonction qui load tout l'index (blight se dump aussi)
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
	} else {
		auto minitigs_file = new zstr::ifstream("_blmonocolor.fa.gz");
		string minitig;
		uint16_t count;
		uint32_t unitigID;
		vector<bool> colors;
		vector<uint16_t> counts;
		while(not minitigs_file->eof())
		{
			getline(*minitigs_file, minitig);
			if(minitig.empty()){continue;}
			vector<int64_t> minitig_id;
			if (minitig[0] == 'A' or minitig[0] == 'C' or minitig[0] == 'G' or minitig[0] == 'T')
			{
				minitig_id=ksl.get_rank_query(minitig); //normalement un seul id car on query tous les kmers d'un minitig
			//~ } else {
				//~ if (record_counts)
				//~ {
					//~ counts = get_counts_minitigs(minitig);
				//~ } else {
					//~ colors = get_colors_minitigs(minitig);
				//~ }
			//~ }
				for(uint64_t i(0);i<minitig_id.size();++i){ // should be a loop over 1 element
					if (not record_counts)
					{
						for (uint c(0); c < colors.size(); ++c)
						{
							if (colors[c]) // appears in this color
							{
									if (not record_reads) // just remember presence/absence
									{
										if(minitig_id[i]>=0 and minitig_id[i]<color_me_amaze[0].size())
										{
											color_me_amaze[c][minitig_id[i]]=1;
										}
									} else { // remember the unitig ID
										if(minitig_id[i]>=0 and minitig_id[i]<color_me_amaze[0].size())
										{
											color_me_amaze_reads[c][minitig_id[i]]=unitigID; // WARNING: starts at 1
										}
									}
							}
						}
					}else{
						for (uint c(0); c < counts.size(); ++c)
						{
							if(minitig_id[i]>=0 and minitig_id[i]<color_me_amaze_counts[0].size())
							{
								if (color_me_amaze_counts[c][minitig_id[i]] == 0)
								{
									color_me_amaze_counts[c][minitig_id[i]]=counts[c];
								}
							}
						}
					}
				}
			} else {
				if (record_counts)
				{
					counts = get_counts_minitigs(minitig);
				} else {
					colors = get_colors_minitigs(minitig);
				}
			}
		}
		delete(minitigs_file);

	}
	if (not color_dump_file.empty())
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
}



//~ void doColoring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light& ksl, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact){
	//~ if (color_load_file.empty()){ // use colors from the file of file
		//~ vector <string> file_names;
		//~ if(exists_test(fof)){
			//~ ifstream fofin(fof);
			//~ string file_name;
			//~ while(not fofin.eof()){
				//~ getline(fofin,file_name);
				//~ if(not file_name.empty()){
					//~ if(exists_test(file_name)){
						//~ file_names.push_back(file_name);
					//~ }else{
						//~ cout<<file_name<<" is not here"<<endl;
					//~ }
				//~ }
			//~ }
		//~ }else{
			//~ cout<<"File of file problem"<<endl;
		//~ }
		//~ color_number = file_names.size();
		//~ if (not record_counts)
		//~ {
			//~ // color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.total_nb_minitigs,0));
			//~ color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.number_super_kmer,0));
		//~ }
		//~ else 
		//~ {
			//~ color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl.number_super_kmer,0));
			//~ // color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl.total_nb_minitigs,0));
		//~ }
		//~ // FOR EACH LINE OF EACH INDEXED FILE
		//~ uint i_file;
		//~ #pragma omp parallel num_threads(nb_threads)
		//~ {
			//~ #pragma omp parallel for
			//~ for(i_file=0;i_file<file_names.size();++i_file){
				//~ ifstream in(file_names[i_file]);
				//~ string line;
				//~ uint16_t count;
				//~ uint32_t unitigID;
				//~ while(not in.eof()){
						//~ getline(in,line);
						//~ if(line.empty()){continue;}
						//~ if(line[0]=='A' or line[0]=='C' or line[0]=='G' or line[0]=='T'){
							//~ if (line.size() >= k){
								//~ vector<int64_t> minitig_ids;
							//~ // I GOT THE IDENTIFIER OF EACH KMER
								//~ minitig_ids=ksl.get_rank_query(line);
								//~ for(uint64_t i(0);i<minitig_ids.size();++i){
									//~ if (not record_counts){
										//~ if (not record_reads) // just remember presence/absence
										//~ {
											//~ if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze[0].size())
												//~ {
													//~ color_me_amaze[i_file][minitig_ids[i]]=1;
												//~ }
										//~ } else { // remember the unitig ID
											//~ if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze[0].size())
												//~ {
													//~ color_me_amaze_reads[i_file][minitig_ids[i]]=unitigID; // WARNING: starts at 1
												//~ }
										//~ }
									//~ }else{
										//~ if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze_counts[0].size())
										//~ {
											//~ if (color_me_amaze_counts[i_file][minitig_ids[i]] == 0)
											//~ {
												//~ color_me_amaze_counts[i_file][minitig_ids[i]]=count;
											//~ }
										//~ }
									//~ }
								//~ }
							//~ }
							//~ unitigID ++;
						//~ }
						//~ else{
							//~ if (record_counts){
								
								//~ count = (uint16_t) parseCoverage_utils(line);
							//~ }
						//~ }
					//~ }
			//~ }
		//~ }
		//~ if (not color_dump_file.empty())
		//~ cout << "dump" << endl;
		//~ cin.get();
		//~ {	
			//~ if (not record_counts)
			//~ {
				//~ if (not record_reads)
				//~ {
					//~ write_color_matrix(color_dump_file, color_me_amaze);

				//~ }
				//~ else
				//~ {
					//~ write_color_matrix_reads(color_dump_file, color_me_amaze_reads);

				//~ }
			//~ } 
			//~ else 
			//~ {
				//~ write_color_matrix_counts(color_dump_file, color_me_amaze_counts);
			//~ }
		//~ }
	//~ }else{ // use color from file on disk
		//~ if (not record_counts)
		//~ {
			//~ if (not record_reads)
			//~ {
				//~ color_me_amaze = load_written_matrix(color_load_file);
				//~ color_number = color_me_amaze.size();
			//~ }
			//~ else
			//~ {
				//~ color_me_amaze_reads = load_written_matrix_reads(color_load_file);
				//~ color_number = color_me_amaze.size();
			//~ }
		//~ }
		//~ else 
		//~ {
			//~ color_me_amaze_counts = load_written_matrix_counts(color_load_file);
			//~ color_number = color_me_amaze_counts.size();
		//~ }
	//~ }
	
//~ }



void build_index(uint k, uint m1,uint m2,uint m3, uint c, uint bit, string& color_load_file, string& color_dump_file, string& fof, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint64_t& color_number, kmer_Set_Light& ksl, uint nb_threads, bool exact, string& output)
{
	//~ ksl.construct_index(input);
	ksl.construct_index_fof(fof);
	// PARSE THE FILE OF FILE
	// ALLOCATE THE COLOR VECTOR
	//~ doColoring(color_load_file, color_dump_file, fof, ksl, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, k, color_number, nb_threads, exact);
	do_coloring(color_load_file, color_dump_file, fof, ksl, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, k, color_number, nb_threads, exact, output);
	string cmd("rm _blmonocolor.fa.gz");
	int sysRet(system(cmd.c_str()));

}


#endif
