
void doColoring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light& ksl, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint k, uint64_t& color_number){
	if (color_load_file.empty()){ // use colors from the file of file
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
			color_me_amaze=vector<vector<uint8_t>>(color_number,vector<uint8_t>(ksl.total_nb_minitigs,0));
		}
		else 
		{
			color_me_amaze_counts=vector<vector<uint16_t>>(color_number,vector<uint16_t>(ksl.total_nb_minitigs,0));
		}
		// FOR EACH LINE OF EACH INDEXED FILE
		uint i_file;
		#pragma omp parallel for
		for(i_file=0;i_file<file_names.size();++i_file){
			ifstream in(file_names[i_file]);
			string line;
			uint16_t count;
			uint32_t unitigID;
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
							for(uint64_t i(0);i<minitig_ids.size();++i){
								if (not record_counts){
									if (not record_reads) // just remember presence/absence
									{
										if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze[0].size())
											{
												color_me_amaze[i_file][minitig_ids[i]]=1;
											}
									} else { // remember the unitig ID
										if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze[0].size())
											{
												color_me_amaze_reads[i_file][minitig_ids[i]]=unitigID; // WARNING: starts at 1
											}
									}
								}else{
									if(minitig_ids[i]>=0 and minitig_ids[i]<color_me_amaze_counts[0].size())
									{
										if (color_me_amaze_counts[i_file][minitig_ids[i]] == 0)
										{
											color_me_amaze_counts[i_file][minitig_ids[i]]=count;
										}
									}
								}
							}
						}
						unitigID ++;
					}
					else{
						if (record_counts){
							
							count = (uint16_t) parseCoverage(line);
						}
					}
				}
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
	}else{ // use color from file on disk
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
		}
		else 
		{
			color_me_amaze_counts = load_written_matrix_counts(color_load_file);
			color_number = color_me_amaze_counts.size();
		}
	}
	
}



void build_index(uint k, uint m1,uint m2,uint m3, uint c, uint bit, uint ex, string& input, string& color_load_file, string& color_dump_file, string& fof, vector<vector<uint8_t>>& color_me_amaze, vector<vector<uint16_t>>& color_me_amaze_counts, vector<vector<uint32_t>>& color_me_amaze_reads, bool record_counts, bool record_reads, uint64_t& color_number, kmer_Set_Light& ksl)
{
	ksl.construct_index(input);
	// PARSE THE FILE OF FILE
	// ALLOCATE THE COLOR VECTOR
	doColoring(color_load_file, color_dump_file, fof, ksl, color_me_amaze, color_me_amaze_counts, color_me_amaze_reads, record_counts, record_reads, k, color_number );
}
