#ifndef BUI
#define BUI

#include "eq_classes.hpp"


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
	uint pred(0);
	for (uint c(1); c < line.size();++c) // convert the bit string to a bit vector
	{
		if(line[c]==':'){
			if (pred!=0){
				counts.push_back(stoi(line.substr(pred,c-pred)));
			}
			pred=c+1;
		}

	}
	counts.push_back(stoi(line.substr(pred)));
	return counts;
}




// build color/count matrix and dump it
void build_matrix(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file)
{
	//~ auto minitigs_file = new zstr::ifstream(output +"/_blmonocolor.fa");
	ifstream minitigs_file(output +"/_blmonocolor.fa");
	uint64_t nb_minitigs(0);
	//~ ofstream out(output_file);
	auto out=new zstr::ofstream(output_file);
	ofstream out_nb(output_file+"_minitig_nb");
	mutex mm;
	//~ out.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t)); // number of minitigs - wait before the real value
	//~ out.write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	//~ out.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t)); // number of minitigs - wait before the real value
	out->write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	#pragma omp parallel num_threads(nb_threads)
	{
		vector<int64_t> minitig_id;
		vector<uint16_t> counts;
		string header,minitig, buffer;
		unsigned char *in;
		//~ while(not minitigs_file->eof())
		while(not minitigs_file.eof())
		{
			#pragma omp critical(infile)
			{
				//~ getline(*minitigs_file, header);
				//~ getline(*minitigs_file, minitig);
				getline(minitigs_file, header);
				getline(minitigs_file, minitig);

			}
			if(minitig.empty() or header.empty()){continue;}
			counts = get_counts_minitigs(header);
			minitig_id.clear();
			if (minitig[0] == 'A' or minitig[0] == 'C' or minitig[0] == 'G' or minitig[0] == 'T')
			{
				mm.lock();
				++nb_minitigs;
				mm.unlock();
				minitig_id=ksl->get_rank_query(minitig.substr(0,k)); // all kmers have the same id so we only query one
				if((not minitig_id.empty()) and minitig_id.back()>=0)
				{
					mm.lock();
					dump_compressed_vector_buff(counts, minitig_id.back(), buffer, in);
					mm.unlock();
					if(buffer.size()>40000)
					{
						#pragma omp critical(outfile)
						{
							//~ out.write(&buffer[0], buffer.size());
							out->write(&buffer[0], buffer.size());
						}
						buffer.clear();
					}
				}
			}
		}
		if(buffer.size()>0)
		{
			#pragma omp critical(outfile)
			{
				//~ out.write(&buffer[0], buffer.size());
				out->write(&buffer[0], buffer.size());
			}
			buffer.clear();
		}
		delete(in);
	}
	//~ out.seekp(0, ios::beg);
	out_nb.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t));
	//~ out.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t));
	//~ delete(minitigs_file);
	minitigs_file.close();
	//~ out.close();
	delete out;
	out_nb.close();
	string cmd("sort -u " + output + "/reindeer_matrix > "+ output + "/reindeer_matrix_sorted"  );
	int systemRet;
	systemRet = system((cmd).c_str());
	//~ cout << cmd << endl;
	//~ remove( (output + "/reindeer_matrix").c_str());
}


// sort the matrix lines to cluster similar counts/color together
// read this sorted matrix, group similar lines and rewrite only a single occurrence
//~ void compute_equivalence_class(minitig_nb)
//~ {
	//~ //sort matrix
	//~ string cmd("sort -t $'\xff' " + output + "/reindeer_matrix > "+ output + "/reindeer_matrix_sorted"  );
	//~ int systemRet;
	//~ systemRet = system((cmd).c_str());
	//~ //read matrix
	//~ ifstream sorted_mat_file(output + "/reindeer_matrix_sorted");
	//~ sorted_mat_file.read(reinterpret_cast<char *>(&color_number), sizeof(uint64_t));
	
	//~ vector<long>& position_in_original_file;
	//~ while (not sorted_mat_file.eof())
	//~ {
		//~ unsigned char* color;
		//~ color = new unsigned char[2*color_number+1204];
		//~ unsigned& line_size; 
		//~ in.read(reinterpret_cast<char *>(&line_size), sizeof(unsigned));
		//~ in.read(reinterpret_cast<char *>(&rank), sizeof(int64_t));
		//~ in.read((char*)(color) , line_size);
		
	//~ }
	//~ sorted_mat_file.close();

//~ }


void build_matrix_for_disk_query(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file)
{
	//~ auto minitigs_file = new zstr::ifstream(output +"/_blmonocolor.fa");
	ifstream  minitigs_file(output +"/_blmonocolor.fa");
	uint64_t nb_minitigs(0);
	//~ ofstream out(output_file);
	ofstream out(output_file);
	ofstream out_nb(output_file+"_minitig_nb");
	ofstream out_position(output_file+"_position");
	mutex mm;
	out.write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	uint nb_treated_minitigs(0);
	#pragma omp parallel num_threads(nb_threads)
	{
		vector<int64_t> minitig_id;
		vector<uint16_t> counts;
		string header,minitig, buffer;
		unsigned char *in;
		//~ while(not minitigs_file->eof())
		while(not minitigs_file.eof())
		{
			#pragma omp critical(infile)
			{
				//~ getline(*minitigs_file, header);
				//~ getline(*minitigs_file, minitig);
				getline(minitigs_file, header);
				getline(minitigs_file, minitig);

			}
			if(minitig.empty() or header.empty()){continue;}
			counts = get_counts_minitigs(header);
			minitig_id.clear();
			if (minitig[0] == 'A' or minitig[0] == 'C' or minitig[0] == 'G' or minitig[0] == 'T')
			{
				mm.lock();
				++nb_minitigs;
				minitig_id=ksl->get_rank_query(minitig.substr(0,k)); // all kmers have the same id so we only query one
				mm.unlock();
				if((not minitig_id.empty()) and minitig_id.back()>=0)
				{
					mm.lock();
					dump_compressed_vector(counts, minitig_id.back(), out, in, out_position);
					nb_treated_minitigs++;
					mm.unlock();
				}
			}
		}
		delete(in);
	}
	out_nb.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t));
	//~ delete(minitigs_file);
	minitigs_file.close();
	out.close();
	out_nb.close();
	out_position.close();
}



// dispatch count vectors in files. Similar counts go in similar files
void write_matrix_in_bucket_files(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file )
{

	vector<ofstream*> all_files;// todo move outside
	for (uint i(0); i < 60; ++i)
	{
		ofstream* out = new ofstream(output + "/matrix_bucket_"+ to_string(i)); //TODO zstr??
		all_files.push_back(out); 
	}
	//~ cout << "ici" << endl;
	string output_file_name;
	//~ auto minitigs_file = new zstr::ifstream(output +"/_blmonocolor.fa"); //TODO zstr
	string minitigs_fn(output +"/_blmonocolor.fa");
	//~ if (not exists_test(minitigs_fn))
		//~ cout << "pbl";
	//~ cout<<minitigs_fn<<endl;
	ifstream minitigs_file(minitigs_fn);
	//~ if(not minitigs_file.good()){
		//~ cout<<"fail open"<<endl;
	//~ }
	//~ cout << output +"/_blmonocolor.fa" << endl;
	//~ cin.get();
	uint64_t nb_minitigs(0);
	//~ ofstream out(output_file);
	ofstream out_nb(output_file+"_eqc_minitig_nb");
	ofstream out_position(output_file+"_position");
	mutex mm;
	//~ out.write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	uint nb_treated_minitigs(0);
	//~ cout << "ici3" << endl;
	//~ cout << output_file << endl;

	#pragma omp parallel num_threads(nb_threads)
	{
		vector<int64_t> minitig_id;
		vector<uint16_t> counts;
		string header,minitig, buffer;
		unsigned char *in;
		//~ while(not minitigs_file->eof())
		while(not minitigs_file.eof() and minitigs_file.good())
		{
			//~ cout << "here" << endl;
			//~ #pragma omp critical(infile)
			//~ {
				//~ getline(*minitigs_file, header);
				//~ getline(*minitigs_file, minitig);
				getline(minitigs_file, header);
				getline(minitigs_file, minitig);

			//~ }
			//~ cout << "here " << minitig << endl;
			//~ cout << "ici 2 " << endl;
			if(minitig.empty() or header.empty()){continue;}
			counts = get_counts_minitigs(header);
			minitig_id.clear();
			//~ cout << minitig[0] << endl;
			if (minitig[0] == 'A' or minitig[0] == 'C' or minitig[0] == 'G' or minitig[0] == 'T')
			{
				mm.lock();
				++nb_minitigs;
				minitig_id=ksl->get_rank_query(minitig.substr(0,k)); // all kmers have the same id so we only query one
				//~ cout << minitig_id << endl;
				mm.unlock();
				if((not minitig_id.empty()) and minitig_id.back()>=0)
				{
					mm.lock();
					//choose the output file + write compressed vector in it
					//~ cout << "HERE*********************" << endl;
					dump_compressed_vector_bucket(counts, minitig_id.back(), in, out_position, all_files);
					nb_treated_minitigs++;
					mm.unlock();
				}
			}
		}
		delete(in);
	}
	//~ cout << nb_minitigs << endl;
	out_nb.write(reinterpret_cast<char*>(&nb_minitigs),sizeof(uint64_t)); //TODO NB COLOR + /!\ relecture
	//~ delete(minitigs_file);
	//~ minitigs_file.close();
	//~ out.close();
	out_nb.close();
	out_position.close();
	for (uint i(0); i < all_files.size(); ++i)
	{
		all_files[i]->close();
	}
	write_eq_class_matrix(output, all_files, nb_minitigs, color_number);
	for (uint i(0); i < all_files.size(); ++i)
	{
		//~ all_files[i]->close();
		string name(output + "/matrix_bucket_"+ to_string(i));
		remove(&name[0]);
		delete all_files[i];
	}
}






// color using minitig file: either build and dump the color matrix during the index construction, or load it during the query
void do_coloring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, bool record_reads, uint k, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*>& compr_minitig_color,vector<unsigned>& compr_minitig_color_size, bool do_query_on_disk, string& nb_eq_class_file, long& eq_class_nb)
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
	if (not color_load_file.empty()) //query
	{
		ifstream nb_eq_f(nb_eq_class_file);
			//~ long nb_eq_class(0);
		nb_eq_f.read(reinterpret_cast<char *>(&eq_class_nb), sizeof(long));
		nb_eq_f.close();
		if (not do_query_on_disk)
		{
			uint64_t color_number;
			uint64_t minitig_number;
			compr_minitig_color = load_compressed_vectors(color_load_file, compr_minitig_color_size, color_number, minitig_number, eq_class_nb);
		}
	} 
	else  //indexing
	{
		if (do_query_on_disk)
			build_matrix_for_disk_query(color_load_file,  color_dump_file,  fof,  ksl,  record_counts, record_reads,  k, color_number, nb_threads,  exact, output, compr_minitig_color, compr_minitig_color_size, color_dump_file);
		else 
			//~ build_matrix(color_load_file,  color_dump_file,  fof,  ksl,  record_counts, record_reads,  k, color_number, nb_threads,  exact, output, compr_minitig_color, compr_minitig_color_size, color_dump_file);
			write_matrix_in_bucket_files(color_load_file,  color_dump_file,  fof,  ksl,  record_counts, record_reads,  k, color_number, nb_threads,  exact, output, compr_minitig_color, compr_minitig_color_size, color_dump_file);
	}
}


// load dumped index(+colors)
kmer_Set_Light* load_rle_index(uint k, string& color_load_file, string& color_dump_file, string& fof, bool record_counts, bool record_reads, uint64_t& color_number, uint nb_threads, bool exact, string& output, vector<unsigned char*> &compr_minitig_color,  vector<unsigned>& compr_minitig_color_sizes, bool do_query_on_disk, string& nb_eq_class_file, long& eq_class_nb)
{
	kmer_Set_Light* ksl= new kmer_Set_Light(output + "/reindeer_index.gz");
	do_coloring(color_load_file, color_dump_file, fof, ksl, record_counts, record_reads, k, color_number, nb_threads, exact, output, compr_minitig_color, compr_minitig_color_sizes, do_query_on_disk, nb_eq_class_file, eq_class_nb);
	string cmd("rm -f " + output +"/_blmonocolor.fa");
	int sysRet(system(cmd.c_str()));
	return ksl;
}



// build index from new file
void build_index(uint k, uint m1,uint m2,uint m3, uint c, uint bit, string& color_load_file, string& color_dump_file, string& fof, bool record_counts, bool record_reads, uint64_t& color_number, kmer_Set_Light& ksl, uint nb_threads, bool exact, string& output, bool do_query_on_disk, bool quantize, bool do_log)
{
	cout << "Minitig coloring..."<< endl;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	// apply minitig merge (-> MMM) with rule regarding colors or counts
	if (record_counts)
	{
		if (quantize)
		{
			ksl.construct_index_fof(fof, output, 2, 0);
		}
		else
		{
			if (do_log)
				ksl.construct_index_fof(fof, output, 3, 0);
			else
				ksl.construct_index_fof(fof, output, 1, 0);
		}
	}
	else
		ksl.construct_index_fof(fof, output, 0, 0);
	vector<unsigned char*>compr_minitig_color;
	vector<unsigned> compr_minitig_color_size;
	long eq_class_nb(0);
	string s("");
	do_coloring(color_load_file, color_dump_file, fof, &ksl, record_counts, record_reads, k, color_number, nb_threads, exact, output, compr_minitig_color, compr_minitig_color_size, do_query_on_disk, s, eq_class_nb);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Coloration done: "<< time_span12.count() << " seconds."<<endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	
	cout << "Dump index..."<< endl;
	ksl.dump_disk(output + "/reindeer_index.gz");
	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	string cmd("rm -f " + output +"/_blmonocolor.fa");
	int sysRet(system(cmd.c_str()));
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t2);
	cout<<"Index written on disk: "<< time_span13.count() << " seconds."<<endl;

}


#endif
