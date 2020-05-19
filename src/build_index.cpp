
#include "build_index.hpp"

using namespace chrono;
using namespace std;

bool DELE_MONOTIG_FILE(true);

// read colors from bcalm headers
vector<uint8_t>  get_colors_monotigs(string& line)
{
	vector<uint8_t> colors;
	vector<string> colors_monotig = split_utils(line,':');
	uint value;
	for (uint c(1); c < colors_monotig.size(); ++c) // convert the bit string to a bit vector
	{
		value = stoi(colors_monotig[c]) ;
		if (value > 0)
			colors.push_back(1);
		else
			colors.push_back(0);
	}
	return colors;
}

// read counts from bcalm headers
vector<uint16_t> get_counts_monotigs(string& line)
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




void get_header_monotig_file(ifstream& in, string& header)
{
	unsigned compressed_header_size;
	header = "";
	in.seekg(5, in.cur);// ">" and minimizer
	//~ cout << "position in file: " << in.tellg() << endl;
	//~ cout << "header capacity: " << header.capacity() << " header size before: " << header.size() << endl;
	in.read(reinterpret_cast<char *>(&compressed_header_size), sizeof(unsigned)); // size of colors/counts with rle
	header.append(reinterpret_cast<char *>(&compressed_header_size), sizeof(unsigned)); 
	//~ cout << "header capacity: " << header.capacity() << " header size before: " << header.size() << endl;
	if (header.capacity() <  sizeof(unsigned) + compressed_header_size)
		header.reserve(header.size() + sizeof(unsigned) + compressed_header_size +1024); //todo ok ?
	header.resize(sizeof(unsigned) + compressed_header_size );
	in.read((char*)(&header[sizeof(unsigned)]), compressed_header_size); //rle
	//~ cout << "header size after: " << header.size() << endl;
	//~ cout << "bytes header: " << (int)header[0] << " " << (int)header[1] <<" " <<(int)header[2]<< " " << (int)header[3]<<" " <<(int)header[4] << " " << (int)header[5] << endl; 
	//// debug
	//~ cout << "rd read comp header size: " << compressed_header_size << endl;
	//~ unsigned char *comp_written;
	//~ comp_written = new unsigned char[compressed_header_size];
	//~ cout << "size of unsigned: " << sizeof(unsigned) << endl;
	//~ comp_written = (unsigned char *)&header[sizeof(unsigned)];
	//~ unsigned char *decoded2;
	//~ decoded2 = new unsigned char [2*2 + 4096];
	//~ unsigned sz2 = trled(comp_written, compressed_header_size, decoded2, 2*2);
	//~ vector <uint16_t> decomp_count2 = count_string_to_count_vector(decoded2, sz2);
	//~ cout << "rd rle read (decoded): " ; 
	//~ for (uint i(0); i < decomp_count2.size() ; ++i)
		//~ cout << decomp_count2[i] << " " ;
	//~ cout << endl;
	//// /debug
	char c = in.get(); // \n
}



// dispatch count vectors in files. Similar counts go in similar files
void write_matrix_in_bucket_files(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, uint k, uint nb_threads,  string& output, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size, string& output_file, bool do_query_on_disk, uint64_t nb_colors, bool quantize, bool log )
{
	//create bucket files for partitionning the compressed counts -> finding eq classes
	vector<ofstream*> all_files;
	for (uint i(0); i <1000; ++i)
	{
		ofstream* out = new ofstream(output + "/matrix_bucket_"+ to_string(i)); //TODO zstr??
		all_files.push_back(out); 
	}
	string output_file_name;
	//~ string monotigs_fn(output +"/_blmonocolor.fa"); //monotigs
	//~ ifstream monotigs_file(monotigs_fn);
	vector<string> monotig_files_names;
	string monotig_folder(output + "/monotig_files/");
	get_all_blout(monotig_folder, monotig_files_names);
	uint64_t nb_monotigs(0);
	ofstream out_info(output_file+"_eqc_info");
	mutex mm;
	uint nb_treated_monotigs(0);
	//~ #pragma omp parallel num_threads(nb_threads)
	string fname;
	for (uint i(0); i < monotig_files_names.size(); ++i)//loop on _blout files
	{
		fname = monotig_files_names[i];
		if (fname != "." and fname != "..")
		{
			fname = output + "/monotig_files/" + fname;
			ifstream monotigs_file(fname);
			if(not exists_test(fname)){cerr << "File problem\n"; continue;}
			vector<int64_t> monotig_id;
			vector<uint16_t> counts;
			vector<uint8_t> colors;
			string header(5, ' '),monotig, buffer;
			header.reserve(4096);
			//~ unsigned char *in = nullptr;
			// record count vector for each monotig at index given by the mphf
			monotigs_file.peek();
			while(not monotigs_file.eof())
			{
				get_header_monotig_file(monotigs_file, header);
				getline(monotigs_file, monotig);
				//~ cout << "sequence monotig: " << monotig << endl;
				if(monotig.empty() or header.empty()){continue;}
				monotig_id.clear();
				if (monotig[0] == 'A' or monotig[0] == 'C' or monotig[0] == 'G' or monotig[0] == 'T')
				{
					mm.lock();
					++nb_monotigs;
					// get index form MPHF
					monotig_id=ksl->get_rank_query(monotig.substr(0,k)); // all kmers have the same id so we only query one
					mm.unlock();
					if((not monotig_id.empty()) and monotig_id.back()>=0)
					{
						mm.lock();
						//~ dump_compressed_vector_bucket(counts, monotig_id.back(), in, all_files,  colors, record_counts);
						dump_compressed_vector_bucket(monotig_id.back(), all_files, header);
						nb_treated_monotigs++;
						mm.unlock();
					}
				}
				monotigs_file.peek();
			}
			//~ remove(fname.c_str()); //TODO uncomment
		}
	}
	for (uint i(0); i < all_files.size(); ++i)
	{
		all_files[i]->flush();
		all_files[i]->close();
	}
	//~ out_nb.write(reinterpret_cast<char*>(&nb_monotigs),sizeof(uint64_t)); 
	out_info.write(reinterpret_cast<char*>(&nb_monotigs),sizeof(uint64_t)); 
	out_info.write(reinterpret_cast<char*>(&k),sizeof(uint));  // in info: 1/nb monotigs, 2/k 3/record option 4/nb eq classes, 5/nb colors
	uint val(1);
	if (! record_counts)
	{
		val = 0;
	} 
	else 
	{
		if (quantize)
		{
			val = 2;
		}
		else if (log)
		{
			val = 3;
		}
		else
		{
			val = 1;
		}
	}
	out_info.write(reinterpret_cast<char*>(&val),sizeof(uint));
	// compute final equivalence class and write them
	write_eq_class_matrix(output, all_files, nb_monotigs, do_query_on_disk, nb_colors, &out_info);
	// remove bucket files
	for (uint i(0); i < all_files.size(); ++i)
	{
		string name(output + "/matrix_bucket_"+ to_string(i));
		remove(&name[0]); 
		delete all_files[i];
	}
	out_info.close();
}




// color using monotig file: either build and dump the color matrix during the index construction, or load it during the query
void do_coloring(string& color_load_file, string& color_dump_file, string& fof, kmer_Set_Light* ksl, bool record_counts, uint k, uint nb_threads,  string& output, vector<unsigned char*>& compr_monotig_color,vector<unsigned>& compr_monotig_color_size, bool do_query_on_disk,  long& eq_class_nb, uint64_t& nb_colors, bool quantize, bool log)
{
	vector <string> file_names;
	if (not color_load_file.empty()) //query
	{
		if (not do_query_on_disk)
		{
			uint64_t color_number;
			uint64_t monotig_number;
			compr_monotig_color = load_compressed_vectors(color_load_file, compr_monotig_color_size, color_number, monotig_number, eq_class_nb);
		}
	} 
	else  //indexing
	{
		if(exists_test(fof))
		{
			ifstream fofin(fof);
			string file_name;
			while(not fofin.eof())
			{
				getline(fofin,file_name);
				if(not file_name.empty())
				{
					if(exists_test(file_name))
					{
						file_names.push_back(file_name);
					}
					else
					{
						cerr<< "[ERROR] " << file_name<<" is not here"<<endl;
					}
				}
			}
		}
		else
		{
			cerr<<"[ERROR] File of file problem"<<endl;
		}
		write_matrix_in_bucket_files(color_load_file,  color_dump_file,  fof,  ksl,  record_counts, k,  nb_threads, output, compr_monotig_color, compr_monotig_color_size, color_dump_file, do_query_on_disk, nb_colors,  quantize, log); 
	}
}


// load dumped index(+colors)
kmer_Set_Light* load_rle_index(uint k, string& color_load_file, string& color_dump_file, string& fof, bool record_counts,  uint nb_threads,  string& output, vector<unsigned char*> &compr_monotig_color,  vector<unsigned>& compr_monotig_color_sizes, bool do_query_on_disk, long& eq_class_nb, uint64_t& nb_colors,  bool quantize, bool log)
{
	kmer_Set_Light* ksl= new kmer_Set_Light(output + "/reindeer_index.gz");
	do_coloring(color_load_file, color_dump_file, fof, ksl, record_counts, k,  nb_threads, output, compr_monotig_color, compr_monotig_color_sizes, do_query_on_disk,  eq_class_nb, nb_colors,   quantize,  log);
	if (DELE_MONOTIG_FILE)
	{
		string cmd("rm -f " + output +"/_blmonocolor.fa"); //todo change
		int sysRet(system(cmd.c_str()));
	}
	return ksl;
}



// build index from new file
void build_index(uint k, uint m1,uint m2,uint m3, uint bit, string& color_load_file, string& color_dump_file, string& fof, bool record_counts,  kmer_Set_Light* ksl, uint nb_threads,  string& output, bool do_query_on_disk, bool quantize, bool do_log, uint64_t nb_colors)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	bool dont_dump(false);
	string in_name(output +"/monotig_files");
	int color_mode;
	if (not dirExists(in_name)) //todo
	{
			cout << "#Monotigs and index constuction..."<< endl;
			// apply monotig merge (-> MMM) with rule regarding colors or counts
			// color 0, count 1, quantize 2, log 3
			if (record_counts)
			{
				if (quantize)
				{
					color_mode = 2;
					ksl->construct_index_fof(fof, output, 2);
				}
				else
				{
					if (do_log)
					{
						color_mode = 3;
						ksl->construct_index_fof(fof, output, 3);
					}
					else
					{
						color_mode = 1; 
						ksl->construct_index_fof(fof, output, 1);
					}
				}
			}
			else
			{
				color_mode = 0;
				ksl->construct_index_fof(fof, output, 0);
			}
	} 
	else 
	{
		cerr << "[Warning] monotig files (monotig_files) were found in output dir, I will use them and I won't delete them" << endl;
		DELE_MONOTIG_FILE = false;
		if (not exists_test(output +"/reindeer_index.gz"))
		{
			// color 0, count 1, quantize 2, log 3
			ksl->construct_index(output +"/monotig_files",output);//todo 
		}
		else
		{
			dont_dump = true;
			cerr << "[Warning] index file (reindeer_index.gz) was found in output dir, I will use it and I won't delete it" << endl;
			ksl = new kmer_Set_Light(output + "/reindeer_index.gz");
		}
	}
	vector<unsigned char*>compr_monotig_color;
	vector<unsigned> compr_monotig_color_size;
	long eq_class_nb(0);
	
	if (! dont_dump)
	{
		cout << "#Dumping index..."<< endl;
		ksl->dump_disk(output + "/reindeer_index.gz");
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> time_span12 = duration_cast<duration<double>>(t2 - t1);
		cout<<"Index written on disk: "<< time_span12.count() << " seconds."<<endl;
	}
	cout << "#Building colors and equivalence classes matrix to be written on disk..." << endl;
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	do_coloring(color_load_file, color_dump_file, fof, ksl, record_counts,  k,  nb_threads,  output, compr_monotig_color, compr_monotig_color_size, do_query_on_disk, eq_class_nb, nb_colors,   quantize,  do_log);
	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	duration<double> time_span34 = duration_cast<duration<double>>(t4 - t3);
	cout<<"Matrix done: "<< time_span34.count() << " seconds."<<endl;
	if (DELE_MONOTIG_FILE)
	{
		string cmd("rm -r " + output +"/monotig_files"); 
		int sysRet(system(cmd.c_str()));
	}
}

