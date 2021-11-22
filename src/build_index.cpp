
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




void get_header_monotig_file(zstr::ifstream& in, string& header)
{
	unsigned compressed_header_size;
	header = "";
	unsigned char* trash; //todo optim
	trash = new unsigned char[5];
	in.read((char*)trash, 5);// ">" and minimizer
	delete [] trash ;
	//~ cout << "position in file: " << in.tellg() << endl;
	//~ cout << "header capacity: " << header.capacity() << " header size before: " << header.size() << endl;
	in.read(reinterpret_cast<char *>(&compressed_header_size), sizeof(unsigned)); // size of colors/counts with rle
	header.append(reinterpret_cast<char *>(&compressed_header_size), sizeof(unsigned)); 
	//~ cout << "header capacity: " << header.capacity() << " header size before: " << header.size() << endl;
	if (header.capacity() <  sizeof(unsigned) + compressed_header_size)
		header.reserve(header.size() + sizeof(unsigned) + compressed_header_size +1024); //todo ok ?
	header.resize(sizeof(unsigned) + compressed_header_size );
	in.read((char*)(&header[sizeof(unsigned)]), compressed_header_size); //rle
	char c = in.get(); // \n
}



// dispatch count vectors in files. Similar counts go in similar files
void write_matrix_in_bucket_files(reindeer_index& index_values, kmer_Set_Light* ksl, vector <unsigned char*>& compressed_colors, vector <unsigned>& compressed_colors_size )
{
	//create bucket files for partitionning the compressed counts -> finding eq classes
	vector<ofstream*> all_files;
	for (uint i(0); i <1000; ++i)
	{
		ofstream* out = new ofstream(index_values.output + "/matrix_bucket_"+ to_string(i) + ".gz"); //todo changer
		all_files.push_back(out); 
	}
	string output_file_name;
	vector<string> monotig_files_names;
	string monotig_folder(index_values.output + "/monotig_files/");//todo changer
	get_all_blout(monotig_folder, monotig_files_names); 
	uint64_t nb_monotigs(0);
	ofstream out_info(index_values.output  + "/reindeer_matrix_eqc_info");//todo changer
	mutex mm;
	uint nb_treated_monotigs(0);
	uint i;
	#pragma omp parallel for num_threads(index_values.threads)
	for (i=0; i < monotig_files_names.size(); ++i)//loop on _blout files
	{
		string fname;
		fname = monotig_files_names[i];
		if (fname != "." and fname != "..")
		{
			fname = index_values.output + "/monotig_files/" + fname;
			zstr::ifstream monotigs_file(fname);
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
				if(monotig.empty() or header.empty()){continue;}
				monotig_id.clear();
				if (monotig[0] == 'A' or monotig[0] == 'C' or monotig[0] == 'G' or monotig[0] == 'T')
				{
					#pragma omp atomic
					++nb_monotigs;
					// get index from MPHF
					monotig_id=ksl->get_rank_query(monotig.substr(0,index_values.k)); // all kmers have the same id so we only query one
					if((not monotig_id.empty()) and monotig_id.back()>=0)
					{
						mm.lock();
						dump_compressed_vector_bucket(monotig_id.back(), all_files, header);
						nb_treated_monotigs++;
						mm.unlock();
					}
				}
				monotigs_file.peek();
			}
			if (DELE_MONOTIG_FILE)
				remove(fname.c_str()); 
		}
	}
	index_values.nb_monotig = nb_monotigs;
	for (uint i(0); i < all_files.size(); ++i)
	{
		all_files[i]->flush();
		all_files[i]->close();
	}
	out_info.write(reinterpret_cast<char*>(&index_values.nb_monotig),sizeof(uint64_t)); 
	out_info.write(reinterpret_cast<char*>(&index_values.k),sizeof(uint));  // in info: 1/nb monotigs, 2/k 3/record option 4/nb eq classes, 5/nb colors
	uint val(1);
	//~ if (! index_values.record_counts)//todo changer, calculé plusieurs fois
	//~ {
		//~ val = 0;
	//~ } 
	//~ else 
	//~ {
		//~ if (index_values.quantize)
		//~ {
			//~ val = 2;
		//~ }
		//~ else if (index_values.do_log)
		//~ {
			//~ val = 3;
		//~ }
		//~ else
		//~ {
			//~ val = 1;
		//~ }
	//~ }
	out_info.write(reinterpret_cast<char*>(&index_values.color_mode),sizeof(uint));
	// compute final equivalence class and write them
	write_eq_class_matrix(index_values.output, all_files, index_values.nb_monotig, index_values.do_query_on_disk, index_values.nb_colors, &out_info);//todo changer
	// remove bucket files
	for (uint i(0); i < all_files.size(); ++i)
	{
		string name(index_values.output + "/matrix_bucket_"+ to_string(i) + ".gz");
		remove(&name[0]); 
		delete all_files[i];
	}
	out_info.close();
}




// color using monotig file: either build and dump the color matrix during the index construction, or load it during the query
void do_coloring(reindeer_index& index_values, kmer_Set_Light* ksl,  vector<unsigned char*>& compr_monotig_color,vector<unsigned>& compr_monotig_color_size)
{
	vector <string> file_names;
	if (not index_values.color_load_file.empty()) //query
	{
		if (not index_values.do_query_on_disk)
		{
			compr_monotig_color = load_compressed_vectors(index_values, compr_monotig_color_size);//todo changer
		}
	} 
	else  //indexing
	{
		if(exists_test(index_values.fof))
		{
			ifstream fofin(index_values.fof);
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
		write_matrix_in_bucket_files(index_values,  ksl,  compr_monotig_color, compr_monotig_color_size); 
	}
}


// load dumped index(+colors)
kmer_Set_Light* load_rle_index(reindeer_index& index_values, vector<unsigned char*> &compr_monotig_color,  vector<unsigned>& compr_monotig_color_sizes)
{
	// build index from disk data
	kmer_Set_Light* ksl= new kmer_Set_Light(index_values.dumped_index);
	// compute colors
	do_coloring(index_values, ksl, compr_monotig_color, compr_monotig_color_sizes);
	//cleaning
	if (DELE_MONOTIG_FILE)
	{
		string cmd("rm -rf " + index_values.output +"/monotig_files"); 
		int sysRet(system(cmd.c_str()));
	}
	return ksl;
}


//todo can be moved before and int color mode become an index_values attribute
void get_color_mode(reindeer_index& index_values)
{
	if (index_values.record_counts)
		if (index_values.quantize)
			index_values.color_mode = 2;
		else
			if (index_values.do_log)
				index_values.color_mode = 3;
			else
				index_values.color_mode = 1;
	else 
		index_values.color_mode = 0;
}

// build index from new file
void build_index(reindeer_index& index_values, kmer_Set_Light* ksl)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
	bool dont_dump(false);
	
	// compute the number of distinct colors = number of samples
	get_color_mode(index_values);

	// no monotig file exist: build index from scratch
	if (not dirExists(index_values.monotig_files))
	{
			cout << "#Monotigs and index constuction..."<< endl;
			// index construction
			// apply monotig merge (-> MMM) with rule regarding colors (presence/absence) or counts: color 0, count 1, quantize 2, log 3
			cout << "mmmmmmmmmmmm" << (int)index_values.color_mode << endl;
			ksl->construct_index_fof(index_values.fof, index_values.output, (int)index_values.color_mode);
	} 
	else // monotig file exists
	{
		cerr << "[Warning] monotig files (monotig_files) were found in output dir, I will use them and I won't delete them" << endl;
		DELE_MONOTIG_FILE = false;
		// monotigs exist but the index was not built
		if (not exists_test(index_values.dumped_index))
		{
			// build file of file for index construction
			string fof_blout(do_fof(index_values.monotig_files, index_values.output)); 
			// index construction
			// apply monotig merge (-> MMM) with rule regarding colors (presence/absence) or counts: color 0, count 1, quantize 2, log 3
			ksl->construct_index_fof(fof_blout, index_values.output,(int)index_values.color_mode); 
			//todo could be optimized: blight needs to update some initialized variables and could take these directly as input super kmers
			remove((index_values.output + "/fof_blout").c_str());
		}
		else // monotig AND the index exist on disk
		{
			dont_dump = true;
			cerr << "[Warning] index file (reindeer_index.gz) was found in output dir, I will use it and I won't delete it" << endl;
			ksl = new kmer_Set_Light(index_values.dumped_index); //(BLight function)
		}
	}
	vector<unsigned char*>compr_monotig_color;
	vector<unsigned> compr_monotig_color_size;
	
	if (! dont_dump) // write index on disk
	{
		cout << "#Dumping index..."<< endl;
		
		ksl->dump_disk(index_values.dumped_index); //(BLight function to dump MPHF on disk)
		
		// log
		high_resolution_clock::time_point t2 = high_resolution_clock::now();
		duration<double> time_span12 = duration_cast<duration<double>>(t2 - t1);
		cout<<"Index written on disk: "<< time_span12.count() << " seconds."<<endl;
	}
	cout << "#Building colors and equivalence classes matrix to be written on disk..." << endl;
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	
	// compute colors for monotigs
	do_coloring(index_values, ksl,  compr_monotig_color, compr_monotig_color_size);
	
	// log and cleaning
	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	duration<double> time_span34 = duration_cast<duration<double>>(t4 - t3);
	cout<<"Matrix done: "<< time_span34.count() << " seconds."<<endl;
	if (DELE_MONOTIG_FILE)
	{
		string cmd("rm -r " + index_values.monotig_files); 
		int sysRet(system(cmd.c_str()));
	}
}

