
#include "eq_classes.hpp"
#include "reindeer.hpp"

using namespace std;


void sort_vectors(vector<count_vector>& matrix_lines)
{
	sort(matrix_lines.begin(), matrix_lines.end(), compare_vec());
}


// write final matrix of equivalence classes
//~ void get_eq_classes(string& output, robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>>& bucket_class, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, zstr::ofstream* out)
void get_eq_classes(string& output, robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>>& bucket_class, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, ofstream* out)
{
	string prev_comp("");
	// for each representant, write the compressed size, index and compressed counts on disk
	// write a file:  integer i at position n shows which row i of the matrix should be looked up for monotig n (surjection monotigs -> class representants)
	for (auto rk(bucket_class.begin()); rk != bucket_class.end(); ++rk)
	{
		++prev_pos;
		count_vector v (rk->second.first);
		char* comp (&v.compressed[0]);
		out->write(reinterpret_cast<char*>(&v.compressed_size),sizeof(unsigned));
		out->write(reinterpret_cast<char*>(&v.monotig_rank),sizeof(int64_t));
		out->write((const char*)comp,(v.compressed_size));
		for (auto && rank: rk->second.second)
		{
			if (rank < final_positions.size())
				final_positions[rank] = prev_pos;
		}
	}
}

void get_eq_classes_disk_query(string& output, robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>>& bucket_class, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, ofstream* out)
{
	string prev_comp("");
	// for each representant, write the compressed size, index and compressed counts on disk
	// write a file:  integer i at position n shows which row i of the matrix should be looked up for monotig n (surjection monotigs -> class representants)
	for (auto rk(bucket_class.begin()); rk != bucket_class.end(); ++rk)
	{
		++prev_pos;
		count_vector v (rk->second.first);
		char* comp (&v.compressed[0]); //homemade RLE compression
		unsigned compr_vector_size = v.compressed.size();
		long position(out->tellp());
		out->write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned));
		int64_t tw(v.monotig_rank);
		out->write(reinterpret_cast<char*>(&tw),sizeof(int64_t));
		out->write((const char*)&comp[0],(compr_vector_size));
		char* returnc = new char[1];
		returnc[0] = (uint8_t) 255;
		out->write(&returnc[0], sizeof(char)); // line separator	
		delete [] returnc;
		for (auto && rank: rk->second.second)
		{
			if (rank < final_positions.size())
				final_positions[rank] = position;//here position is the position in file given by tellp
		}
	}
}

//sort count vectors by file, write one occurence per count in a new matrix file
//~ void write_eq_class_matrix(string& output, vector<ofstream*>& all_files, uint64_t nb_unitigs, bool do_query_on_disk, uint64_t nb_colors, ofstream* out_info)
template <class T>
void Reindeer_Index<T>::write_eq_class_matrix(vector<ofstream*>& all_files, ofstream* out_info)
{
	cout << "Sorting datasets to find equivalence classes..." << endl;
	vector<long> final_positions(nb_monotig);
	//~ zstr::ofstream * outz;
	ofstream *out = nullptr;
	//~ if (do_query_on_disk)
	//~ {
		//~ out = new ofstream(output + "/reindeer_matrix_eqc");
			
		out = new ofstream(matrix_eqc_file);
		out->write(reinterpret_cast<char*>(&nb_colors),sizeof(uint64_t)); // number of colors
	//~ }
	//~ else
	//~ {
		//~ outz = new zstr::ofstream(output + "/reindeer_matrix_eqc.gz");
		//~ outz->write(reinterpret_cast<char*>(&nb_colors),sizeof(uint64_t)); // number of colors

	//~ }

	long prev_pos(-1);
	uint i(0);
	mutex mm; 
	//~ #pragma omp parallel for num_threads(4)
	// in each bucket file
	for (i=0; i < all_files.size(); ++i)
	{
		int64_t rank;
		uint64_t sizecp(nb_colors*2 + 1024);
		char comp[sizecp];
		unsigned comp_size;
		vector<count_vector> count_vecs;
		robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>> bucket_class;
		string compressed;
		//~ zstr::ifstream in(output + "/matrix_bucket_"+ to_string(i) + ".gz"); //TODO zstr??
		zstr::ifstream in(reindeer_index_files + "/matrix_bucket_"+ to_string(i) + ".gz"); //TODO zstr??
		if (not is_empty_zfile(in)){
			in.peek();
			while (not in.eof() and not in.fail())
			{
				//read file and store each vector in struct
				read_matrix_compressed_line(in, rank, comp, comp_size);
				compressed.assign(comp, comp_size);
				count_vector v({comp_size, rank, compressed});
				// insert a single class representant in hash map
				if (not bucket_class.count(v.compressed))
				{
					vector<uint64_t> vv;
					vv.push_back(v.monotig_rank);
					pair<count_vector, vector<uint64_t>> p (v,vv);
					bucket_class.insert({v.compressed,p});
				}
				else
				{
					// for count vectors identical to a representant, rememember their index is associated to this representant
					bucket_class[v.compressed].second.push_back(v.monotig_rank);
				}
				in.peek();
			}
			mm.lock();
			if (do_query_on_disk)
			{
				//~ get_eq_classes_disk_query(output, bucket_class,  nb_unitigs, nb_colors, final_positions, prev_pos, out);
				get_eq_classes_disk_query(reindeer_index_files, bucket_class,  nb_monotig, nb_colors, final_positions, prev_pos, out);
			}
			else
			{
				//~ get_eq_classes(output, bucket_class,  nb_unitigs, nb_colors, final_positions, prev_pos, outz);
				get_eq_classes(reindeer_index_files, bucket_class,  nb_monotig, nb_colors, final_positions, prev_pos, out);
			}
			mm.unlock();
		} 
			
		//~ in.close();
	}
	long nb_eq_class(prev_pos + 1);
	cout << "Number of equivalence classes found: " << prev_pos + 1<< endl;
	//~ zstr::ofstream * out_position =  new zstr::ofstream(output + "/reindeer_matrix_eqc_position.gz");
	zstr::ofstream * out_position =  new zstr::ofstream(matrix_eqc_position_file);
	uint nb(0);
	for (uint i(0); i < final_positions.size(); ++i)
	{
		out_position->write(reinterpret_cast<char*>(&final_positions[i]), sizeof(long)); 
		++nb;
	}
	delete out_position;
	out_info->write(reinterpret_cast<char*>(&nb_eq_class), sizeof(long));
	out_info->write(reinterpret_cast<char*>(&nb_colors), sizeof(uint64_t));
	
	//~ if (do_query_on_disk)
		delete out;
	//~ else
		//~ delete outz;
}



