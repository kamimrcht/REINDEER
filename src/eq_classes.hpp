#ifndef EQ
#define EQ


#include "../trle/trle.h"
#include "../blight/utils.h"
#include "../blight/robin_hood.h"

using namespace std;


struct count_vector
{
	unsigned compressed_size;
	int64_t minitig_rank;
	string compressed;
};


struct compare_vec
{
	bool operator()(count_vector& vec1, count_vector& vec2)
	{
		return vec1.compressed < vec2.compressed;
	}
};


void sort_vectors(vector<count_vector>& matrix_lines)
{
	sort(matrix_lines.begin(), matrix_lines.end(), compare_vec());
}


//~ uint32_t xorshift ( uint32_t x ) {
	//~ x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	//~ x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	//~ x = ( ( x >> 16 ) ^ x );
	//~ return x;
//~ }

  uint64_t xorshift ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}


void dump_compressed_vector_bucket_disk_query(vector<uint16_t>& counts, int64_t minitig_id, unsigned char *in, ofstream& out_positions, vector<ofstream*>& bucket_files)
{

	uint n(counts.size()*2);
	uint nn(counts.size()*2 + 1024);
	vector<unsigned char> comp;
	//~ in = (unsigned char*)&counts[0];
	unsigned compr_vector_size;
	int64_t tw(minitig_id);

	comp = RLE16C(counts); //homemade RLE with escape character 255
	compr_vector_size = comp.size();
	cout << comp.size() << endl;
	uint32_t bucket_nb;
	// hash bit from the compressed counts to choose the bucket
	if (compr_vector_size >= 8){
		bucket_nb = (xorshift((uint64_t)comp[0] + (uint64_t)comp[1]*256 +  (uint64_t)comp[2]*256*256 + (uint64_t)comp[3]*256*256*256 + (uint64_t)comp[4]*256*256*256*256 + (uint64_t)comp[5]*256*256*256*256*256 + (uint64_t)comp[6]*256*256*256*256*256*256 + (uint64_t)comp[7]*256*256*256*256*256*256*256 ) % bucket_files.size());
	
	}else if (compr_vector_size >= 4){
		bucket_nb = (xorshift((uint64_t)comp[0] + (uint64_t)comp[1]*256 +  (uint64_t)comp[2]*256*256 + (uint64_t)comp[3]*256*256*256  ) % bucket_files.size());
	}else{
		bucket_nb = (xorshift((uint64_t)comp[0]) % bucket_files.size());
	}
	//write info in the bucket
	long position(bucket_files[bucket_nb]->tellp());
	bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned)); //size of the compressed information
	bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&tw),sizeof(int64_t)); //index
	bucket_files[bucket_nb]->write((const char*)&comp[0],(compr_vector_size)); // compressed count vector	
	out_positions.write(reinterpret_cast<char*>(&position), sizeof(long)); //position in vector of count vectors
}

void dump_compressed_vector_bucket(vector<uint16_t>& counts, int64_t minitig_id, unsigned char *in, ofstream& out_positions, vector<ofstream*>& bucket_files, vector<uint8_t>& colors, bool record_counts )
{
	//get the bucket number
	uint n;
	if (record_counts)
		n = counts.size()*2;
	else
		n = colors.size();
	uint nn(n  + 1024);
	unsigned char comp[nn];
	if (record_counts)
		in = (unsigned char*)&counts[0];
	else
		in = (unsigned char*)&colors[0];
	unsigned compr_vector_size;
	compr_vector_size = trlec(in, n, comp) ;

	uint32_t bucket_nb;
	// hash bit from the compressed counts to choose the bucket
	if (compr_vector_size >= 8){
		bucket_nb = (xorshift((uint64_t)comp[0] + (uint64_t)comp[1]*256 +  (uint64_t)comp[2]*256*256 + (uint64_t)comp[3]*256*256*256 + (uint64_t)comp[4]*256*256*256*256 + (uint64_t)comp[5]*256*256*256*256*256 + (uint64_t)comp[6]*256*256*256*256*256*256 + (uint64_t)comp[7]*256*256*256*256*256*256*256 ) % bucket_files.size());
	
	}else if (compr_vector_size >= 4){
		bucket_nb = (xorshift((uint64_t)comp[0] + (uint64_t)comp[1]*256 +  (uint64_t)comp[2]*256*256 + (uint64_t)comp[3]*256*256*256  ) % bucket_files.size());
	}else{
		bucket_nb = (xorshift((uint64_t)comp[0]) % bucket_files.size());
	}
	//write info in the bucket
	long position(bucket_files[bucket_nb]->tellp());
	bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned)); //size of the compressed information
	bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&minitig_id),sizeof(int64_t)); //index
	bucket_files[bucket_nb]->write((const char*)comp,(compr_vector_size)); // compressed count vector
	out_positions.write(reinterpret_cast<char*>(&position), sizeof(long)); //position in vector of count vectors
}





void read_matrix_compressed_line(ifstream& in, int64_t& rank, char* comp, unsigned& comp_size)
{
	in.read(reinterpret_cast<char *>(&comp_size), sizeof(unsigned));
	in.read(reinterpret_cast<char *>(&rank), sizeof(int64_t));
	in.read((char*)(comp) , comp_size);

}


// write final matrix of equivalence classes
void get_eq_classes(string& output, robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>>& bucket_class, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, zstr::ofstream* out)
{
	string prev_comp("");
	// for each representant, write the compressed size, index and compressed counts on disk
	// write a file:  integer i at position n shows which row i of the matrix should be looked up for minitig n (surjection minitigs -> class representants)
	for (auto rk(bucket_class.begin()); rk != bucket_class.end(); ++rk)
	{
		++prev_pos;
		count_vector v (rk->second.first);
		char* comp (&v.compressed[0]);
		out->write(reinterpret_cast<char*>(&v.compressed_size),sizeof(unsigned));
		out->write(reinterpret_cast<char*>(&v.minitig_rank),sizeof(int64_t));
		out->write((const char*)comp,(v.compressed_size));
		for (auto && rank: rk->second.second)
		{
			if (rank < final_positions.size())
				final_positions[rank] = prev_pos;
		}
	}
}


void get_eq_classes_disk_query(string& output, robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>>& bucket_class, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, zstr::ofstream* out)
{
	string prev_comp("");
	// for each representant, write the compressed size, index and compressed counts on disk
	// write a file:  integer i at position n shows which row i of the matrix should be looked up for minitig n (surjection minitigs -> class representants)
	for (auto rk(bucket_class.begin()); rk != bucket_class.end(); ++rk)
	{
		//~ ++prev_pos;
		count_vector v (rk->second.first);
		char* comp (&v.compressed[0]); //homemade RLE compression
		//~ char* decomp = new char[color_number * 2 +1024];
		//~ auto sz = trled(comp, v.compressed.size() , decomp, color_number * 2);
		// 1 -decompress to compress with homemade RLE
		// 2-compress
		//~ vector<unsigned char> comp2 = RLE16C(counts); //homemade RLE with escape character 255 // TODO change and use homemade RLE before if disk query

		unsigned compr_vector_size = v.compressed.size();
		long position(out->tellp());
		out->write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned));
		int64_t tw(v.minitig_rank);
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
void write_eq_class_matrix(string& output, vector<ofstream*>& all_files, uint64_t nb_unitigs, uint64_t color_number, bool do_query_on_disk)
{
	cout << "Sorting datasets to find equivalence classes..." << endl;
	vector<long> final_positions(nb_unitigs);
	//todo parallel
	//~ ofstream out(output + "/reindeer_matrix_eqc");
	auto out = new zstr::ofstream(output + "/reindeer_matrix_eqc");
	out->write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	long prev_pos(-1);
	uint i(0);
	mutex mm; 
	//~ #pragma omp parallel for num_threads(4)
	// in each bucket file
	for (i=0; i < all_files.size(); ++i)
	{
		int64_t rank;
		uint sizecp(color_number*2 + 1024);
		char comp[sizecp];
		unsigned comp_size;
		vector<count_vector> count_vecs;
		robin_hood::unordered_map<string, pair<count_vector, vector<uint64_t>>> bucket_class;
		string compressed;
		ifstream in(output + "/matrix_bucket_"+ to_string(i)); //TODO zstr??
		if (not is_empty_file(in)){
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
					vv.push_back(v.minitig_rank);
					pair<count_vector, vector<uint64_t>> p (v,vv);
					bucket_class.insert({v.compressed,p});
				}
				else
				{
					// for count vectors identical to a representant, rememember their index is associated to this representant
					bucket_class[v.compressed].second.push_back(v.minitig_rank);
				}
				in.peek();

			}
			
			mm.lock();
			if (do_query_on_disk)
				get_eq_classes_disk_query(output, bucket_class,  nb_unitigs, color_number, final_positions, prev_pos, out);
			else
				get_eq_classes(output, bucket_class,  nb_unitigs, color_number, final_positions, prev_pos, out);
			mm.unlock();
		}
		in.close();
	}
	long nb_eq_class(prev_pos + 1);
	cout << "Number of equivalence classes found: " << prev_pos + 1<< endl;
	ofstream out_position(output + "/reindeer_matrix_eqc_position");
	uint nb(0);
	for (uint i(0); i < final_positions.size(); ++i)
	{
		out_position.write(reinterpret_cast<char*>(&final_positions[i]), sizeof(long)); 
		++nb;
	}
	out_position.close();
	//~ out.close();
	//~ out->close();
	delete out;
	ofstream out_nbc(output + "/reindeer_matrix_eqc_nb_class");
	out_nbc.write(reinterpret_cast<char*>(&nb_eq_class), sizeof(long));
	out_nbc.close();
}


#endif
