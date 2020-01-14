#ifndef EQ
#define EQ


#include "../trle/trle.h"
#include "../blight/utils.h"

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


//~ uint32_t hashnadine ( uint32_t x ) {
	//~ x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	//~ x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	//~ x = ( ( x >> 16 ) ^ x );
	//~ return x;
//~ }

  uint64_t hashnadine ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}

void dump_compressed_vector_bucket(vector<uint16_t>& counts, int64_t minitig_id, unsigned char *in, ofstream& out_positions, vector<ofstream*>& bucket_files)
{
	//get the bucket number
	uint n(counts.size()*2);
	uint nn(counts.size()*2 + 1024);
	unsigned char comp[nn];
	in = (unsigned char*)&counts[0];
	unsigned compr_vector_size = trlec(in, n, comp) ;
	uint64_t bucket_nb;
	if (compr_vector_size >= 8)
		bucket_nb = (hashnadine((uint64_t)comp[0]) % bucket_files.size());
	else
		bucket_nb = (hashnadine((uint8_t)comp[0]) % bucket_files.size());

	//~ uint32_t bucket_nb (hashnadine((uint32_t)comp[0] + (uint32_t)comp[1]*256) % bucket_files.size());
		//~ cout << "here "<< " " << bucket_nb << endl;

	//write info in the bucket
	long position(bucket_files[bucket_nb]->tellp());
	bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned));
	int64_t tw(minitig_id);
	bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&tw),sizeof(int64_t));
	bucket_files[bucket_nb]->write((const char*)comp,(compr_vector_size));
	out_positions.write(reinterpret_cast<char*>(&position), sizeof(long));
}





void read_matrix_compressed_line(ifstream& in, int64_t& rank, char* comp, unsigned& comp_size)
{
	in.read(reinterpret_cast<char *>(&comp_size), sizeof(unsigned));
	in.read(reinterpret_cast<char *>(&rank), sizeof(int64_t));
	in.read((char*)(comp) , comp_size);

}


void get_eq_classes(string& output, vector<count_vector>& count_vecs, uint64_t unitig_nb, uint64_t color_number, vector<long>& final_positions, long& prev_pos, ofstream& out)
{
	char* comp;
	string prev_comp("");
	

	//todo parallel
	//~ long position_in_outfile;
	//~ vector<long> final_positions(unitig_nb);
	for (uint i(0); i < count_vecs.size(); ++i)
	{
		if (count_vecs[i].compressed == prev_comp)
		{
			final_positions[count_vecs[i].minitig_rank] = prev_pos;
		}
		else
		{
			++prev_pos;
			prev_comp = count_vecs[i].compressed;
			//write the new vector in the matrix
			//~ position_in_outfile = out.tellp();
			//~ prev_pos = out.tellp();
			final_positions[count_vecs[i].minitig_rank] = prev_pos;
			char* comp (&count_vecs[i].compressed[0]);
			out.write(reinterpret_cast<char*>(&count_vecs[i].compressed_size),sizeof(unsigned));
			out.write(reinterpret_cast<char*>(&count_vecs[i].minitig_rank),sizeof(int64_t));
			out.write((const char*)comp,(count_vecs[i].compressed_size));
			
		}
	}
	//~ cin.get();
	//write a position file
	//~ ofstream out_position(output + "/reindeer_matrix_eqc_position");
	//~ uint nb(0);
	//~ for (uint i(0); i < final_positions.size(); ++i)
	//~ {
		//~ out_position.write(reinterpret_cast<char*>(&final_positions[i]), sizeof(long));
		//~ if (final_positions[i] != 0)
		//~ cout << final_positions[i] << endl;
		//~ ++nb;
	//~ }
	//~ out_position.close();
	//~ cout << "nb positions " << nb << endl;
}


//sort count vectors by file, write one occurence per count in a new matrix file
void write_eq_class_matrix(string& output, vector<ofstream*>& all_files, uint64_t nb_unitigs, uint64_t color_number)
{
	cout << "Sorting datasets to find equivalence classes..." << endl;
	vector<long> final_positions(nb_unitigs);
	//todo parallel
	ofstream out(output + "/reindeer_matrix_eqc");
	out.write(reinterpret_cast<char*>(&color_number),sizeof(uint64_t)); // number of colors
	long prev_pos(-1);
	uint i(0);
	mutex mm; 
	#pragma omp parallel for num_threads(4)
	for (i=0; i < all_files.size(); ++i)
	{
		int64_t rank;
		uint sizecp(color_number*2 + 1024);
		char comp[sizecp];
		unsigned comp_size;
		vector<count_vector> count_vecs;
		string compressed;
		ifstream in(output + "/matrix_bucket_"+ to_string(i)); //TODO zstr??
		if (not is_empty_file(in)){
			while (not in.eof())
			{
				
				//read file and store each vector in struct
				read_matrix_compressed_line(in, rank, comp, comp_size);
				//~ unsigned char* lo = decode_vector((unsigned char*)comp, comp_size, color_number);
				//~ vector<uint16_t> cc = count_string_to_count_vector(lo, comp_size);
				//~ cout << "*************" << cc[0] << " " << cc[1] << endl;
				compressed.assign(comp, comp_size);
				count_vector v({comp_size, rank, compressed});
				//~ cout << v.compressed << endl;
				unsigned char* lo = decode_vector((unsigned char*)&v.compressed[0], v.compressed_size, color_number);
				vector<uint16_t> cc = count_string_to_count_vector(lo, v.compressed_size);
				count_vecs.push_back(v);
			}
			
			//sort compressed counts
			sort_vectors(count_vecs);
			// go through sorted vectors and write the final matrix
			mm.lock();
			get_eq_classes(output, count_vecs,  nb_unitigs, color_number, final_positions, prev_pos, out);
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
		//~ out_position.write(reinterpret_cast<char*>(&final_positions[i]), sizeof(long)); // ???
		out_position.write(reinterpret_cast<char*>(&final_positions[i]), sizeof(long)); // ???
		++nb;
	}
	out_position.close();
	out.close();
	ofstream out_nbc(output + "/reindeer_matrix_eqc_nb_class");
	out_nbc.write(reinterpret_cast<char*>(&nb_eq_class), sizeof(long));
	out_nbc.close();
}


#endif
