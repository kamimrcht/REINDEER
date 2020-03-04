#ifndef MAT
#define MAT


#include "../trle/trle.h"
#include "query.hpp"
#include "utils.hpp"

using namespace std;


// dump rle vector on buffer
void dump_compressed_vector_buff(vector<uint16_t>& counts, int64_t minitig_id, string& buffer, unsigned char *in)
{
	// one vector corresponding to the minitig count/colors
	// convert to string for the compression
	uint n(counts.size()*2);
	uint nn(counts.size()*2 + 1024);
	unsigned char comp[nn];
	in = (unsigned char*)&counts[0];
	unsigned compr_vector_size = trlec(in, n, comp) ; //todo can I write it now on the disk
	buffer.append(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned));
	int64_t tw(minitig_id);
	buffer.append(reinterpret_cast<char*>(&tw),sizeof(int64_t));
	buffer.append((const char*)comp,(compr_vector_size));
}

// load rle encoded matrix from disk (keep compressed in ram)
vector<unsigned char*> load_compressed_vectors(const string& input_file, vector<unsigned>&vector_sizes, uint64_t& color_number, uint64_t& minitig_number, long eq_class_nb)
{
	ifstream in_nb(input_file + "_minitig_nb");
	auto in = new zstr::ifstream(input_file);
	if(not exists_test(input_file))
	{
		cout<<"File problem"<<endl;
		exit(0);
	}
	int64_t rank;
	unsigned line_size;
	in_nb.read(reinterpret_cast<char *>(&minitig_number), sizeof(uint64_t));
	in->read(reinterpret_cast<char *>(&color_number), sizeof(uint64_t));
	vector<unsigned char*> out; //todo resize when loading eq class
	uint t(0);
	for (uint i(0); i < eq_class_nb; ++i) //lire le nb de classes d'eq
	{
		in->read(reinterpret_cast<char *>(&line_size), sizeof(unsigned));
		in->read(reinterpret_cast<char *>(&rank), sizeof(int64_t));
		vector_sizes.push_back(line_size); 
		unsigned char* c = new unsigned char[line_size];
		in->read((char*)(c) , line_size);
		out.push_back(c);
		++t;
	}
	in_nb.close();
	delete in;
	return out;
}


void dump_compressed_vector_bucket_disk_query(vector<uint16_t>& counts, int64_t minitig_id, unsigned char *in,vector<ofstream*>& bucket_files,  vector<uint8_t>& colors, bool record_counts)
{

	vector<unsigned char> comp;
	unsigned compr_vector_size;
	int64_t tw(minitig_id);
	if (record_counts)
	{
		uint n(counts.size()*2);
		uint nn(counts.size()*2 + 1024);
		comp = RLE16C(counts); //homemade RLE with escape character 255
	}
	else
	{
		comp = RLE8C(colors);
	}
	compr_vector_size = comp.size();
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
}

void dump_compressed_vector_bucket(vector<uint16_t>& counts, int64_t minitig_id, unsigned char *in,  vector<ofstream*>& bucket_files, vector<uint8_t>& colors, bool record_counts )
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
}





void read_matrix_compressed_line(ifstream& in, int64_t& rank, char* comp, unsigned& comp_size)
{
	in.read(reinterpret_cast<char *>(&comp_size), sizeof(unsigned));
	in.read(reinterpret_cast<char *>(&rank), sizeof(int64_t));
	in.read((char*)(comp) , comp_size);

}


#endif
