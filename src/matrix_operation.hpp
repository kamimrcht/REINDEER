#ifndef MAT
#define MAT


#include "../trle/trle.h"
#include "query.hpp"

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
	//~ in->close();
	delete in;
	return out;
}


// dump rle vector on disk
void dump_compressed_vector(vector<uint16_t>& counts, int64_t minitig_id, ofstream& out, unsigned char *in, ofstream& out_positions)
{
	// one vector corresponding to the minitig count/colors
	// convert to string for the compression
	uint n(counts.size()*2);
	uint nn(counts.size()*2 + 1024);
	vector<unsigned char> comp = RLE16C(counts); //homemade RLE with escape character 255
	unsigned compr_vector_size = comp.size();
	long position(out.tellp());
	out.write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned));
	int64_t tw(minitig_id);
	out.write(reinterpret_cast<char*>(&tw),sizeof(int64_t));
	out.write((const char*)&comp[0],(compr_vector_size));
	char* returnc = new char[1];
	returnc[0] = (uint8_t) 255;
	out.write(&returnc[0], sizeof(char)); // line separator
	out_positions.write(reinterpret_cast<char*>(&position), sizeof(long));
	delete [] returnc;
}






#endif
