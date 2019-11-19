#ifndef MAT
#define MAT


#include "../trle/trle.h"


using namespace std;


// load rle encoded matrix from disk (keep compressed in ram)
vector<unsigned char*> load_compressed_vectors(const string& input_file, vector<unsigned>&vector_sizes, uint64_t& color_number, uint64_t& minitig_number)
{
	ifstream in(input_file);
	if(not exists_test(input_file))
	{
		cout<<"File problem"<<endl;
		exit(0);
	}
	int64_t rank;
	unsigned line_size;
	in.read(reinterpret_cast<char *>(&minitig_number), sizeof(uint64_t));
	in.read(reinterpret_cast<char *>(&color_number), sizeof(uint64_t));
	vector_sizes.resize(minitig_number, 0);
	vector<unsigned char*> out(minitig_number);
	for (uint minit(0); minit < out.size(); ++minit)
	{
		in.read(reinterpret_cast<char *>(&line_size), sizeof(unsigned));
		in.read(reinterpret_cast<char *>(&rank), sizeof(int64_t));
		vector_sizes[rank] = line_size; 
		out[rank] = new unsigned char [line_size];
		
		in.read((char*)((out[rank])) , line_size);
	}
	return out;
}


// dump rle vector on disk
void dump_compressed_vector(vector<uint16_t>& counts, int64_t minitig_id, ofstream& out, unsigned char *in)
{
	// one vector corresponding to the minitig count/colors
	// convert to string for the compression
	uint n(counts.size()*2);
	uint nn(counts.size()*2 + 1024);
	unsigned char comp[nn];
	in = (unsigned char*)&counts[0];
	unsigned compr_vector_size = trlec(in, n, comp) ; //todo can I write it now on the disk
	out.write(reinterpret_cast<char*>(&compr_vector_size),sizeof(unsigned));
	int64_t tw(minitig_id);
	out.write(reinterpret_cast<char*>(&tw),sizeof(int64_t));
	out.write((const char*)comp,(compr_vector_size));
}






#endif
