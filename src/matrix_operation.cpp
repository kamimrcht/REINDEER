#include "matrix_operation.hpp"
using namespace std;

// dump rle vector on buffer
void dump_compressed_vector_buff(vector<uint16_t>& counts, int64_t monotig_id, string& buffer, unsigned char* in)
{
    // one vector corresponding to the monotig count/colors
    // convert to string for the compression
    uint n(counts.size() * 2);
    uint nn(counts.size() * 2 + 1024);
    unsigned char comp[nn];
    in = (unsigned char*)&counts[0];
    unsigned compr_vector_size = trlec(in, n, comp); //todo can I write it now on the disk
    buffer.append(reinterpret_cast<char*>(&compr_vector_size), sizeof(unsigned));
    int64_t tw(monotig_id);
    buffer.append(reinterpret_cast<char*>(&tw), sizeof(int64_t));
    buffer.append((const char*)comp, (compr_vector_size));
}

// load rle encoded matrix from disk (keep compressed in ram)
vector<unsigned char*> load_compressed_vectors(const string& input_file, vector<unsigned>& vector_sizes, uint64_t& color_number, uint64_t& monotig_number, long eq_class_nb)
{
    //~ ifstream in_nb(input_file + "_monotig_nb");
    auto in = new zstr::ifstream(input_file);
    if (not exists_test(input_file)) {
        cout << "File problem" << endl;
        exit(0);
    }
    int64_t rank;
    unsigned line_size;
    in->read(reinterpret_cast<char*>(&color_number), sizeof(uint64_t));
    vector<unsigned char*> out; //todo resize when loading eq class
    uint t(0);
    for (uint i(0); i < eq_class_nb; ++i) //lire le nb de classes d'eq
    {
        in->read(reinterpret_cast<char*>(&line_size), sizeof(unsigned));
        in->read(reinterpret_cast<char*>(&rank), sizeof(int64_t));
        vector_sizes.push_back(line_size);
        unsigned char* c = new unsigned char[line_size];
        in->read((char*)(c), line_size);
        out.push_back(c);
        ++t;
    }
    delete in;
    return out;
}

void dump_compressed_vector_bucket_disk_query(vector<uint16_t>& counts, int64_t monotig_id, unsigned char* in, vector<ofstream*>& bucket_files, vector<uint8_t>& colors, bool record_counts)
{
    vector<unsigned char> comp;
    unsigned compr_vector_size;
    int64_t tw(monotig_id);
    if (record_counts) {
        uint n(counts.size() * 2);
        uint nn(counts.size() * 2 + 1024);
        comp = RLE16C(counts); //homemade RLE with escape character 255
        vector<uint16_t> decomp = RLE16D(comp);

    } else {
        comp = RLE8C(colors);
    }
    compr_vector_size = comp.size();
    uint32_t bucket_nb;
    // hash bit from the compressed counts to choose the bucket
    if (compr_vector_size >= 8) {
        bucket_nb = (xorshift((uint64_t)comp[0] + (uint64_t)comp[1] * 256 + (uint64_t)comp[2] * 256 * 256 + (uint64_t)comp[3] * 256 * 256 * 256 + (uint64_t)comp[4] * 256 * 256 * 256 * 256 + (uint64_t)comp[5] * 256 * 256 * 256 * 256 * 256 + (uint64_t)comp[6] * 256 * 256 * 256 * 256 * 256 * 256 + (uint64_t)comp[7] * 256 * 256 * 256 * 256 * 256 * 256 * 256) % bucket_files.size());
    } else if (compr_vector_size >= 4) {
        bucket_nb = (xorshift((uint64_t)comp[0] + (uint64_t)comp[1] * 256 + (uint64_t)comp[2] * 256 * 256 + (uint64_t)comp[3] * 256 * 256 * 256) % bucket_files.size());
    } else {
        bucket_nb = (xorshift((uint64_t)comp[0]) % bucket_files.size());
    }
    //write info in the bucket
    long position(bucket_files[bucket_nb]->tellp());
    bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&compr_vector_size), sizeof(unsigned)); //size of the compressed information
    bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&tw), sizeof(int64_t)); //index
    bucket_files[bucket_nb]->write((const char*)&comp[0], (compr_vector_size)); // compressed count vector
}

//write count vectors in buckets to compare them and de-duplicate the count matrix
void dump_compressed_vector_bucket(int64_t monotig_id, vector<ofstream*>& bucket_files, string& header)
{
    //get the bucket number
    uint32_t bucket_nb;
    // hash bit from the compressed counts to choose the bucket
    if (header.size() >= 12) {
        bucket_nb = (xorshift((uint8_t)(header[4]) + (uint8_t)(header[5]) * 256 + (uint8_t)(header[6]) * 256 * 256 + (uint8_t)(header[7]) * 256 * 256 * 256 + (uint8_t)(header[8]) * 256 * 256 * 256 * 256 + (uint8_t)(header[9]) * 256 * 256 * 256 * 256 * 256 + (uint8_t)(header[10]) * 256 * 256 * 256 * 256 * 256 * 256 + (uint8_t)(header[11]) * 256 * 256 * 256 * 256 * 256 * 256 * 256) % bucket_files.size());

    } else if (header.size() >= 8) {
        //~ #pragma GCC diagnostic ignored "-fpermissive"
        bucket_nb = (xorshift((uint8_t)(header[4]) + (uint8_t)(header[5]) * 256 + (uint8_t)(header[6]) * 256 * 256 + (uint8_t)(header[7]) * 256 * 256 * 256) % bucket_files.size());
    } else {
        bucket_nb = (xorshift((uint8_t)(header[4])) % bucket_files.size());
    }
    //write info in the bucket
    bucket_files[bucket_nb]->write((const char*)(&header[0]), header.size()); // size of compressed rle + rle
    bucket_files[bucket_nb]->write(reinterpret_cast<char*>(&monotig_id), sizeof(int64_t)); //index
}

void read_matrix_compressed_line(zstr::ifstream& in, int64_t& rank, char* comp, unsigned& comp_size)
{
    in.read(reinterpret_cast<char*>(&comp_size), sizeof(unsigned));
    in.read((char*)(comp), comp_size);
    in.read(reinterpret_cast<char*>(&rank), sizeof(int64_t));
}
