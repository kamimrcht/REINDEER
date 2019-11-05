#ifndef KSL
#define KSL



#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "common.h"
#include "bbhash.h"
#include "bm64.h"
#include "encoding.h"
#include "zstr.hpp"
#include "bmserial.h"



using namespace std;


// FOR k<32
#define kmer uint64_t
// FOR k<64
//~ #define kmer __uint128_t




#define minimizer_type uint32_t




//~ typedef SingleHashFunctor128 hasher_t;
typedef boomphf::SingleHashFunctor<kmer>  hasher_t;



//~ typedef boomphf::mphf<  hasher_t,kmer  > MPHF;
typedef boomphf::mphf<  kmer, hasher_t  > MPHF;



struct bucket_minimizer{
	uint64_t current_pos;
	uint64_t start;
	//~ uint32_t abundance_minimizer;
	uint32_t nuc_minimizer;
	uint32_t skmer_number;
};



struct info_mphf{
	uint64_t mphf_size;
	uint64_t start;
	MPHF* kmer_MPHF;
	uint8_t bit_to_encode;
	bool empty;
};



// Represents the cardinality of a pow2 sized set. Allows div/mod arithmetic operations on indexes.
template<typename T>
struct Pow2 {
	Pow2(uint_fast8_t bits) : _bits(bits) {
		assume(bits < CHAR_BIT*sizeof(T), "Pow2(%u > %u)", unsigned(bits), unsigned(CHAR_BIT*sizeof(T)));
	}

	uint_fast8_t bits() const { return _bits; }
	T value() const { return T(1) << _bits; }
	explicit operator T() const { return value(); }
	T max() const { return value() - T(1); }
	Pow2& operator=(const Pow2&) = default;

	Pow2():_bits(0){}
	Pow2(const Pow2&) = default;

	friend T operator*(const T& x, const Pow2& y) { return x << y._bits; }
	friend T& operator*=(T& x, const Pow2& y) { return x <<= y._bits; }
	friend T operator/(const T& x, const Pow2& y) { return x >> y._bits; }
	friend T& operator/=(T& x, const Pow2& y) { return x >>= y._bits; }
	friend T operator%(const T& x, const Pow2& y) { return x & y.max(); }
	friend T& operator%=(T& x, const Pow2& y) { return x &= y.max(); }
	Pow2& operator>>=(uint_fast8_t d) { _bits -= d; return *this; }
	Pow2& operator<<=(uint_fast8_t d) { _bits += d; return *this; }
	friend bool operator<(const T& x, const Pow2& y) { return x < y.value(); }
	friend bool operator<=(const T& x, const Pow2& y) { return x < y.value(); }
	friend T operator+(const T& x, const Pow2& y) { return x + y.value(); }
	friend T& operator+=(T& x, const Pow2& y) { return x += y.value(); }
	friend T operator-(const T& x, const Pow2& y) { return x - y.value(); }
	friend T& operator-=(T& x, const Pow2& y) { return x -= y.value(); }
private:
	uint_fast8_t _bits;
};



class kmer_Set_Light{
public:
	uint64_t k, m1, m2, m3, minimizer_size_graph;
	uint64_t coreNumber;
	uint64_t bit_saved_sub;

	Pow2<kmer> offsetUpdateAnchor;
	Pow2<kmer> offsetUpdateMinimizer;
	Pow2<uint> mphf_number;
	Pow2<uint> number_superbuckets;
	Pow2<uint> minimizer_number;
	Pow2<uint> minimizer_number_graph;
	Pow2<uint> number_bucket_per_mphf;
	Pow2<uint> bucket_per_superBuckets;
	Pow2<uint> positions_to_check;

	vector<bool> bucketSeq;
	vector<bool> positions;
	bm::bvector<> position_super_kmers= bm::bvector<>(bm::BM_GAP);
	bm::bvector<>::rs_index_type* position_super_kmers_RS;
	mutex positions_mutex[4096];

	bucket_minimizer* all_buckets;
	uint32_t* abundance_minimizer_temp;
	uint8_t* abundance_minimizer;
	info_mphf* all_mphf;

	uint64_t number_kmer = 0;
	uint64_t number_super_kmer = 0;
	uint64_t largest_MPHF = 0;
	uint64_t positions_total_size = 0;
	uint64_t positions_int = 0;
	atomic<uint64_t> number_query;

	uint64_t total_nb_minitigs = 0;
	//~ uint64_t hell_bucket;
	double bit_per_kmer = 0;
	uint32_t largest_bucket_nuc_all = 0;
	uint64_t gammaFactor=2;
	kmer_Set_Light(uint64_t k_val,uint64_t m1_val, uint64_t m2_val, uint64_t m3_val, uint64_t coreNumber_val, uint64_t bit_to_save)
		: k(k_val)
		, m1(m1_val)
		, m2(m1_val)
		, m3(m3_val)
		, minimizer_size_graph(m2_val < m1 ? max(m1+2,(uint64_t)11) : m2_val)
		, coreNumber(coreNumber_val)
		, bit_saved_sub(bit_to_save)

		, offsetUpdateAnchor(2*k)
		, offsetUpdateMinimizer(2*m1)
		, mphf_number(2*m2)
		, number_superbuckets(2*m3)
		, minimizer_number(2*m1)
		, minimizer_number_graph(2*minimizer_size_graph)
		, number_bucket_per_mphf(2*(m1-m2))
		, bucket_per_superBuckets(2*(m1-m3))
		, positions_to_check(bit_to_save)
	{
		all_buckets=new bucket_minimizer[minimizer_number.value()]();
		all_mphf=new info_mphf[mphf_number.value()];
		for(uint64_t i(0);i<mphf_number;++i){
			all_mphf[i].mphf_size=0;
			all_mphf[i].bit_to_encode=0;
			all_mphf[i].start=0;
			all_mphf[i].empty=true;

		}
		for(uint64_t i(0);i<minimizer_number.value();++i){
			all_buckets[i].current_pos=0;
			all_buckets[i].start=0;
			all_buckets[i].nuc_minimizer=0;
			all_buckets[i].skmer_number=0;
		}
		number_query=0;
		struct rlimit rl;
		getrlimit (RLIMIT_NOFILE, &rl);
		rl.rlim_cur = number_superbuckets.value()+5;
		setrlimit (RLIMIT_NOFILE, &rl);
	}

	kmer_Set_Light(const string& index_file);

	~kmer_Set_Light () {
		delete position_super_kmers_RS;
		delete[] all_buckets;
		for(uint64_t i(0);i<mphf_number.value();++i){
			if( not all_mphf[i].empty){
				delete all_mphf[i].kmer_MPHF;
			}
		}
		delete[] all_mphf;
		struct stat buffer;
		string filerm;
		for(uint64_t i(0);i<number_superbuckets.value();++i){
			filerm=("_out"+to_string(i));
			if((stat (filerm.c_str(), &buffer) == 0)){
				remove(filerm.c_str());
			}
		}
	}

	bool exists(const kmer& query);
	void create_super_buckets(const string& input_file);
	void read_super_buckets(const string& input_file);
	void create_mphf_mem(uint64_t beg,uint64_t end);
	void create_mphf_disk(uint64_t beg,uint64_t end);
	void updateK(kmer& min, char nuc);
	void updateRCK(kmer& min, char nuc);
	void updateM(kmer& min, char nuc);
	void updateRCM(kmer& min, char nuc);
	void fill_positions(uint64_t beg,uint64_t end);
	bool exists(const string& query);
	void multiple_query(const string& query);
	kmer minimizer_according_xs(kmer seq);
	void abundance_minimizer_construct(const string& input_file);
	int64_t correct_pos(kmer mini, uint64_t p);
	void str2bool(const string& str,uint64_t mini);
	kmer update_kmer(uint64_t pos,kmer mini,kmer input);
	kmer get_kmer(uint64_t pos,uint64_t mini);
	void print_kmer(kmer num,uint64_t n=100);
	int32_t query_get_pos_unitig(const kmer canon,kmer minimizer);
	uint64_t multiple_query_serial(const uint64_t minimizerV, const vector<kmer>& kmerV);
	void file_query(const string& query_file);
	uint64_t bool_to_int(uint64_t n_bits_to_encode,uint64_t pos,uint64_t start);
	uint64_t multiple_query_optimized(kmer minimizerV, const vector<kmer>& kmerV);
	void int_to_bool(uint64_t n_bits_to_encode,uint64_t X, uint64_t pos,uint64_t start);
	kmer update_kmer_local(uint64_t pos,const vector<bool>& V,kmer input);
	vector<bool> get_seq(kmer mini,uint64_t pos,uint64_t n);
	kmer minimizer_graph(kmer seq);
	bool single_query(const kmer minimizer, kmer kastor);
	bool multiple_minimizer_query_bool(const kmer minimizer,  kmer kastor,uint64_t prefix_length,uint64_t suffix_length);
	int64_t multiple_minimizer_query_hash(const kmer minimizer,  kmer kastor,uint64_t prefix_length,uint64_t suffix_length);
	bool query_kmer_bool(kmer canon);
	pair<uint64_t,uint64_t> query_sequence_bool(const string& query);
	string kmer2str(kmer num);
	kmer regular_minimizer(kmer seq);
	void create_super_buckets_regular(const string&,int dbg_id=0);
	int64_t query_kmer_hash(kmer canon);
	int64_t query_get_hash(const kmer canon,kmer minimizer);
	vector<int64_t> query_sequence_hash(const string& query);
	void construct_index(const string& input_file);
	//~ void report_memusage(boomphf::memreport_t& report, const std::string& prefix="blight", bool add_struct=true);
	vector<int64_t> query_sequence_minitig(const string& query);
	int64_t query_get_rank_minitig(const kmer canon,uint64_t minimizer);
	int64_t query_kmer_minitig(kmer canon);
	kmer mantis(uint64_t n);
	vector<kmer> kmer_to_superkmer(const kmer canon,kmer minimizer, int64_t& rank, int64_t& hash);
	int64_t hash_to_rank(const int64_t hash,kmer minimizer);
	int64_t kmer_to_hash(const kmer canon,kmer minimizer);
	void get_monocolor_minitigs(const  vector<string>& minitigs, const vector<int64_t>& color,zstr::ofstream* out, const string& mini, uint64_t number_color);
	void merge_super_buckets(const string& input_file, uint64_t number_color);
	string compaction(const string& seq1,const string& seq2,bool);
	void construct_index_fof(const string& input_file);
	void reset();


	void dump_disk(const string& output_file);
	vector<bool> get_presence_query(const string& seq);
	vector<int64_t> get_rank_query(const string& seq);
	vector<int64_t>get_hashes_query(const string& seq);
	void file_query_presence(const string& query_file);
	void file_query_hases(const string& query_file);
	void file_query_rank(const string& query_file);
	void file_query_all_test(const string& query_file,bool);


};



// iterator from disk file of u_int64_t with buffered read,   todo template
class bfile_iterator : public std::iterator<std::forward_iterator_tag, kmer>{
	public:

	bfile_iterator()
	: _is(nullptr)
	, _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (kmer *) malloc(_buffsize*sizeof(kmer));
	}

	bfile_iterator(const bfile_iterator& cr)
	{
		_buffsize = cr._buffsize;
		_pos = cr._pos;
		_is = cr._is;
		_buffer = (kmer *) malloc(_buffsize*sizeof(kmer));
		 memcpy(_buffer,cr._buffer,_buffsize*sizeof(kmer) );
		_inbuff = cr._inbuff;
		_cptread = cr._cptread;
		_elem = cr._elem;
	}

	bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (kmer *) malloc(_buffsize*sizeof(kmer));
		int reso = fseek(_is,0,SEEK_SET);
		advance();
	}

	~bfile_iterator()
	{
		if(_buffer!=NULL)
			free(_buffer);
	}


	kmer const& operator*()  {  return _elem;  }

	bfile_iterator& operator++()
	{
		advance();
		return *this;
	}

	friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
	{
		if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
		assert(lhs._is == rhs._is);
		return rhs._pos == lhs._pos;
	}

	friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
	private:
	void advance()
	{
		_pos++;

		if(_cptread >= _inbuff)
		{
			int res = fread(_buffer,sizeof(kmer),_buffsize,_is);
			_inbuff = res; _cptread = 0;

			if(res == 0)
			{
				_is = nullptr;
				_pos = 0;
				return;
			}
		}

		_elem = _buffer[_cptread];
		_cptread ++;
	}
	kmer _elem;
	FILE * _is;
	unsigned long _pos;

	kmer * _buffer; // for buffered read
	int _inbuff, _cptread;
	int _buffsize;
};


class file_binary{
	public:

	file_binary(const char* filename)
	{
		_is = fopen(filename, "rb");
		if (!_is) {
			throw std::invalid_argument("Error opening " + std::string(filename));
		}
	}

	~file_binary()
	{
		fclose(_is);
	}

	bfile_iterator begin() const
	{
		return bfile_iterator(_is);
	}

	bfile_iterator end() const {return bfile_iterator(); }

	size_t        size () const  {  return 0;  }//todo ?

	private:
	FILE * _is;
};







#endif
