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
#include "common.h"
#include "bbhash.h"



using namespace std;


// FOR k<32
#define kmer uint64_t
// FOR k<64
//#define kmer __uint128_t




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
};



struct info_mphf{
	uint64_t mphf_size;
	uint64_t start;
	MPHF* kmer_MPHF;
	uint8_t bit_to_encode;
	bool empty;
};



struct extended_minimizer{
	uint prefix_fragile;
	uint suffix_fragile;
	uint fragile;
	uint32_t mini;
	uint32_t original;
	uint32_t extended_mini;
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
	const uint k, m1, m2, m3, extension_minimizer, minimizer_size_graph;
	const uint coreNumber;
	const uint bit_saved_sub;

	const Pow2<kmer> offsetUpdateAnchor;
	const Pow2<kmer> offsetUpdateMinimizer;
	const Pow2<uint> mphf_number;
	const Pow2<uint> number_superbuckets;
	const Pow2<uint> minimizer_number;
	const Pow2<uint> minimizer_number_graph;
	const Pow2<uint> number_bucket_per_mphf;
	const Pow2<uint> bucket_per_superBuckets;
	const Pow2<uint> positions_to_check;

	vector<bool> bucketSeq;
	vector<bool> positions;
	vector<bool>* Valid_kmer;
	bucket_minimizer* all_buckets;
	uint32_t* abundance_minimizer_temp;
	uint8_t* abundance_minimizer;
	info_mphf* all_mphf;

	uint64_t number_kmer = 0;
	uint64_t number_super_kmer = 0;
	uint64_t largest_MPHF = 0;
	uint64_t positions_total_size = 0;
	uint64_t number_query=0;


	//~ uint hell_bucket;
	double bit_per_kmer = 0;
	uint largest_bucket_nuc_all = 0;
	const uint gammaFactor=2;
	const bool light_mode=true;

	kmer_Set_Light(uint k_val,uint m1_val, uint m2_val, uint m3_val, uint coreNumber_val, uint bit_to_save,uint ex)
		: k(k_val)
		, m1(m1_val)
		, m2(m2_val ? m2_val : m1)
		, m3(m3_val)
		, extension_minimizer(ex)
		, minimizer_size_graph(m1-2*extension_minimizer)
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
		, positions_to_check(bit_saved_sub)
	{
		all_buckets=new bucket_minimizer[minimizer_number.value()]();
		all_mphf=new info_mphf[mphf_number.value()];
		for(uint i(0);i<mphf_number;++i){
			all_mphf[i].mphf_size=0;
			all_mphf[i].bit_to_encode=0;
			all_mphf[i].start=0;
			all_mphf[i].empty=true;
		}
	}

	~kmer_Set_Light () {
		delete[] all_buckets;
		for(uint i(0);i<mphf_number;++i){
			delete all_mphf[i].kmer_MPHF;
		}
		delete[] all_mphf;
		//~ delete[] abundance_minimizer;;
	}

	bool exists(const kmer& query);
	void create_super_buckets(const string& input_file);
	void read_super_buckets(const string& input_file);
	void create_mphf(uint32_t beg,uint32_t end);
	void updateK(kmer& min, char nuc);
	void updateRCK(kmer& min, char nuc);
	void updateM(kmer& min, char nuc);
	void updateRCM(kmer& min, char nuc);
	void fill_positions(uint32_t beg,uint32_t end);
	bool exists(const string& query);
	void multiple_query(const string& query);
	uint32_t minimizer_according_xs(kmer seq);
	void abundance_minimizer_construct(const string& input_file);
	int64_t correct_pos(uint32_t mini, uint64_t p);
	void str2bool(const string& str,uint mini);
	kmer update_kmer(uint64_t pos,uint32_t mini,kmer input);
	kmer get_kmer(uint64_t pos,uint64_t mini);
	void print_kmer(kmer num,uint n=100);
	int32_t query_get_pos_unitig(const kmer canon,uint minimizer);
	uint32_t get_anchors(const string& query,uint& minimizer, vector<kmer>& kmerV,uint pos);
	uint multiple_query_serial(const uint minimizerV, const vector<kmer>& kmerV);
	void file_query(const string& query_file);
	uint32_t bool_to_int(uint n_bits_to_encode,uint64_t pos,uint64_t start);
	uint multiple_query_optimized(uint minimizerV, const vector<kmer>& kmerV);
	void int_to_bool(uint n_bits_to_encode,uint64_t X, uint64_t pos,uint64_t start);
	kmer update_kmer_local(uint64_t pos,const vector<bool>& V,kmer input);
	vector<bool> get_seq(uint32_t mini,uint64_t pos,uint32_t n);
	uint32_t minimizer_graph(kmer seq);
	uint32_t minimizer_extended(kmer seq);
	pair<uint32_t,uint32_t> minimizer_and_more(kmer seq, uint& prefix_fragile, uint& suffix_fragile);
	bool single_query(const uint minimizer, kmer kastor);
	bool multiple_minimizer_query_bool(const uint minimizer,  kmer kastor,uint prefix_length,uint suffix_length);
	int64_t multiple_minimizer_query_hash(const uint minimizer,  kmer kastor,uint prefix_length,uint suffix_length);
	bool query_kmer_bool(kmer canon);
	pair<uint32_t,uint32_t> query_sequence_bool(const string& query);
	string kmer2str(kmer num);
	extended_minimizer minimizer_and_more(kmer seq);
	extended_minimizer get_extended_minimizer_from_min(kmer seq, uint32_t mini, uint position_minimizer);
	void print_extended(extended_minimizer);
	uint32_t regular_minimizer(kmer seq);
	void create_super_buckets_extended(const string&);
	void create_super_buckets_regular(const string&);
	int64_t query_kmer_hash(kmer canon);
	int64_t query_get_hash(const kmer canon,uint32_t minimizer);
	vector<int64_t> query_sequence_hash(const string& query);
	void construct_index(const string& input_file);
	void report_memusage(boomphf::memreport_t& report, const std::string& prefix="blight", bool add_struct=true);
};







#endif
