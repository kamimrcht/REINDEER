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
#include "bbhash.h"



using namespace std;



//~ #define kmer uint64_t
#define kmer __uint128_t
#define minimizer_type uint32_t




//~ typedef SingleHashFunctor128 hasher_t;
typedef boomphf::SingleHashFunctor<kmer>  hasher_t;


//~ typedef boomphf::mphf<  hasher_t,kmer  > MPHF;
typedef boomphf::mphf<  kmer, hasher_t  > MPHF;




struct bucket_minimizer{
	//~ uint32_t abundance_minimizer;
	uint32_t nuc_minimizer;
	uint32_t current_pos;
	uint32_t start;
};


struct info_mphf{
	uint32_t mphf_size;
	uint32_t start;
	uint8_t bit_to_encode;
	MPHF* kmer_MPHF;
};



class kmer_Set_Light{
public:
	uint k;
	uint m1,m2,m3;
	uint number_superbuckets;
	uint minimizer_number;
	uint bucket_per_superBuckets;
	uint coreNumber;
	uint gammaFactor;
	uint bit_saved_sub;
	uint positions_to_check;
	uint64_t number_kmer;
	uint64_t number_super_kmer;
	uint64_t largest_MPHF;
	uint64_t positions_total_size;
	uint largest_bucket_nuc_all;
	uint64_t number_query;
	uint number_bucket_per_mphf;
	kmer offsetUpdateAnchor=1;
	kmer offsetUpdateMinimizer=1;
	double bit_per_kmer;
	bool light_mode;

	vector<bool> bucketSeq;
	vector<bool>* Valid_kmer;
	//~ vector<uint> abundance_minimizer;

	//~ vector<MPHF> kmer_MPHF;
	vector<bool> positions;
	//~ vector<uint> mphf_size;
	//~ vector<uint> bit_to_encode;
	bucket_minimizer* all_buckets;
	uint32_t* abundance_minimizer_temp;
	uint8_t* abundance_minimizer;
	info_mphf* all_mphf;
	kmer_Set_Light(uint k_val,uint m1_val, uint m2_val, uint m3_val, uint coreNumber_val, uint bit_to_save){
		k=k_val;
		m1=m1_val;
		m2=m2_val;
		m3=m3_val;
		light_mode=true;
		//~ light_mode=false;
		bit_saved_sub=bit_to_save;
		number_kmer=0;
		number_query=0;
		number_super_kmer=0;
		largest_bucket_nuc_all=0;
		positions_total_size=0;
		largest_MPHF=0;
		bit_per_kmer=0;
		positions_to_check=1<<bit_saved_sub;
		offsetUpdateAnchor<<=2*k;
		offsetUpdateMinimizer<<=2*m1;
		number_superbuckets=1<<(2*m3);
		minimizer_number=1<<(2*m1);
		number_bucket_per_mphf=1<<(2*(m1-m2));
		bucket_per_superBuckets=minimizer_number/number_superbuckets;
		//~ cout<<minimizer_number<<endl;
		coreNumber=coreNumber_val;
		gammaFactor=2;
		//~ bucketSeq.resize(minimizer_number);
		//~ Valid_kmer.resize(bucket_per_superBuckets);
		Valid_kmer=new vector<bool>[bucket_per_superBuckets];
		for(uint i(0);i<bucket_per_superBuckets;++i){
			Valid_kmer[i]={};
		}
		//~ abundance_minimizer.resize(minimizer_number);
		all_buckets=new bucket_minimizer[minimizer_number];
		//~ abundance_minimizer_temp=new uint32_t[minimizer_number];
		//~ abundance_minimizer=new uint8_t[minimizer_number];
		//~ memset(abundance_minimizer, 0, minimizer_number*sizeof(uint8_t));
		//~ memset(abundance_minimizer_temp, 0, minimizer_number*sizeof(uint32_t));
		for(uint i(0);i<minimizer_number;++i){
			all_buckets[i]={0,0,0};
		}
		//~ all_mphf.resize(minimizer_number/number_bucket_per_mphf);
		all_mphf=new info_mphf[minimizer_number/number_bucket_per_mphf];
		for(uint i(0);i<minimizer_number/number_bucket_per_mphf;++i){
			all_mphf[i].mphf_size=0;
			all_mphf[i].bit_to_encode=0;
			all_mphf[i].start=0;
		}
		//~ mphf_size.resize(minimizer_number/number_bucket_per_mphf,0);
		//~ kmer_MPHF.resize(minimizer_number/number_bucket_per_mphf);
		//~ positions.resize(minimizer_number/number_bucket_per_mphf);
		//~ bit_to_encode.resize(minimizer_number/number_bucket_per_mphf,1);
	}

	~kmer_Set_Light () {
		delete[] all_buckets;
		for(uint i(0);i<minimizer_number/number_bucket_per_mphf;++i){
			delete all_mphf[i].kmer_MPHF;
		}
		delete[] all_mphf;
		delete[] Valid_kmer;
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
	kmer update_kmer(uint32_t pos,uint32_t mini,kmer input);
	kmer get_kmer(uint32_t pos,uint32_t mini);
	void print_kmer(kmer num);
	int32_t query_get_pos_unitig(const kmer canon,uint minimizer);
	uint32_t get_anchors(const string& query,uint& minimizer, vector<kmer>& kmerV,uint pos);
	uint multiple_query_serial(const uint minimizerV, const vector<kmer>& kmerV);
	void file_query(const string& query_file,bool bi);
	uint32_t bool_to_int(uint n_bits_to_encode,uint pos,uint64_t start);
	uint multiple_query_optimized(uint minimizerV, const vector<kmer>& kmerV);
	void int_to_bool(uint n_bits_to_encode,uint32_t X, uint32_t pos,uint64_t start);
	kmer update_kmer_local(uint32_t pos,const vector<bool>& V,kmer input);
	vector<bool> get_seq(uint32_t mini,uint32_t pos,uint32_t n);
};







#endif
