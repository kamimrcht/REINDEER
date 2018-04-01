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



#define kmer uint64_t
//~ #define kmer __uint128_t
#define minimizer_type uint32_t



typedef boomphf::SingleHashFunctor<kmer>  hasher;
typedef boomphf::mphf<  kmer, hasher  > MPHF;


struct Nadine_la_cuisine_francaise{
	vector<bool> bucketSeq;
	vector<bool> valid_kmers;
	uint abundance_minimizer;
};


struct bucket_minimizer{
	uint abundance_minimizer;
	uint nuc_minimizer;
	uint current_pos;
	vector<bool> bucketSeq;

};


struct info_mphf{
	uint mphf_size;
	uint bit_to_encode;
	vector<bool> positions;
	MPHF kmer_MPHF;
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
	uint64_t largest_bucket_nuc_all;
	uint number_bucket_per_mphf;
	kmer offsetUpdateAnchor=1;
	double bit_per_kmer;
	bool light_mode;

	//~ vector<vector<bool>> bucketSeq;
	vector<bool>* Valid_kmer;
	//~ vector<uint> abundance_minimizer;

	//~ vector<MPHF> kmer_MPHF;
	//~ vector<vector<bool>> positions;
	//~ vector<uint> mphf_size;
	//~ vector<uint> bit_to_encode;
	bucket_minimizer** all_buckets;
	info_mphf** all_mphf;
	kmer_Set_Light(uint k_val,uint m1_val, uint m2_val, uint m3_val, uint coreNumber_val){
		k=k_val;
		m1=m1_val;
		m2=m2_val;
		m3=m3_val;
		light_mode=true;
		//~ light_mode=false;
		bit_saved_sub=5;
		number_kmer=0;
		number_super_kmer=0;
		largest_bucket_nuc_all=0;
		positions_total_size=0;
		largest_MPHF=0;
		bit_per_kmer=0;
		positions_to_check=1<<bit_saved_sub;
		offsetUpdateAnchor<<=2*k;
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
		//~ abundance_minimizer.resize(minimizer_number);
		all_buckets=new bucket_minimizer*[minimizer_number];
		//~ all_mphf.resize(minimizer_number/number_bucket_per_mphf);
		all_mphf=new info_mphf*[minimizer_number/number_bucket_per_mphf];
		//~ mphf_size.resize(minimizer_number/number_bucket_per_mphf,0);
		//~ kmer_MPHF.resize(minimizer_number/number_bucket_per_mphf);
		//~ positions.resize(minimizer_number/number_bucket_per_mphf);
		//~ bit_to_encode.resize(minimizer_number/number_bucket_per_mphf,1);
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
	void get_anchors(const string& query,vector<uint>& minimizerV, vector<kmer>& kmerV);
	uint multiple_query_serial(const vector<uint>& minimizerV, const vector<kmer>& kmerV);
	void file_query(const string& query_file);
	uint32_t bool_to_int(uint n_bits_to_encode,uint pos,const vector<bool>& V);









};







#endif
