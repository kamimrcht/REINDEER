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




class kmer_Set_Light{
public:
	uint k;
	uint m;
	uint n;
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
	kmer offsetUpdateAnchor=1;
	double bit_per_kmer;
	bool light_mode;

	//~ vector<string> 	buckets;
	vector<vector<bool>> bucketSeq;
	vector<vector<bool>> Valid_kmer;
	vector<MPHF> kmer_MPHF;
	vector<vector<bool>> positions;
	vector<uint> buckets_size;
	vector<uint> abundance_minimizer;
	kmer_Set_Light(uint k_val,uint m_val, uint n_val,uint coreNumber_val){
		k=k_val;
		m=m_val;
		n=n_val;
		light_mode=true;
		//~ light_mode=false;
		bit_saved_sub=8;
		number_kmer=0;
		number_super_kmer=0;
		positions_total_size=0;
		largest_MPHF=0;
		bit_per_kmer=0;
		positions_to_check=1<<bit_saved_sub;
		offsetUpdateAnchor<<=2*k;
		number_superbuckets=1<<(2*n);
		minimizer_number=1<<(2*m);
		bucket_per_superBuckets=minimizer_number/number_superbuckets;
		coreNumber=coreNumber_val;
		gammaFactor=5;
		bucketSeq.resize(minimizer_number);
		//~ buckets.resize(minimizer_number);
		Valid_kmer.resize(minimizer_number);
		buckets_size.resize(minimizer_number,0);
		abundance_minimizer.resize(minimizer_number);
		kmer_MPHF.resize(minimizer_number);
		positions.resize(minimizer_number);
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






};







#endif
