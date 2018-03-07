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



typedef boomphf::SingleHashFunctor<kmer>  hasher;
typedef boomphf::mphf<  kmer, hasher  > MPHF;




class kmer_Set_Light{
public:
	uint k;
	uint m;
	uint n;
	uint number_superbuckets;
	uint minimizer_number;
	uint coreNumber;
	uint gammaFactor;
	uint bit_saved_sub;
	uint positions_to_check;
	kmer offsetUpdateAnchor=1;


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
		bit_saved_sub=1;
		positions_to_check=1<<bit_saved_sub;
		offsetUpdateAnchor<<=2*k;
		number_superbuckets=1<<n;
		minimizer_number=1<<(2*m);
		coreNumber=coreNumber_val;
		gammaFactor=100;
		bucketSeq.resize(minimizer_number);
		//~ buckets.resize(minimizer_number);
		Valid_kmer.resize(minimizer_number);
		buckets_size.resize(minimizer_number);
		abundance_minimizer.resize(minimizer_number);
		kmer_MPHF.resize(minimizer_number);
		positions.resize(minimizer_number);
	}

	bool exists(const kmer& query);
	void create_super_buckets(const string& input_file);
	void read_super_buckets(const string& input_file);
	void create_mphf();
	void updateK(kmer& min, char nuc);
	void updateRCK(kmer& min, char nuc);
	void updateM(kmer& min, char nuc);
	void updateRCM(kmer& min, char nuc);
	void fill_positions();
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
