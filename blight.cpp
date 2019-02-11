
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <tmmintrin.h>

#include "bbhash.h"
#include "blight.h"
#include "zstr.hpp"
#include "common.h"




using namespace std;
using namespace chrono;



static inline kmer nuc2int(char c){
	switch(c){
		/*
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		*/
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
		case 'N': return 0;
	}
	cout<<"Unknow nucleotide: "<<c<<"!"<<endl;
	exit(0);
	return 0;
}



static inline kmer nuc2intrc(char c){
	switch(c){
		/*
		case 'a': return 0;
		case 'c': return 1;
		case 'g': return 2;
		case 't': return 3;
		*/
		case 'A': return 3;
		case 'C': return 2;
		case 'G': return 1;
		case 'T': return 0;
		case 'N': return 0;
	}
	cout<<"Unknow nucleotide: "<<c<<"!"<<endl;
	exit(0);
	return 0;
}



inline uint number_miss(const string str1,const string str2){
	uint res(0);
	for(uint i(0);i<str1.size();++i){
		if(str1[i]!=str2[i]){
			res++;
		}
	}
	return res;
}



inline string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}



inline char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



inline  string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



inline string getCanonical(const string& str){
	return (min(str,revComp(str)));
}



inline kmer str2num(const string& str){
	kmer res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			case 'T':res+=3;break;
			case 'N':res+=0;break;
			default:cout<<"bug"<<"!"<<endl;
	exit(0);
		}
	}
	return res;
}


inline uint32_t revhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



inline uint32_t unrevhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ( ( x >> 16 ) ^ x ) * 0x64ea2d65;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



inline uint64_t revhash ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



inline uint64_t unrevhash ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



template<typename T>
inline T xs(const T& x) { return revhash(x); }





// It's quite complex to bitshift mmx register without an immediate (constant) count
// See: https://stackoverflow.com/questions/34478328/the-best-way-to-shift-a-m128i
inline __m128i mm_bitshift_left(__m128i x, unsigned count)
{
	assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_slli_si128(x, 8);
	if (count >= 64) //TODO: bench: Might be faster to skip this fast-path branch
		return _mm_slli_epi64(carry, count-64);  // the non-carry part is all zero, so return early
	// else
	carry = _mm_srli_epi64(carry, 64-count);

	x = _mm_slli_epi64(x, count);
	return _mm_or_si128(x, carry);
}



inline __m128i mm_bitshift_right(__m128i x, unsigned count)
{
	assume(count < 128, "count=%u >= 128", count);
	__m128i carry = _mm_srli_si128(x, 8);
	if (count >= 64)
		return _mm_srli_epi64(carry, count-64);  // the non-carry part is all zero, so return early
	// else
	carry = _mm_slli_epi64(carry, 64-count);

	x = _mm_srli_epi64(x, count);
	return _mm_or_si128(x, carry);
}


inline __uint128_t rcb(const __uint128_t& in, uint n){
	assume(n <= 64, "n=%u > 64", n);
	union kmer_u { __uint128_t k; __m128i m128i; uint64_t u64[2]; uint8_t u8[16];};
	kmer_u res = { .k = in };
	static_assert(sizeof(res) == sizeof(__uint128_t), "kmer sizeof mismatch");

	// Complement
	res.m128i = ~res.m128i;

	// Swap byte order
	kmer_u shuffidxs = { .u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0} };
	res.m128i = _mm_shuffle_epi8 (res.m128i, shuffidxs.m128i);

	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	for(uint64_t& x : res.u64) {
		x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
		x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc
	}

	// Realign to the right
	res.m128i = mm_bitshift_right(res.m128i, 128 - 2*n);
	return res.k;
}



inline uint64_t rcb(uint64_t in, uint n) {
	assume(n <= 32, "n=%u > 32", n);
	// Complement, swap byte order
	uint64_t res = __builtin_bswap64(~in);
	// Swap nuc order in bytes
	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
	const uint64_t c2 = 0x3333333333333333;
	res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

	// Realign to the right
	res >>= 64 - 2*n;

	return res;
}



inline uint32_t rcb(uint32_t in, uint n) {
	assume(n <= 16, "n=%u > 16", n);
	// Complement, swap byte order
	uint32_t res = __builtin_bswap32(~in);

	// Swap nuc order in bytes
	const uint32_t c1 = 0x0f0f0f0f;
	const uint32_t c2 = 0x33333333;
	res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
	res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

	// Realign to the right
	res >>= 32 - 2*n;

	return res;
}



inline void kmer_Set_Light::updateK(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchor;
}



inline void kmer_Set_Light::updateM(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateMinimizer;
}



static inline kmer min_k (const kmer& k1,const kmer& k2){
	if(k1<=k2){
		return k1;
	}
	return k2;
}



inline void kmer_Set_Light::updateRCK(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-2));
}



inline void kmer_Set_Light::updateRCM(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*m1-2));
}



static inline uint32_t knuth_hash (uint32_t x){
	return x*2654435761;
}



static inline size_t hash2(int i1)
{
	size_t ret = i1;
	ret *= 2654435761U;
	return ret ^ 69;
}



static inline kmer get_int_in_kmer(kmer seq,uint64_t pos,uint number_nuc){
	seq>>=2*pos;
	return  ((seq)%(1<<(2*number_nuc)));
}



extended_minimizer kmer_Set_Light::get_extended_minimizer_from_min(kmer seq, uint32_t mini, uint position_minimizer){
	extended_minimizer res;
	res.mini=mini;
	res.suffix_fragile=res.prefix_fragile=0;
	if(position_minimizer>=extension_minimizer){
		res.extended_mini=get_int_in_kmer(seq,position_minimizer-extension_minimizer,m1);
	}else{
		res.suffix_fragile=extension_minimizer-position_minimizer;
		res.extended_mini=get_int_in_kmer(seq,0,m1-extension_minimizer+position_minimizer);
	}
	if(position_minimizer-extension_minimizer+m1>k){
		res.prefix_fragile=m1+position_minimizer-extension_minimizer-k;
	}
	res.fragile=res.suffix_fragile+res.prefix_fragile;
	if(not res.fragile){
		res.extended_mini=min(res.extended_mini,rcb(res.extended_mini,m1));
	}
	return res;
}



void kmer_Set_Light::print_extended(extended_minimizer min){
	print_kmer(min.mini);
	print_kmer(min.extended_mini);
}



extended_minimizer kmer_Set_Light::minimizer_and_more(kmer seq){
	extended_minimizer res;
	uint horrible_counter(0);
	kmer seq2(seq);
	uint32_t mini,mmer;
	mmer=seq%minimizer_number_graph;
	mini=min(mmer,rcb(mmer,minimizer_size_graph));
	res=get_extended_minimizer_from_min(seq2,mini,0);
	for(uint i(1);i<=k-minimizer_size_graph;i++){
		seq>>=2;
		mmer=seq%minimizer_number_graph;
		mmer=min(mmer,rcb(mmer,minimizer_size_graph));
		if((xs(mini)>xs(mmer))){
			horrible_counter=0;
			mini=mmer;
			res=get_extended_minimizer_from_min(seq2,mini,i);
		}else if((xs(mini)==xs(mmer))){
			extended_minimizer res2(get_extended_minimizer_from_min(seq2,mini,i));
			if(not res2.fragile){
				if(res.fragile){
					res=res2;
				}else if(res.extended_mini>res2.extended_mini){
					res=res2;
				}
			}else{
				if(res.fragile){
					horrible_counter++;
					if(horrible_counter>extension_minimizer-1){
						res.mini=0;
						res.extended_mini=0;
						res.fragile=res.prefix_fragile=res.suffix_fragile=0;

						return res;
					}
				}
			}
		}
	}
	return res;
}



uint32_t kmer_Set_Light::regular_minimizer(kmer seq){
	uint32_t mini,mmer;
	mmer=seq%minimizer_number_graph;
	mini=min(mmer,rcb(mmer,minimizer_size_graph));
	uint64_t hash_mini = xs(mini);
	for(uint i(1);i<=k-minimizer_size_graph;i++){
		seq>>=2;
		mmer=seq%minimizer_number_graph;
		mmer=min(mmer,rcb(mmer,minimizer_size_graph));
		uint64_t hash = xs(mmer);
		if(hash_mini>hash){
			mini=mmer;
			hash_mini=hash;
		}
	}
	return mini;
}



uint32_t kmer_Set_Light::minimizer_extended(kmer seq){
	kmer seq2(seq);
	uint32_t mini,mmer,position_minimizer(0);
	mini=seq%minimizer_number;
	mini=min(mini,rcb(mini,m1));
	for(uint i(0);i<k-m1;i++){
		seq>>=2;
		mmer=seq%minimizer_number;
		mmer=min(mmer,rcb(mmer,m1));
		if((xs(mini)>xs(mmer))?mini:mmer){
			mini=mmer;
			position_minimizer=i;
		}
	}
	if(position_minimizer>=extension_minimizer){
		mini=get_int_in_kmer(seq2,position_minimizer-extension_minimizer,m1+2*extension_minimizer);
	}else{
		mini=get_int_in_kmer(seq2,0,m1+position_minimizer+extension_minimizer);
		mini<<=(2*(extension_minimizer-position_minimizer));
	}
	mini=min(mini,rcb(mini,m1));
	return mini;
}




void kmer_Set_Light::abundance_minimizer_construct(const string& input_file){
	auto inUnitigs=new zstr::ifstream(input_file);
	if( not inUnitigs->good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	string ref,useless;
	while(not inUnitigs->eof()){
		getline(*inUnitigs,useless);
		getline(*inUnitigs,ref);
		//FOREACH UNITIG
		if(not ref.empty() and not useless.empty()){
			//FOREACH KMER
			kmer seq(str2num(ref.substr(0,m1))),rcSeq(rcb(seq,m1)),canon(min_k(seq,rcSeq));
				abundance_minimizer_temp[canon]++;
			uint i(0);
			for(;i+m1<ref.size();++i){
				updateM(seq,ref[i+m1]);
				updateRCM(rcSeq,ref[i+m1]);
				canon=(min_k(seq,rcSeq));
					abundance_minimizer_temp[canon]++;
			}
		}
	}
	for(uint i(0);i<minimizer_number;++i){
		abundance_minimizer[i]=(uint8_t)(log2(abundance_minimizer_temp[i])*8);
	}
	delete[] abundance_minimizer_temp;
	delete inUnitigs;
}



static inline int64_t round_eight(int64_t n){
	return n+8;
}



void kmer_Set_Light::create_super_buckets_extended(const string& input_file){
	uint64_t total_nuc_number(0);
	auto inUnitigs=new zstr::ifstream(input_file);
	if( not inUnitigs->good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	vector<ostream*> out_files;
	for(uint i(0);i<number_superbuckets;++i){
		//~ auto out =new ofstream("_out"+to_string(i));
		auto out =new zstr::ofstream("_out"+to_string(i));
		out_files.push_back(out);
	}
	omp_lock_t lock[number_superbuckets.value()];
	for (uint i=0; i<number_superbuckets; i++){
		omp_init_lock(&(lock[i]));
	}
	#pragma omp parallel num_threads(coreNumber)
	{
		string ref,useless;
		minimizer_type old_minimizer,minimizer,precise_minimizer,old_precise_minimizer;
		while(not inUnitigs->eof()){
			#pragma omp critical(dataupdate)
			{
				getline(*inUnitigs,useless);
				getline(*inUnitigs,ref);
			}
			//FOREACH UNITIG
			if(not ref.empty() and not useless.empty()){
				old_minimizer=minimizer=minimizer_number_graph.value();
				uint last_position(0);
				//FOREACH KMER
				kmer seq(str2num(ref.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				auto nadine(minimizer_and_more(canon));
				uint fragile=nadine.fragile;
				minimizer=nadine.mini;
				precise_minimizer=nadine.extended_mini;
				old_minimizer=minimizer;
				old_precise_minimizer=precise_minimizer;
				uint i(0);
				for(;i+k<ref.size();++i){
					updateK(seq,ref[i+k]);
					updateRCK(rcSeq,ref[i+k]);
					canon=(min_k(seq, rcSeq));
					//COMPUTE KMER MINIMIZER
					nadine = minimizer_and_more(canon);
					uint new_fragile(nadine.fragile);
					minimizer=nadine.mini;
					precise_minimizer=nadine.extended_mini;
					if(old_minimizer!=minimizer or ( new_fragile==fragile and old_precise_minimizer!=precise_minimizer)){
						omp_set_lock(&(lock[((old_precise_minimizer))/bucket_per_superBuckets]));
						*(out_files[((old_precise_minimizer))/bucket_per_superBuckets])<<">"+to_string(old_precise_minimizer)+"\n"<<ref.substr(last_position,i-last_position+k)<<"\n";
						omp_unset_lock(&(lock[((old_precise_minimizer))/bucket_per_superBuckets]));
						#pragma omp atomic
						all_buckets[old_precise_minimizer].nuc_minimizer+=(i-last_position+k);
						#pragma omp atomic
						all_mphf[old_precise_minimizer/number_bucket_per_mphf].mphf_size+=(i-last_position+k)-k+1;
						#pragma omp atomic
						total_nuc_number+=(i-last_position+k);
						last_position=i+1;
						old_minimizer=minimizer;
						old_precise_minimizer=precise_minimizer;
						fragile=new_fragile;
					}else{
						if(fragile > new_fragile){
							old_precise_minimizer=precise_minimizer;
							fragile=new_fragile;
						}
					}
				}
				if(ref.size()-last_position>k-1){
					omp_set_lock(&(lock[((old_precise_minimizer))/bucket_per_superBuckets]));
					*(out_files[((old_precise_minimizer))/bucket_per_superBuckets])<<">"+to_string(old_precise_minimizer)+"\n"<<ref.substr(last_position)<<"\n";
					omp_unset_lock(&(lock[((old_precise_minimizer))/bucket_per_superBuckets]));
					#pragma omp atomic
					all_buckets[old_precise_minimizer].nuc_minimizer+=(ref.substr(last_position)).size();
					#pragma omp atomic
					total_nuc_number+=(ref.substr(last_position)).size();
					#pragma omp atomic
					all_mphf[old_precise_minimizer/number_bucket_per_mphf].mphf_size+=(ref.substr(last_position)).size()-k+1;
				}
			}
		}
	}
	delete inUnitigs;
	for(uint i(0);i<number_superbuckets;++i){
		*out_files[i]<<flush;
		delete(out_files[i]);
	}
	bucketSeq.resize(total_nuc_number*2);
	bucketSeq.shrink_to_fit();
	uint64_t i(0),total_pos_size(0);
	uint max_bucket_mphf(0);
	for(uint BC(0);BC<minimizer_number;++BC){
		all_buckets[BC].start=i;
		all_buckets[BC].current_pos=i;
		i+=all_buckets[BC].nuc_minimizer;
		max_bucket_mphf=max(all_buckets[BC].nuc_minimizer,max_bucket_mphf);
		if((BC+1)%number_bucket_per_mphf==0){
			int n_bits_to_encode((ceil(log2(max_bucket_mphf+1))-bit_saved_sub));
			if(n_bits_to_encode<1){n_bits_to_encode=1;}
			all_mphf[BC/number_bucket_per_mphf].bit_to_encode=n_bits_to_encode;
			all_mphf[BC/number_bucket_per_mphf].start=total_pos_size;
			total_pos_size+=round_eight(n_bits_to_encode*all_mphf[BC/number_bucket_per_mphf].mphf_size);
			all_mphf[BC/number_bucket_per_mphf].empty=false;
			if(BC>0){
				all_mphf[BC/number_bucket_per_mphf].mphf_size+=all_mphf[(BC/number_bucket_per_mphf)-1].mphf_size;
			}
			max_bucket_mphf=0;
		}
	}
	positions.resize(total_pos_size);
	positions.shrink_to_fit();
}



void kmer_Set_Light::construct_index(const string& input_file){
	if(m1<m2){
		cout<<"n should be inferior to m"<<endl;
		exit(0);
	}
	if(m2<m3){
		cout<<"s should be inferior to n"<<endl;
		exit(0);
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	if(extension_minimizer==0){
		create_super_buckets_regular(input_file);
	}else{
		create_super_buckets_extended(input_file);
	}

	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Super bucket created: "<< time_span12.count() << " seconds."<<endl;

	read_super_buckets("_out");

	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);
	cout<<"Indexes created: "<< time_span13.count() << " seconds."<<endl;
	duration<double> time_spant = duration_cast<duration<double>>(t13 - t1);
	cout << "The whole indexing took me " << time_spant.count() << " seconds."<< endl;
}



void kmer_Set_Light::create_super_buckets_regular(const string& input_file){
	uint64_t total_nuc_number(0);
	auto inUnitigs=new zstr::ifstream(input_file);
	if( not inUnitigs->good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	vector<ostream*> out_files;
	for(uint i(0);i<number_superbuckets;++i){
		auto out =new zstr::ofstream("_out"+to_string(i));
		out_files.push_back(out);
	}
	omp_lock_t lock[number_superbuckets.value()];
	for (uint i=0; i<number_superbuckets; i++){
		omp_init_lock(&(lock[i]));
	}
	#pragma omp parallel num_threads(coreNumber)
	{
		string ref,useless;
		minimizer_type old_minimizer,minimizer;
		while(not inUnitigs->eof()){
			#pragma omp critical(dataupdate)
			{
				getline(*inUnitigs,useless);
				getline(*inUnitigs,ref);
			}
			//FOREACH UNITIG
			if(not ref.empty() and not useless.empty()){
				old_minimizer=minimizer=minimizer_number_graph.value();
				uint last_position(0);
				//FOREACH KMER
				kmer seq(str2num(ref.substr(0,k)));
				minimizer=regular_minimizer(seq);
				old_minimizer=minimizer;
				uint i(0);
				for(;i+k<ref.size();++i){
					updateK(seq,ref[i+k]);
					//COMPUTE KMER MINIMIZER
					minimizer=regular_minimizer(seq);
					if(old_minimizer!=minimizer){
						omp_set_lock(&(lock[((old_minimizer))/bucket_per_superBuckets]));
						*(out_files[((old_minimizer))/bucket_per_superBuckets])<<">"+to_string(old_minimizer)+"\n"<<ref.substr(last_position,i-last_position+k)<<"\n";
						omp_unset_lock(&(lock[((old_minimizer))/bucket_per_superBuckets]));
						#pragma omp atomic
						all_buckets[old_minimizer].nuc_minimizer+=(i-last_position+k);
						#pragma omp atomic
						all_mphf[old_minimizer/number_bucket_per_mphf].mphf_size+=(i-last_position+k)-k+1;
						all_mphf[old_minimizer/number_bucket_per_mphf].empty=false;
						#pragma omp atomic
						total_nuc_number+=(i-last_position+k);
						last_position=i+1;
						old_minimizer=minimizer;
					}
				}
				if(ref.size()-last_position>k-1){
					omp_set_lock(&(lock[((old_minimizer))/bucket_per_superBuckets]));
					*(out_files[((old_minimizer))/bucket_per_superBuckets])<<">"+to_string(old_minimizer)+"\n"<<ref.substr(last_position)<<"\n";
					omp_unset_lock(&(lock[((old_minimizer))/bucket_per_superBuckets]));
					#pragma omp atomic
					all_buckets[old_minimizer].nuc_minimizer+=(ref.substr(last_position)).size();
					#pragma omp atomic
					total_nuc_number+=(ref.substr(last_position)).size();
					#pragma omp atomic
					all_mphf[old_minimizer/number_bucket_per_mphf].mphf_size+=(ref.substr(last_position)).size()-k+1;
					all_mphf[old_minimizer/number_bucket_per_mphf].empty=false;
				}
			}
		}
	}
	delete inUnitigs;
	for(uint i(0);i<number_superbuckets;++i){
		*out_files[i]<<flush;
		delete(out_files[i]);
	}
	bucketSeq.resize(total_nuc_number*2);
	bucketSeq.shrink_to_fit();
	uint64_t i(0),total_pos_size(0);
	uint max_bucket_mphf(0);
	uint64_t hash_base(0),old_hash_base(0);
	for(uint BC(0);BC<minimizer_number;++BC){
		all_buckets[BC].start=i;
		all_buckets[BC].current_pos=i;
		i+=all_buckets[BC].nuc_minimizer;
		max_bucket_mphf=max(all_buckets[BC].nuc_minimizer,max_bucket_mphf);
		if((BC+1)%number_bucket_per_mphf==0){
			int n_bits_to_encode((ceil(log2(max_bucket_mphf+1))-bit_saved_sub));
			if(n_bits_to_encode<1){n_bits_to_encode=1;}
			all_mphf[BC/number_bucket_per_mphf].bit_to_encode=n_bits_to_encode;
			all_mphf[BC/number_bucket_per_mphf].start=total_pos_size;
			total_pos_size+=round_eight(n_bits_to_encode*all_mphf[BC/number_bucket_per_mphf].mphf_size);
			hash_base+=all_mphf[(BC/number_bucket_per_mphf)].mphf_size;
			all_mphf[BC/number_bucket_per_mphf].mphf_size=old_hash_base;
			old_hash_base=hash_base;
			max_bucket_mphf=0;
		}
	}
	positions.resize(total_pos_size);
	positions.shrink_to_fit();
}



void kmer_Set_Light::str2bool(const string& str,uint mini){
	for(uint i(0);i<str.size();++i){
		Valid_kmer[mini%bucket_per_superBuckets].push_back(true);
		switch (str[i]){
			case 'A':
				bucketSeq[(all_buckets[mini].current_pos+i)*2]=(false);
				bucketSeq[(all_buckets[mini].current_pos+i)*2+1]=(false);
				break;
			case 'C':
				bucketSeq[(all_buckets[mini].current_pos+i)*2]=(false);
				bucketSeq[(all_buckets[mini].current_pos+i)*2+1]=(true);
				break;
			case 'G':
				bucketSeq[(all_buckets[mini].current_pos+i)*2]=(true);
				bucketSeq[(all_buckets[mini].current_pos+i)*2+1]=(false);
				break;
			case 'T':
				bucketSeq[(all_buckets[mini].current_pos+i)*2]=(true);
				bucketSeq[(all_buckets[mini].current_pos+i)*2+1]=(true);
				break;
			default:
				cout<<"nope"<<endl;
			}
	}
	all_buckets[mini].current_pos+=(str.size());
	for(uint i(0);i<k-1;++i){
		Valid_kmer[mini%bucket_per_superBuckets][Valid_kmer[mini%bucket_per_superBuckets].size()-k+i+1]=(false);
	}
}



void kmer_Set_Light::read_super_buckets(const string& input_file){
	uint64_t total_size(0);
	//#pragma omp parallel num_threads(1)
	Valid_kmer=new vector<bool>[bucket_per_superBuckets.value()]();
	{
		string useless,line;
		//#pragma omp for
		for(uint SBC=0;SBC<number_superbuckets.value();++SBC){
			uint BC(SBC*bucket_per_superBuckets);
			zstr::ifstream in((input_file+to_string(SBC)));
			while(not in.eof() and in.good()){
				useless="";
				getline(in,useless);
				getline(in,line);
				if(not useless.empty()){
					useless=useless.substr(1);
					uint minimizer(stoi(useless));
					str2bool(line,minimizer);
					//#pragma omp atomic
					number_kmer+=line.size()-k+1;
					//#pragma omp atomic
					number_super_kmer++;
				}
			}
			remove((input_file+to_string(SBC)).c_str());
			create_mphf(BC,BC+bucket_per_superBuckets);
			fill_positions(BC,BC+bucket_per_superBuckets);
			BC+=bucket_per_superBuckets;
			cout<<"-"<<flush;
		}
	}
	delete[] Valid_kmer;
	cout<<endl;
	cout<<"----------------------INDEX RECAP----------------------------"<<endl;
	cout<<"Kmer in graph: "<<intToString(number_kmer)<<endl;
	cout<<"Super Kmer in graph: "<<intToString(number_super_kmer)<<endl;
	cout<<"Average size of Super Kmer: "<<intToString(number_kmer/(number_super_kmer))<<endl;
	cout<<"Total size of the partitionned graph: "<<intToString(bucketSeq.size()/2)<<endl;
	cout<<"Total size of the partitionned graph: "<<intToString(bucketSeq.capacity()/2)<<endl;
	cout<<"Largest MPHF: "<<intToString(largest_MPHF)<<endl;
	cout<<"Largest Bucket: "<<intToString(largest_bucket_nuc_all)<<endl;

	cout<<"Size of the partitionned graph (MBytes): "<<intToString(bucketSeq.size()/(8*1024*1024))<<endl;
	if(not light_mode){
		cout<<"Space used for separators (MBytes): "<<intToString(total_size/(8*1024*1024))<<endl;
	}
	cout<<"Total Positions size (MBytes): "<<intToString(positions.size()/(8*1024*1024))<<endl;
	cout<<"Size of the partitionned graph (bit per kmer): "<<((double)(bucketSeq.size())/(number_kmer))<<endl;
	bit_per_kmer+=((double)(bucketSeq.size())/(number_kmer));
	cout<<"Total Positions size (bit per kmer): "<<((double)positions.size()/number_kmer)<<endl;
	bit_per_kmer+=((double)positions.size()/number_kmer);
	cout<<"TOTAL Bits per kmer (without bbhash): "<<bit_per_kmer<<endl;
	cout<<"TOTAL Bits per kmer (with bbhash): "<<bit_per_kmer+4<<endl;
	cout<<"TOTAL Size estimated (MBytes): "<<(bit_per_kmer+4)*number_kmer/(8*1024*1024)<<endl;

	boomphf::memreport_t report;
	report_memusage(report);
	boomphf::print_memreport(report);
}



inline kmer kmer_Set_Light::get_kmer(uint64_t mini,uint64_t pos){
	kmer res(0);
	uint64_t bit = (all_buckets[mini].start+pos)*2;
	const uint64_t bitlast = bit + 2*k;
	for(;bit<bitlast;bit+=2){
		res<<=2;
		res |= bucketSeq[bit]*2 | bucketSeq[bit+1];
	}
	return res;
}



vector<bool> kmer_Set_Light::get_seq(uint32_t mini,uint64_t pos,uint32_t n){
	return vector<bool>(bucketSeq.begin()+(all_buckets[mini].start+pos)*2,bucketSeq.begin()+(all_buckets[mini].start+pos+n)*2);
}



inline kmer kmer_Set_Light::update_kmer(uint64_t pos,uint32_t mini,kmer input){
	return update_kmer_local(all_buckets[mini].start+pos, bucketSeq, input);
}



inline kmer kmer_Set_Light::update_kmer_local(uint64_t pos,const vector<bool>& V,kmer input){
	input<<=2;
	uint64_t bit0 = pos*2;
	input |= V[bit0]*2 | V[bit0+1];
	return input%offsetUpdateAnchor;
}



void kmer_Set_Light::print_kmer(kmer num,uint n){
	Pow2<kmer> anc(2*(k-1));
	for(uint i(0);i<k and i<n;++i){
		uint nuc=num/anc;
		num=num%anc;
		if(nuc==3){
			cout<<"T";
		}
		if(nuc==2){
			cout<<"G";
		}
		if(nuc==1){
			cout<<"C";
		}
		if(nuc==0){
			cout<<"A";
		}
		if (nuc>=4){
			cout<<nuc<<endl;
			cout<<"WTF"<<endl;
		}
		anc>>=2;
	}
	cout<<endl;
}



inline string kmer_Set_Light::kmer2str(kmer num){
	string res;
	Pow2<kmer> anc(2*(k-1));
	for(uint i(0);i<k;++i){
		uint nuc=num/anc;
		num=num%anc;
		if(nuc==3){
			res+="T";
		}
		if(nuc==2){
			res+="G";
		}
		if(nuc==1){
			res+="C";
		}
		if(nuc==0){
			res+="A";
		}
		if (nuc>=4){
			cout<<nuc<<endl;
			cout<<"WTF"<<endl;
		}
		anc>>=2;
	}
	return res;
}






void kmer_Set_Light::create_mphf(uint begin_BC,uint end_BC){
	#pragma omp parallel  num_threads(coreNumber)
		{
		vector<kmer> anchors;
		uint largest_bucket_anchor(0);
		uint largest_bucket_nuc(0);
		#pragma omp for schedule(dynamic, number_bucket_per_mphf.value())
		for(uint BC=(begin_BC);BC<end_BC;++BC){
			if(all_buckets[BC].nuc_minimizer!=0){
				largest_bucket_nuc=max(largest_bucket_nuc,all_buckets[BC].nuc_minimizer);
				largest_bucket_nuc_all=max(largest_bucket_nuc_all,all_buckets[BC].nuc_minimizer);
				uint bucketSize(1);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				anchors.push_back(canon);
				for(uint j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
					if(not Valid_kmer[BC%bucket_per_superBuckets][j+1]){
					//~ if(false){
						j+=k-1;
						if((j+k)<all_buckets[BC].nuc_minimizer){
							seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
							anchors.push_back(canon);
							bucketSize++;
						}
					}else{
						seq=update_kmer(j+k,BC,seq);
						rcSeq=(rcb(seq,k));
						canon=(min_k(seq, rcSeq));
						anchors.push_back(canon);
						bucketSize++;
					}
				}
				largest_bucket_anchor=max(largest_bucket_anchor,bucketSize);
			}
			if((BC+1)%number_bucket_per_mphf==0 and not anchors.empty()){
				largest_MPHF=max(largest_MPHF,anchors.size());
				all_mphf[BC/number_bucket_per_mphf].kmer_MPHF= new boomphf::mphf<kmer,hasher_t>(anchors.size(),anchors,gammaFactor);
				anchors.clear();
				largest_bucket_anchor=0;
				largest_bucket_nuc=(0);
			}
		}
	}
}



void kmer_Set_Light::int_to_bool(uint n_bits_to_encode,uint64_t X, uint64_t pos,uint64_t start){
	for(uint64_t i(0);i<n_bits_to_encode;++i){
		positions[i+pos*n_bits_to_encode+start]=X%2;
		X>>=1;
	}
}



uint32_t kmer_Set_Light::bool_to_int(uint n_bits_to_encode,uint64_t pos,uint64_t start){
	uint32_t res(0);
	uint32_t acc(1);
	for(uint64_t i(0);i<n_bits_to_encode;++i, acc<<=1){
		if(positions[i+pos*n_bits_to_encode+start]){
			res |= acc;
		}
	}
	return res*positions_to_check;
}



void kmer_Set_Light::fill_positions(uint begin_BC,uint end_BC){
	#pragma omp parallel for num_threads(coreNumber)
	for(uint BC=(begin_BC);BC<end_BC;++BC){
			if(all_buckets[BC].nuc_minimizer>0){
				int n_bits_to_encode(all_mphf[BC/number_bucket_per_mphf].bit_to_encode);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				for(uint j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
					if(not Valid_kmer[BC%bucket_per_superBuckets][j+1]){
						j+=k-1;
						if((j+k)<all_buckets[BC].nuc_minimizer){
							seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
							//~ #pragma omp critical(dataupdate)
							{
								int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
							}
						}
					}else{
						seq=update_kmer(j+k,BC,seq);
						rcSeq=(rcb(seq,k));
						canon=(min_k(seq, rcSeq));
						//~ #pragma omp critical(dataupdate)
						{
							int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
						}
					}
				}
			}
	}
	for(uint BC=(begin_BC);BC<end_BC;++BC){
		Valid_kmer[BC%bucket_per_superBuckets].clear();
	}
}



int64_t kmer_Set_Light::correct_pos(uint32_t mini, uint64_t p){
	if(Valid_kmer[mini%bucket_per_superBuckets].size()<p+k){
		return p;
	}
	for(uint i(0);i<k;i++){
		if(Valid_kmer[mini%bucket_per_superBuckets][p+i]){
			return (p+i);
		}
	}
	return p;
}



bool kmer_Set_Light::query_kmer_bool(kmer canon){
	if(extension_minimizer>0){
		auto nadine(minimizer_and_more(canon));
		uint fragile(nadine.fragile);
		uint32_t minimizer=nadine.extended_mini;
		if(fragile){
			if(multiple_minimizer_query_bool(minimizer,  canon, nadine.prefix_fragile,nadine.suffix_fragile)){
				return true;
			}else{
				return false;
			}
		}else{
			if(single_query(minimizer,canon)){
				return true;
			}else{
				return false;
			}
		}
	}else{
		uint32_t min(regular_minimizer(canon));
		return single_query(min,canon);
	}
}


int64_t kmer_Set_Light::query_kmer_hash(kmer canon){
	if(extension_minimizer>0){
		auto nadine(minimizer_and_more(canon));
		uint fragile(nadine.fragile);
		uint32_t minimizer=nadine.extended_mini;
		if(fragile){
			return multiple_minimizer_query_hash(minimizer,  canon, nadine.prefix_fragile,nadine.suffix_fragile);
		}else{
			return query_get_hash(canon,minimizer);
		}
	}else{
		return query_get_hash(canon,regular_minimizer(canon));
	}
}


pair<uint32_t,uint32_t> kmer_Set_Light::query_sequence_bool(const string& query){
	uint res(0);
	uint fail(0);
	if(query.size()<k){
		return make_pair(0,0);
	}
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
	uint i(0);
	canon=(min_k(seq, rcSeq));
	if(query_kmer_bool(canon)){++res;}else{++fail;}
	for(;i+k<query.size();++i){
		updateK(seq,query[i+k]);
		updateRCK(rcSeq,query[i+k]);
		canon=(min_k(seq, rcSeq));
		if(query_kmer_bool(canon)){++res;}else{++fail;}
	}
	return make_pair(res,fail);
}



vector<int64_t> kmer_Set_Light::query_sequence_hash(const string& query){
	vector<int64_t> res;
	if(query.size()<k){
		return res;
	}
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
	uint i(0);
	canon=(min_k(seq, rcSeq));
	res.push_back(query_kmer_hash(canon));
	for(;i+k<query.size();++i){
		updateK(seq,query[i+k]);
		updateRCK(rcSeq,query[i+k]);
		canon=(min_k(seq, rcSeq));
		res.push_back(query_kmer_hash(canon));
	}
	return res;
}




uint kmer_Set_Light::multiple_query_serial(const uint minimizer, const vector<kmer>& kmerV){
	uint res(0);
	for(uint i(0);i<kmerV.size();++i){
		if(single_query(minimizer,kmerV[i])){
			++res;
		}
	}
	return res;
}



bool kmer_Set_Light::multiple_minimizer_query_bool(uint minimizer, kmer kastor,uint prefix_length,uint suffix_length){
	if(suffix_length>0){
		const Pow2<uint32_t> max_completion(2*suffix_length);
		minimizer*=max_completion;
		for(uint i(0);i<max_completion;++i){
			uint32_t poential_min(min(minimizer|i,rcb(minimizer|i,m1)));
			if(single_query(poential_min,kastor)){
				return true;
			}
		}
	}
	if(prefix_length>0){
		const Pow2<uint32_t> max_completion(2*(prefix_length));
		const Pow2<uint32_t> mask(2*(m1-prefix_length));
		for(uint i(0);i<max_completion;++i){
			uint32_t poential_min(min(minimizer|i*mask,rcb(minimizer|i*mask,m1)));
			if(single_query(poential_min,kastor)){
				return true;
			}
		}
	}
	return false;
}


int64_t kmer_Set_Light::multiple_minimizer_query_hash(uint minimizer, kmer kastor,uint prefix_length,uint suffix_length){
	if(suffix_length>0){
		uint32_t max_completion(1);
		max_completion<<=(2*suffix_length);
		minimizer<<=(2*suffix_length);
		for(uint i(0);i<max_completion;++i){
			uint32_t poential_min(min(minimizer+i,rcb(minimizer+i,m1)));
			return query_get_hash(kastor,poential_min);
		}
	}
	if(prefix_length>0){
		uint32_t max_completion(1);
		uint32_t mask(1);
		max_completion<<=(2*(prefix_length));
		mask<<=(2*(m1-prefix_length));
		for(uint i(0);i<max_completion;++i){
			uint32_t poential_min(min(minimizer+i*mask,rcb(minimizer+i*mask,m1)));
			return query_get_hash(kastor,poential_min);
		}
	}
	return -1;
}



bool kmer_Set_Light::single_query(const uint minimizer, kmer kastor){
	return (query_get_pos_unitig(kastor,minimizer)>=0);
}



static inline uint next_different_value(const vector<uint>& minimizerV,uint start, uint m){
	uint i(0);
	for(;i+start<minimizerV.size();++i){
		if(minimizerV[i+start]!=m){
			return start+i-1;
		}
	}
	return start+i-1;
}



uint kmer_Set_Light::multiple_query_optimized(uint32_t minimizer, const vector<kmer>& kmerV){
	uint res(0);
	for(uint i(0);i<kmerV.size();++i){
		uint64_t pos=query_get_pos_unitig(kmerV[i],minimizer);
		uint next(kmerV.size()-1);
		if(next!=i){
			uint64_t pos2=query_get_pos_unitig(kmerV[next],minimizer);
			if(pos2-pos==next-i){
				res+=next-i+1;
				i=next;
			}else{
				++res;
			}
		}else{
			++res;
		}
	}
	return res;
}



int32_t kmer_Set_Light::query_get_pos_unitig(const kmer canon,uint minimizer){
	#pragma omp atomic
	number_query++;
	if(unlikely(all_mphf[minimizer/number_bucket_per_mphf].empty))
		return -1;

	uint64_t hash=(all_mphf[minimizer/number_bucket_per_mphf].kmer_MPHF->lookup(canon));
	if(unlikely(hash == ULLONG_MAX))
		return -1;

	int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf].bit_to_encode);
	uint64_t pos(bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf].start));
	if(likely((pos+k-1)<all_buckets[minimizer].nuc_minimizer)){
		kmer seqR=get_kmer(minimizer,pos);
		kmer rcSeqR, canonR;
		for(uint64_t j=(pos);j<pos+positions_to_check;++j){
			rcSeqR=(rcb(seqR,k));
			canonR=(min_k(seqR, rcSeqR));
			if(canon==canonR){
				return j;
			}
			seqR=update_kmer(j+k,minimizer,seqR);//can be avoided
		}
	}
	return -1;
}


int64_t kmer_Set_Light::query_get_hash(const kmer canon,uint minimizer){
	#pragma omp atomic
	number_query++;
	if(unlikely(all_mphf[minimizer/number_bucket_per_mphf].empty))
		return -1;

	uint64_t hash=(all_mphf[minimizer/number_bucket_per_mphf].kmer_MPHF->lookup(canon));
	if(unlikely(hash == ULLONG_MAX))
		return -1;

	int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf].bit_to_encode);
	uint64_t pos(bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf].start));
	if(likely((pos+k-1)<all_buckets[minimizer].nuc_minimizer)){
		kmer seqR=get_kmer(minimizer,pos);
		kmer rcSeqR, canonR;
		for(uint64_t j=(pos);j<pos+positions_to_check;++j){
			rcSeqR=(rcb(seqR,k));
			canonR=(min_k(seqR, rcSeqR));
			if(canon==canonR){
				return hash+all_mphf[minimizer/number_bucket_per_mphf].mphf_size;
			}
			seqR=update_kmer(j+k,minimizer,seqR);//can be avoided
		}
	}
	return -1;
}



void kmer_Set_Light::file_query(const string& query_file){
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	auto in=new zstr::ifstream(query_file);
	uint64_t TP(0),FP(0);
	#pragma omp parallel num_threads(coreNumber)
	{
		vector<kmer> kmerV;
		while(not in->eof() and in->good()){
			string query;
			#pragma omp critical(dataupdate)
			{
				getline(*in,query);
				getline(*in,query);
			}
			if(query.size()>=k){
				pair<uint,uint> pair(query_sequence_bool(query));
				#pragma atomic
				TP+=pair.first;
				#pragma atomic
				FP+=pair.second;
			}
		}
	}
	cout<<"-----------------------QUERY RECAP 2----------------------------"<<endl;
	cout<<"Good kmer: "<<intToString(TP)<<endl;
	cout<<"Erroneous kmers: "<<intToString(FP)<<endl;
	cout<<"Query performed: "<<intToString(number_query)<<endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole QUERY took me " << time_span.count() << " seconds."<< endl;
	delete in;
}




void kmer_Set_Light::report_memusage(boomphf::memreport_t& report, const std::string& prefix, bool add_struct) {
	if(add_struct)
		report[prefix+"::sizeof(struct)"] += sizeof(kmer_Set_Light);
	report[prefix+"::positions"] += positions.size() / CHAR_BIT;
	report[prefix+"::bucketSeq"] += bucketSeq.size() / CHAR_BIT;

	report[prefix+"::sizeof(bucket_minimizer)*minimizer_number"] += sizeof(bucket_minimizer) * minimizer_number;
	report[prefix+"::sizeof(info_mphf)*mphf_number"] += sizeof(info_mphf) * mphf_number;
	for(uint i(0);i<mphf_number;++i){
		if(all_mphf[i].kmer_MPHF)
			all_mphf[i].kmer_MPHF->report_memusage(report, prefix+"::kmer_MPHF");
	}
}





