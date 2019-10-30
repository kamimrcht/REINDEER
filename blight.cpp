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
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include "bbhash.h"
#include "blight.h"
#include "zstr.hpp"
#include "common.h"




using namespace std;
using namespace chrono;



uint64_t asm_log2(const uint64_t x) {
  uint64_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}



static inline kmer nuc2int(char c){
	return (c/2)%4;
}



static inline kmer nuc2intrc(char c){
	return ((c/2)%4)^2;
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
		res+=(str[i]/2)%4;
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



static inline uint32_t knuth_hash (uint32_t x){
	return x*2654435761;
}



kmer xors(kmer y){
	y^=(y<<13); y^=(y>>17);y=(y^=(y<<15)); return y;
}



kmer hash64shift(kmer key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}



static inline size_t hash2(int i1)
{
	size_t ret = i1;
	ret *= 2654435761U;
	return ret ^ 69;
}



template<typename T>
inline T xs(const T& x) { return hash64shift(x); }





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

	// Swap byte order

	kmer_u shuffidxs = { .u8 = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0} };

	res.m128i = _mm_shuffle_epi8 (res.m128i, shuffidxs.m128i);

	// Swap nuc order in bytes

	const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;

	const uint64_t c2 = 0x3333333333333333;

	for(uint64_t& x : res.u64) {

		x = ((x & c1) << 4) | ((x & (c1 << 4)) >> 4); // swap 2-nuc order in bytes

		x = ((x & c2) << 2) | ((x & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

		x ^= 0xaaaaaaaaaaaaaaaa; // Complement;

	}

	// Realign to the right

	res.m128i = mm_bitshift_right(res.m128i, 128 - 2*n);

	return res.k;

}



inline uint64_t rcb(uint64_t in, uint n){
    assume(n <= 32, "n=%u > 32", n);
    // Complement, swap byte order
    uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
    // Swap nuc order in bytes
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
    const uint64_t c2 = 0x3333333333333333;
    res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
    res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

    // Realign to the right
    res >>= 64 - 2 * n;
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



static inline kmer get_int_in_kmer(kmer seq,uint64_t pos,uint number_nuc){
	seq>>=2*pos;
	return  ((seq)%(1<<(2*number_nuc)));
}



kmer canonize(kmer x, uint k){
	return min(x,rcb(x,k));
	return ( (__builtin_popcountll(x) & 1) ? x : rcb(x, k)) >> 1;
}


void print_bin(uint64_t n){
	uint64_t mask=1;
	mask<<=63;
	for(uint i(0);i<64;++i){
		cout<<n/mask;
		if(n/mask==1){n-=mask;}
		mask>>=1;
	}
	cout<<"\n";
}


kmer kmer_Set_Light::mantis(uint64_t n){
	return n;
	if(n==0){return (0);}
	int64_t prefix((asm_log2(n)));
	int64_t exp=prefix;
	int offset=prefix-(minimizer_number.bits()-6);
	offset=max(offset,0);
	uint64_t suffix (n-((uint64_t)1<<prefix));
	suffix>>=offset;
	uint64_t res=suffix;
	res+=((exp)<<((minimizer_number.bits()-6)));
	return res;
}



kmer kmer_Set_Light::regular_minimizer(kmer seq){
	kmer mini,mmer;
	mmer=seq%minimizer_number_graph;
	mini=mmer=canonize(mmer,minimizer_size_graph);
	uint64_t hash_mini = (xs(mmer));
	for(uint i(1);i<=k-minimizer_size_graph;i++){
		seq>>=2;
		mmer=seq%minimizer_number_graph;
		mmer=canonize(mmer,minimizer_size_graph);
		uint64_t hash = (xs(mmer));
		if(hash_mini>hash){
			mini=mmer;
			hash_mini=hash;
		}
	}
	return revhash((uint64_t)mini)%minimizer_number;
	//~ return mini>>(minimizer_number_graph.bits()-minimizer_number.bits());
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
	return n;
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

	create_super_buckets_regular(input_file);

	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Super bucket created: "<< time_span12.count() << " seconds."<<endl;

	read_super_buckets("_blout");

	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);
	cout<<"Indexes created: "<< time_span13.count() << " seconds."<<endl;
	duration<double> time_spant = duration_cast<duration<double>>(t13 - t1);
	cout << "The whole indexing took me " << time_spant.count() << " seconds."<< endl;
}



void kmer_Set_Light::create_super_buckets_regular(const string& input_file,bool clean){
	uint64_t total_nuc_number(0);
	auto inUnitigs=new zstr::ifstream(input_file);
	if( not inUnitigs->good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	vector<ostream*> out_files;
	for(uint i(0);i<number_superbuckets;++i){
		if(clean){
			auto out =new zstr::ofstream("_blout"+to_string(i));
			out_files.push_back(out);
		}else{
			auto out =new zstr::ofstream("_blout"+to_string(i),ofstream::app);
			out_files.push_back(out);
		}

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
				old_minimizer=minimizer=minimizer_number.value();
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
						all_buckets[old_minimizer].skmer_number++;
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
					all_buckets[old_minimizer].skmer_number++;
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
	uint64_t hash_base(0),old_hash_base(0), nb_skmer_before(0), last_skmer_number(0);
	for(uint BC(0);BC<minimizer_number.value();++BC){
		all_buckets[BC].start=i;
		all_buckets[BC].current_pos=i;
		i+=all_buckets[BC].nuc_minimizer;
		max_bucket_mphf=max(all_buckets[BC].skmer_number,max_bucket_mphf);
		if (BC == 0)
		{
			nb_skmer_before = 0; // I replace skmer_number by the total number of minitigs before this bucket
		}
		else
		{
			nb_skmer_before = all_buckets[BC-1].skmer_number;
		}
		uint64_t local_skmercount( all_buckets[BC].skmer_number);
		all_buckets[BC].skmer_number=total_nb_minitigs;
		total_nb_minitigs+=local_skmercount;
		if((BC+1)%number_bucket_per_mphf==0){
			int n_bits_to_encode((ceil(log2(max_bucket_mphf+1)))-bit_saved_sub);
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
	total_nb_minitigs = all_buckets[(uint) minimizer_number - 1].skmer_number + last_skmer_number; // total number of minitigs
	positions.resize(total_pos_size);
	positions.shrink_to_fit();
}



void kmer_Set_Light::str2bool(const string& str,uint mini){
	for(uint i(0);i<str.size();++i){
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
				bucketSeq[(all_buckets[mini].current_pos+i)*2+1]=(true);
				break;
			case 'T':
				bucketSeq[(all_buckets[mini].current_pos+i)*2]=(true);
				bucketSeq[(all_buckets[mini].current_pos+i)*2+1]=(false);
				break;
			default:
				cout<<"nope"<<endl;
			}
	}
	all_buckets[mini].current_pos+=(str.size());
}



void kmer_Set_Light::read_super_buckets(const string& input_file){
	uint64_t total_size(0);
	#pragma omp parallel num_threads(coreNumber)
	uint32_t kmer_number_in_buckets(1);
	{
		string useless,line;
		#pragma omp for
		for(uint SBC=0;SBC<number_superbuckets.value();++SBC){
			vector<uint32_t> number_kmer_accu(bucket_per_superBuckets.value(),0);
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
					position_super_kmers[number_kmer_accu[minimizer%bucket_per_superBuckets]+all_mphf[minimizer].mphf_size]=true;
					number_kmer+=line.size()-k+1;
					number_kmer_accu[minimizer%bucket_per_superBuckets]+=line.size()-k+1;
					number_super_kmer++;
				}
			}
			remove((input_file+to_string(SBC)).c_str());
			//~ create_mphf_mem(BC,BC+bucket_per_superBuckets);
			create_mphf_disk(BC,BC+bucket_per_superBuckets);
			position_super_kmers.optimize();
			position_super_kmers.optimize_gap_size();
			fill_positions(BC,BC+bucket_per_superBuckets);
			BC+=bucket_per_superBuckets;
			//~ cout<<"-"<<flush;
		}
	}
	position_super_kmers.optimize();
	position_super_kmers.optimize_gap_size();
	position_super_kmers_RS=new bm::bvector<>::rs_index_type();
	position_super_kmers.build_rs_index(position_super_kmers_RS);
	cout<<endl;
	cout<<"----------------------INDEX RECAP----------------------------"<<endl;
	cout<<"Kmer in graph: "<<intToString(number_kmer)<<endl;
	cout<<"Super Kmer in graph: "<<intToString(number_super_kmer)<<endl;
	cout<<"Average size of Super Kmer: "<<intToString(number_kmer/(number_super_kmer))<<endl;
	//~ cout<<"Total size of the partitionned graph: "<<intToString(bucketSeq.size()/2)<<endl;
	cout<<"Total size of the partitionned graph: "<<intToString(bucketSeq.capacity()/2)<<endl;
	cout<<"Largest MPHF: "<<intToString(largest_MPHF)<<endl;
	cout<<"Largest Bucket: "<<intToString(largest_bucket_nuc_all)<<endl;

	cout<<"Size of the partitionned graph (MBytes): "<<intToString(bucketSeq.size()/(8*1024*1024))<<endl;
	cout<<"Total Positions size (MBytes): "<<intToString(positions.size()/(8*1024*1024))<<endl;
	cout<<"Size of the partitionned graph (bit per kmer): "<<((double)(bucketSeq.size())/(number_kmer))<<endl;
	bit_per_kmer+=((double)(bucketSeq.size())/(number_kmer));
	cout<<"Total Positions size (bit per kmer): "<<((double)positions.size()/number_kmer)<<endl;
	bit_per_kmer+=((double)positions.size()/number_kmer);
	cout<<"TOTAL Bits per kmer (without bbhash): "<<bit_per_kmer<<endl;
	cout<<"TOTAL Bits per kmer (with bbhash): "<<bit_per_kmer+4<<endl;
	cout<<"TOTAL Size estimated (MBytes): "<<(bit_per_kmer+4)*number_kmer/(8*1024*1024)<<endl;
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



vector<bool> kmer_Set_Light::get_seq(kmer mini,uint64_t pos,uint32_t n){
	return vector<bool>(bucketSeq.begin()+(all_buckets[mini].start+pos)*2,bucketSeq.begin()+(all_buckets[mini].start+pos+n)*2);
}



inline kmer kmer_Set_Light::update_kmer(uint64_t pos,kmer mini,kmer input){
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
		if(nuc==2){
			cout<<"T";
		}
		if(nuc==3){
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






void kmer_Set_Light::create_mphf_mem(uint begin_BC,uint end_BC){
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
				uint32_t i_kmer(1);
				for(uint j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
					if(position_super_kmers[all_mphf[BC].mphf_size+bucketSize]){
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



void kmer_Set_Light::create_mphf_disk(uint begin_BC,uint end_BC){
	#pragma omp parallel  num_threads(coreNumber)
		{
		uint largest_bucket_anchor(0);
		uint largest_bucket_nuc(0);
		#pragma omp for schedule(dynamic, number_bucket_per_mphf.value())
		for(uint BC=(begin_BC);BC<end_BC;++BC){
			uint64_t mphfSize(0);
			string name("_blkmers"+to_string(BC));
			if(all_buckets[BC].nuc_minimizer!=0){
				ofstream out(name,ofstream::binary);
				largest_bucket_nuc=max(largest_bucket_nuc,all_buckets[BC].nuc_minimizer);
				largest_bucket_nuc_all=max(largest_bucket_nuc_all,all_buckets[BC].nuc_minimizer);
				uint bucketSize(1);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				out.write(reinterpret_cast<char*>(&canon),sizeof(canon));
				mphfSize++;
				for(uint j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
					if(position_super_kmers[all_mphf[BC].mphf_size+bucketSize]){
						j+=k-1;
						if((j+k)<all_buckets[BC].nuc_minimizer){
							seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
							out.write(reinterpret_cast<char*>(&canon),sizeof(canon));
							bucketSize++;
							mphfSize++;
						}
					}else{
						seq=update_kmer(j+k,BC,seq);
						rcSeq=(rcb(seq,k));
						canon=(min_k(seq, rcSeq));
						out.write(reinterpret_cast<char*>(&canon),sizeof(canon));
						bucketSize++;
						mphfSize++;
					}
				}
				largest_bucket_anchor=max(largest_bucket_anchor,bucketSize);
			}
			 //~ #pragma omp critical(coute)
			//~ {
				//~ cout<<mphfSize<<	"|"<<flush;
			//~ }
			if((BC+1)%number_bucket_per_mphf==0 and mphfSize!=0){
				largest_MPHF=max(largest_MPHF,mphfSize);
				auto data_iterator = file_binary(name.c_str());
				all_mphf[BC/number_bucket_per_mphf].kmer_MPHF= new boomphf::mphf<kmer,hasher_t>(mphfSize,data_iterator,gammaFactor);
				remove(name.c_str());
				largest_bucket_anchor=0;
				largest_bucket_nuc=(0);
				mphfSize=0;
			}
		}

	}
}



void kmer_Set_Light::int_to_bool(uint n_bits_to_encode,uint64_t X, uint64_t pos,uint64_t start){
	for(uint64_t i(0);i<n_bits_to_encode;++i){
		uint64_t pos_check(i+pos*n_bits_to_encode+start);
		positions_mutex[pos_check*4096/positions.size()].lock();
		positions[pos_check]=X%2;
		positions_mutex[pos_check*4096/positions.size()].unlock();
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
	return res;
}



void kmer_Set_Light::fill_positions(uint begin_BC,uint end_BC){
	#pragma omp parallel for num_threads(coreNumber)
	for(uint BC=(begin_BC);BC<end_BC;++BC){
		uint32_t super_kmer_id(0);
		if(all_buckets[BC].nuc_minimizer>0){
			uint32_t kmer_id(1);
			int n_bits_to_encode(all_mphf[BC/number_bucket_per_mphf].bit_to_encode);
			kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
			int_to_bool(n_bits_to_encode,super_kmer_id/positions_to_check.value(),all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
			for(uint j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
				if(position_super_kmers[all_mphf[BC].mphf_size+kmer_id]){
					j+=k-1;
					super_kmer_id++;
					kmer_id++;
					if((j+k)<all_buckets[BC].nuc_minimizer){
						seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
						{
							int_to_bool(n_bits_to_encode,super_kmer_id/positions_to_check.value(),all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
						}
					}
				}else{
					seq=update_kmer(j+k,BC,seq);
					rcSeq=(rcb(seq,k));
					canon=(min_k(seq, rcSeq));
					kmer_id++;
					{
						int_to_bool(n_bits_to_encode,super_kmer_id/positions_to_check.value(),all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
					}
				}
			}
		}
	}
}



bool kmer_Set_Light::query_kmer_bool(kmer canon){
	kmer min(regular_minimizer(canon));
	return single_query(min,canon);
}



int64_t kmer_Set_Light::query_kmer_hash(kmer canon){
	return query_get_hash(canon,regular_minimizer(canon));
}



int64_t kmer_Set_Light::query_kmer_minitig(kmer canon){
	return query_get_rank_minitig(canon,regular_minimizer(canon));
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


vector<int64_t> kmer_Set_Light::query_sequence_minitig(const string& query){
	vector<int64_t> res;
	if(query.size()<k){
		return res;
	}
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
	uint i(0);
	canon=(min_k(seq, rcSeq));
	res.push_back(query_kmer_minitig(canon));
	for(;i+k<query.size();++i){
		updateK(seq,query[i+k]);
		updateRCK(rcSeq,query[i+k]);
		canon=(min_k(seq, rcSeq));
		res.push_back(query_kmer_minitig(canon));
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



bool kmer_Set_Light::multiple_minimizer_query_bool(kmer minimizer, kmer kastor,uint prefix_length,uint suffix_length){
	if(suffix_length>0){
		const Pow2<kmer> max_completion(2*suffix_length);
		minimizer*=max_completion.value();
		for(uint i(0);i<max_completion;++i){
			kmer poential_min(min(minimizer|i,rcb(minimizer|i,m1)));
			if(single_query(poential_min,kastor)){
				return true;
			}
		}
	}
	if(prefix_length>0){
		const Pow2<kmer> max_completion(2*(prefix_length));
		const Pow2<kmer> mask(2*(m1-prefix_length));
		for(uint i(0);i<max_completion;++i){
			kmer poential_min(min(minimizer|i*mask,rcb(minimizer|i*mask,m1)));
			if(single_query(poential_min,kastor)){
				return true;
			}
		}
	}
	return false;
}


int64_t kmer_Set_Light::multiple_minimizer_query_hash(kmer minimizer, kmer kastor,uint prefix_length,uint suffix_length){
	if(suffix_length>0){
		kmer max_completion(1);
		max_completion<<=(2*suffix_length);
		minimizer<<=(2*suffix_length);
		for(uint i(0);i<max_completion;++i){
			kmer poential_min(min(minimizer+i,rcb(minimizer+i,m1)));
			return query_get_hash(kastor,poential_min);
		}
	}
	if(prefix_length>0){
		kmer max_completion(1);
		kmer mask(1);
		max_completion<<=(2*(prefix_length));
		mask<<=(2*(m1-prefix_length));
		for(uint i(0);i<max_completion;++i){
			kmer poential_min(min(minimizer+i*mask,rcb(minimizer+i*mask,m1)));
			return query_get_hash(kastor,poential_min);
		}
	}
	return -1;
}



bool kmer_Set_Light::single_query(const kmer minimizer, kmer kastor){
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



uint kmer_Set_Light::multiple_query_optimized(kmer minimizer, const vector<kmer>& kmerV){
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



int32_t kmer_Set_Light::query_get_pos_unitig(const kmer canon,kmer minimizer){
	if(unlikely(all_mphf[minimizer/number_bucket_per_mphf].empty)){
		return -1;
	}


	uint64_t hash=(all_mphf[minimizer/number_bucket_per_mphf].kmer_MPHF->lookup(canon));
	if(unlikely(hash == ULLONG_MAX)){
		return -1;
	}

	int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf].bit_to_encode);
	uint64_t read_position(bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf].start)*positions_to_check.value());
	uint64_t rank(all_buckets[minimizer].skmer_number+read_position);
	bm::id_t pos;
	bool found = position_super_kmers.select(rank+1, pos, *(position_super_kmers_RS));
	for(uint check_super_kmer(0);check_super_kmer<positions_to_check.value();++check_super_kmer){
		bm::id_t next_position(position_super_kmers.get_next(pos));
		bm::id_t stop_position;

		if(next_position==0){
			 stop_position=bucketSeq.size();
		}else{
			stop_position =next_position+(rank+check_super_kmer)*(k-1)+1;
		}
		pos+=(rank+check_super_kmer)*(k-1);
		if(likely(((uint64_t)pos+k-1)<bucketSeq.size())){
			kmer seqR=get_kmer(0,pos);
			kmer rcSeqR, canonR;
			for(uint64_t j=(pos);j<stop_position;++j){
				rcSeqR=(rcb(seqR,k));
				canonR=(min_k(seqR, rcSeqR));
				if(canon==canonR){
					return j;
				}
				seqR=update_kmer(j+k,0,seqR);//can be avoided
			}
		}
		pos=next_position;
	}
	return -1;
}



int64_t kmer_Set_Light::query_get_hash(const kmer canon,kmer minimizer){
	number_query++;
	if(unlikely(all_mphf[minimizer/number_bucket_per_mphf].empty))
		return -1;

	uint64_t hash=(all_mphf[minimizer/number_bucket_per_mphf].kmer_MPHF->lookup(canon));
	if(unlikely(hash == ULLONG_MAX))
		return -1;

	int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf].bit_to_encode);

	uint64_t rank(all_buckets[minimizer].skmer_number+bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf].start)*positions_to_check.value());
	bm::id_t pos;

	bool found = position_super_kmers.select(rank+1, pos, *(position_super_kmers_RS));
	for(uint check_super_kmer(0);check_super_kmer<positions_to_check.value();++check_super_kmer){
		bm::id_t next_position(position_super_kmers.get_next(pos));
		bm::id_t stop_position;

		if(next_position==0){
			 stop_position=all_buckets[minimizer].nuc_minimizer-k+1;
		}else{
			stop_position =next_position+(rank+check_super_kmer)*(k-1)-1;
		}
		pos+=(rank+check_super_kmer)*(k-1)-1;
		if(likely(((uint64_t)pos+k-1)<all_buckets[minimizer].nuc_minimizer)){
			kmer seqR=get_kmer(minimizer,pos);
			kmer rcSeqR, canonR;
			for(uint64_t j=(pos);j<stop_position;++j){
				rcSeqR=(rcb(seqR,k));
				canonR=(min_k(seqR, rcSeqR));
				if(canon==canonR){
					return hash+all_mphf[minimizer/number_bucket_per_mphf].mphf_size;
				}
				seqR=update_kmer(j+k,minimizer,seqR);//can be avoided
			}
		}
		pos=next_position;
	}
	return -1;
}




int64_t kmer_Set_Light::query_get_rank_minitig(const kmer canon,uint minimizer){
	number_query++;
	if(unlikely(all_mphf[minimizer/number_bucket_per_mphf].empty))
		return -1;

	uint64_t hash=(all_mphf[minimizer/number_bucket_per_mphf].kmer_MPHF->lookup(canon));
	if(unlikely(hash == ULLONG_MAX))
		return -1;

	int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf].bit_to_encode);

	uint64_t rank(all_buckets[minimizer].skmer_number+bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf].start)*positions_to_check.value());
	bm::id_t pos;

	bool found = position_super_kmers.select(rank+1, pos, *(position_super_kmers_RS));
	for(uint check_super_kmer(0);check_super_kmer<positions_to_check.value();++check_super_kmer){
		bm::id_t next_position(position_super_kmers.get_next(pos));
		bm::id_t stop_position;

		if(next_position==0){
			 stop_position=all_buckets[minimizer].nuc_minimizer-k+1;
		}else{
			stop_position =next_position+(rank+check_super_kmer)*(k-1)-1;
		}
		pos+=(rank+check_super_kmer)*(k-1)-1;
		if(likely(((uint64_t)pos+k-1)<all_buckets[minimizer].nuc_minimizer)){
			kmer seqR=get_kmer(minimizer,pos);
			kmer rcSeqR, canonR;
			for(uint64_t j=(pos);j<stop_position;++j){
				rcSeqR=(rcb(seqR,k));
				canonR=(min_k(seqR, rcSeqR));
				if(canon==canonR){
					return rank+all_buckets[minimizer].skmer_number;
				}
				seqR=update_kmer(j+k,minimizer,seqR);//can be avoided
			}
		}
		pos=next_position;
	}
	return -1;
}





void kmer_Set_Light::file_query(const string& query_file){
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	auto in=new zstr::ifstream(query_file);
	atomic<uint64_t> TP(0),FP(0);
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
				TP+=pair.first;
				FP+=pair.second;
			}
		}
	}
	cout<<"-----------------------QUERY RECAP 2----------------------------"<<endl;
	cout<<"Good kmer: "<<intToString(TP)<<endl;
	cout<<"Erroneous kmers: "<<intToString(FP)<<endl;
	cout<<"Query performed: "<<intToString(FP+TP)<<endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole QUERY took me " << time_span.count() << " seconds."<< endl;
	delete in;
}


void dump_vector_bool(const vector<bool> V, ostream* out ){
	int cmp = 0;
	uint8_t output=0;
	vector<uint8_t> buf;
	for(uint i(0);i<V.size();++i){
		output=output | ((V[i] ? 1 : 0) << cmp);
		cmp++;
		if(cmp==8){
			buf.push_back(output);
			if(buf.size()==4000){
				out->write((char*)buf.data(),buf.size());
				buf.clear();
			}
			cmp = 0;
			output = 0;
		}
	}
	if(V.size()%8!=0){
		buf.push_back(output);
	}
    out->write((char*)buf.data(),buf.size());
}


void read_vector_bool(vector<bool>& V, ifstream* out, uint n_bits ){
	uint size_buffer(4000);
	uint n_bytes(n_bits/8+(n_bits%8==0 ? 0 :1));
	uint position(0);
	vector<uint8_t> buf(size_buffer,0);
	while(position+size_buffer<n_bytes){
		out->read((char*)buf.data(),size_buffer);
		for(uint i(0);i<buf.size();++i){
			V.push_back(buf[i] & 1);
			V.push_back(buf[i] & 2);
			V.push_back(buf[i] & 4);
			V.push_back(buf[i] & 8);
			V.push_back(buf[i] & 16);
			V.push_back(buf[i] & 32);
			V.push_back(buf[i] & 64);
			V.push_back(buf[i] & 128);
		}
		position+=size_buffer;
	}
	buf.resize(n_bytes-position,0);
	out->read((char*)buf.data(),n_bytes-position);
	for(uint i(0);i<buf.size();++i){
		V.push_back(buf[i] & 1);
		V.push_back(buf[i] & 2);
		V.push_back(buf[i] & 4);
		V.push_back(buf[i] & 8);
		V.push_back(buf[i] & 16);
		V.push_back(buf[i] & 32);
		V.push_back(buf[i] & 64);
		V.push_back(buf[i] & 128);
	}
}



void kmer_Set_Light::dump_disk(const string& output_file){
	filebuf fb;
	remove(output_file.c_str());
	fb.open (output_file, ios::out | ios::binary | ios::trunc);
	ostream out(&fb);

	bm::serializer<bm::bvector<> > bvs;
	bvs.byte_order_serialization(false);
	bvs.gap_length_serialization(false);
	bm::serializer<bm::bvector<> >::buffer sbuf;
	{
		unsigned char* buf = 0;
		bvs.serialize(position_super_kmers, sbuf);
		buf = sbuf.data();
		uint64_t sz = sbuf.size();
		auto point2 =&buf[0];
		out.write(reinterpret_cast<const char*>(&sz),sizeof(sz));
		out.write((char*)point2,sz);
	}


	out.write(reinterpret_cast<const char*>(&k),sizeof(k));
	out.write(reinterpret_cast<const char*>(&m1),sizeof(m1));
	out.write(reinterpret_cast<const char*>(&m3),sizeof(m3));
	out.write(reinterpret_cast<const char*>(&minimizer_size_graph),sizeof(minimizer_size_graph));
	out.write(reinterpret_cast<const char*>(&bit_saved_sub),sizeof(bit_saved_sub));

	uint64_t positions_size(positions.size()), bucketSeq_size(bucketSeq.size());
	out.write(reinterpret_cast<const char*>(&positions_size),sizeof(positions_size));
	out.write(reinterpret_cast<const char*>(&bucketSeq_size),sizeof(bucketSeq_size));


	for(uint i(0);i<mphf_number;++i){
		out.write(reinterpret_cast<const char*>(&all_mphf[i].mphf_size),sizeof(all_mphf[i].mphf_size));
		out.write(reinterpret_cast<const char*>(&all_mphf[i].empty),sizeof(all_mphf[i].empty));
		out.write(reinterpret_cast<const char*>(&all_mphf[i].start),sizeof(all_mphf[i].start));
		out.write(reinterpret_cast<const char*>(&all_mphf[i].bit_to_encode),sizeof(all_mphf[i].bit_to_encode));
		if(not all_mphf[i].empty){
			all_mphf[i].kmer_MPHF->save(out);
		}
	}

	dump_vector_bool(bucketSeq,&out);
	dump_vector_bool(positions,&out);


	for(uint i(0);i<minimizer_number;++i){
		out.write(reinterpret_cast<const char*>(&all_buckets[i].skmer_number),sizeof(all_buckets[i].skmer_number));
		out.write(reinterpret_cast<const char*>(&all_buckets[i].nuc_minimizer),sizeof(all_buckets[i].nuc_minimizer));
		out.write(reinterpret_cast<const char*>(&all_buckets[i].current_pos),sizeof(all_buckets[i].current_pos));
		out.write(reinterpret_cast<const char*>(&all_buckets[i].start),sizeof(all_buckets[i].start));
	}

	fb.close();
	cout<<"Index dump"<<endl;
}


kmer_Set_Light::kmer_Set_Light(const string& index_file){

	ifstream out(index_file);


	uint64_t sz;
	out.read(reinterpret_cast< char*>(&sz),sizeof(sz));
	uint8_t* buff= new uint8_t[sz];
	out.read((char*)buff,sz);

	bm::deserialize(position_super_kmers, buff);
	delete buff;


	out.read(reinterpret_cast< char*>(&k),sizeof(k));
	out.read(reinterpret_cast< char*>(&m1),sizeof(m1));
	m2=m1;
	coreNumber=20;
	uint64_t  bucketSeq_size, positions_size;
	out.read(reinterpret_cast< char*>(&m3),sizeof(m3));
	out.read(reinterpret_cast< char*>(&minimizer_size_graph),sizeof(minimizer_size_graph));
	out.read(reinterpret_cast< char*>(&bit_saved_sub),sizeof(bit_saved_sub));
	out.read(reinterpret_cast< char*>(&positions_size),sizeof(positions_size));
	out.read(reinterpret_cast< char*>(&bucketSeq_size),sizeof(bucketSeq_size));
	offsetUpdateAnchor=(2*k);
	offsetUpdateMinimizer=(2*m1);
	mphf_number=(2*m2);
	number_superbuckets=(2*m3);
	minimizer_number=(2*m1);
	minimizer_number_graph=(2*minimizer_size_graph);
	number_bucket_per_mphf=(2*(m1-m2));
	bucket_per_superBuckets=(2*(m1-m3));
	positions_to_check=(bit_saved_sub);
	gammaFactor=2;


	all_buckets=new bucket_minimizer[minimizer_number.value()]();
	all_mphf=new info_mphf[mphf_number.value()];
	for(uint i(0);i<mphf_number;++i){
		out.read(reinterpret_cast< char*>(&all_mphf[i].mphf_size),sizeof(all_mphf[i].mphf_size));
		out.read(reinterpret_cast< char*>(&all_mphf[i].empty),sizeof(all_mphf[i].empty));
		out.read(reinterpret_cast< char*>(&all_mphf[i].start),sizeof(all_mphf[i].start));
		out.read(reinterpret_cast< char*>(&all_mphf[i].bit_to_encode),sizeof(all_mphf[i].bit_to_encode));
		if(not all_mphf[i].empty){
			all_mphf[i].kmer_MPHF=new boomphf::mphf<kmer,hasher_t>;
			all_mphf[i].kmer_MPHF->load(out);
		}
	}


	read_vector_bool(bucketSeq,&out,bucketSeq_size);
	read_vector_bool(positions,&out,positions_size);


	for(uint i(0);i<minimizer_number;++i){
		out.read(reinterpret_cast< char*>(&all_buckets[i].skmer_number),sizeof(all_buckets[i].skmer_number));
		out.read(reinterpret_cast< char*>(&all_buckets[i].nuc_minimizer),sizeof(all_buckets[i].nuc_minimizer));
		out.read(reinterpret_cast< char*>(&all_buckets[i].current_pos),sizeof(all_buckets[i].current_pos));
		out.read(reinterpret_cast< char*>(&all_buckets[i].start),sizeof(all_buckets[i].start));
	}


	position_super_kmers_RS=new bm::bvector<>::rs_index_type();
	position_super_kmers.build_rs_index(position_super_kmers_RS);

	cout<<"Index loaded"<<endl;
}





