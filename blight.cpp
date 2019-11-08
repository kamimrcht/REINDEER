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



uint64_t max_super_kmer_size(0);

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
	for(uint64_t i(0);i<str.size();i++){
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



inline __uint128_t rcb(const __uint128_t& in, uint64_t n){

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



inline uint64_t rcb(uint64_t in, uint64_t n){
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



static inline kmer get_int_in_kmer(kmer seq,uint64_t pos,uint64_t number_nuc){
	seq>>=2*pos;
	return  ((seq)%(1<<(2*number_nuc)));
}



kmer canonize(kmer x, uint64_t k){
	return min(x,rcb(x,k));
	return ( (__builtin_popcountll(x) & 1) ? x : rcb(x, k)) >> 1;
}


void print_bin(uint64_t n){
	uint64_t mask=1;
	mask<<=63;
	for(uint64_t i(0);i<64;++i){
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
	for(uint64_t i(1);i<=k-minimizer_size_graph;i++){
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
			uint64_t i(0);
			for(;i+m1<ref.size();++i){
				updateM(seq,ref[i+m1]);
				updateRCM(rcSeq,ref[i+m1]);
				canon=(min_k(seq,rcSeq));
					abundance_minimizer_temp[canon]++;
			}
		}
	}
	for(uint64_t i(0);i<minimizer_number;++i){
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



void kmer_Set_Light::reset(){
	number_kmer=number_super_kmer=largest_MPHF=positions_total_size=positions_int=total_nb_minitigs=largest_bucket_nuc_all=0;
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
}



bool exists_test (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}



void kmer_Set_Light::construct_index_fof(const string& input_file){
	if(m1<m2){
		cout<<"n should be inferior to m"<<endl;
		exit(0);
	}
	if(m2<m3){
		cout<<"s should be inferior to n"<<endl;
		exit(0);
	}
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	ifstream infof(input_file);
	string file;
	uint64_t i_file(1);
	//STACK SUPERBUCKETS
	while(not infof.eof()){
		file.clear();
		getline(infof,file);
		if(exists_test(file)){
			create_super_buckets_regular(file,i_file++);
		}
	}
	cout<<"Partition created"<<endl;
	{
		zstr::ofstream out("_blmonocolor.fa.gz",ios::app);

		#pragma omp parallel for num_threads(coreNumber)
		for(uint i_superbuckets=0; i_superbuckets<number_superbuckets.value(); ++i_superbuckets){
			//SORT SUPERBUCKETS
			string cmd("sort _blout"+to_string(i_superbuckets)+" --output _blsout"+to_string(i_superbuckets));
			uint res(system(cmd.c_str()));
			if(res!=0){cout<<"Problem with sort command"<<endl;exit(0);}
			merge_super_buckets("_blsout"+to_string(i_superbuckets),i_file-1,&out);
			remove(("_blsout"+to_string(i_superbuckets)).c_str());

		}
	}
	reset();
	cout<<"Monocolor minitig computed, now regular indexing start"<<endl;
	create_super_buckets_regular("_blmonocolor.fa.gz");

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



uint16_t parseCoverage(const string& str){
	size_t pos(str.find("km:f:"));
	if(pos==string::npos){
		pos=(str.find("KM:f:"));
	}
	if(pos==string::npos){
		return 1;
	}
	uint i(1);
	while(str[i+pos+5]!=' '){
		++i;
	}
	return (uint16_t)stof(str.substr(pos+5,i));
}



void kmer_Set_Light::create_super_buckets_regular(const string& input_file,int dbg_id){
	uint64_t total_nuc_number(0);
	auto inUnitigs=new zstr::ifstream(input_file);
	if( not inUnitigs->good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	vector<ostream*> out_files;
	for(uint64_t i(0);i<number_superbuckets;++i){
		if(dbg_id==0){
			auto out =new zstr::ofstream("_blout"+to_string(i));
			out_files.push_back(out);
		}else{
			//FOR SPLITTER
			auto out =new ofstream("_blout"+to_string(i),ofstream::app);
			out_files.push_back(out);
		}

	}
	omp_lock_t lock[number_superbuckets.value()];
	for (uint64_t i=0; i<number_superbuckets; i++){
		omp_init_lock(&(lock[i]));
	}
	#pragma omp parallel num_threads(coreNumber)
	{
		string ref,useless;
		minimizer_type old_minimizer,minimizer;
		while(not inUnitigs->eof()){
			ref=useless="";
			#pragma omp critical(dataupdate)
			{
				getline(*inUnitigs,useless);
				getline(*inUnitigs,ref);
			}
			//FOREACH UNITIG
			if(not ref.empty() and not useless.empty()){
				old_minimizer=minimizer=minimizer_number.value();
				uint64_t last_position(0);
				//FOREACH KMER
				kmer seq(str2num(ref.substr(0,k)));
				minimizer=regular_minimizer(seq);
				old_minimizer=minimizer;
				uint64_t i(0);
				for(;i+k<ref.size();++i){
					updateK(seq,ref[i+k]);
					//COMPUTE KMER MINIMIZER
					minimizer=regular_minimizer(seq);
					if(old_minimizer!=minimizer){
						omp_set_lock(&(lock[((old_minimizer))/bucket_per_superBuckets]));
						if(dbg_id==0){
							*(out_files[((old_minimizer))/bucket_per_superBuckets])<<">"+to_string(old_minimizer)+"\n"<<ref.substr(last_position,i-last_position+k)<<"\n";
						}else{
							*(out_files[((old_minimizer))/bucket_per_superBuckets])<<""+to_string(old_minimizer)+":"<<to_string(dbg_id)<<":"<<ref.substr(last_position,i-last_position+k)<<":"<<parseCoverage(useless)<<"\n";
						}
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
					if(dbg_id==0){
						*(out_files[((old_minimizer))/bucket_per_superBuckets])<<">"+to_string(old_minimizer)+"\n"<<ref.substr(last_position)<<"\n";
					}else{
						*(out_files[((old_minimizer))/bucket_per_superBuckets])<<""+to_string(old_minimizer)+":"<<to_string(dbg_id)<<":"<<ref.substr(last_position)<<":"<<parseCoverage(useless)<<"\n";
					}
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
	for(uint64_t i(0);i<number_superbuckets;++i){
		*out_files[i]<<flush;
		delete(out_files[i]);
	}
	if(dbg_id!=0){return;}
	bucketSeq.resize(total_nuc_number*2);
	bucketSeq.shrink_to_fit();
	uint64_t i(0),total_pos_size(0);
	uint32_t max_bucket_mphf(0);
	uint64_t hash_base(0),old_hash_base(0), nb_skmer_before(0), last_skmer_number(0);
	for(uint64_t BC(0);BC<minimizer_number.value();++BC){
		all_buckets[BC].start=i;
		all_buckets[BC].current_pos=i;
		i+=all_buckets[BC].nuc_minimizer;
		max_bucket_mphf=max(all_buckets[BC].skmer_number,max_bucket_mphf);
		if (BC == 0){
			nb_skmer_before = 0; // I replace skmer_number by the total number of minitigs before this bucket
		}
		else{
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
	positions_int=positions.size()/64+(positions.size()%64==0 ? 0 :1);
	positions.shrink_to_fit();
}



void kmer_Set_Light::str2bool(const string& str,uint64_t mini){
	for(uint64_t i(0);i<str.size();++i){
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


vector<string> split(const string &s, char delim){
	stringstream ss(s);
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)) {
		elems.push_back(move(item));
	}
	return elems;
}

string bool2str(vector<bool> V){
	string result;
	for(uint64_t i(0);i<V.size();++i){
		result+=(V[i]?'1':'0');
	}
	return result;
}



string color_coverage2str(const vector<uint16_t>& V){
	string result;
	for(uint64_t i(0);i<V.size();++i){
		result+=":"+to_string(V[i]);
	}
	return result;
}



string kmer_Set_Light::compaction(const string& seq1,const string& seq2, bool recur=true){
	uint s1(seq1.size()),s2(seq2.size());
	if(s1==0 or s2==0){return "";}
	string rc2(revComp(seq2)),end1(seq1.substr(s1-k+1,k-1)), beg2(seq2.substr(0,k-1));
	if(end1==beg2){return seq1+(seq2.substr(k-1));}
	string begrc2(rc2.substr(0,k-1));
	if(end1==begrc2){return seq1+(rc2.substr(k-1));}
	if(recur){
		return compaction(revComp(seq1),seq2,false	);
	}else{
		//~ cout<<seq1<<" "<<seq2<<" "<<revComp(seq1)<<" "<<revComp(seq2)<<endl;
		//~ cout<<"STRANGER THINGS"<<endl;//TODO SOME COMPACTION ARE MISSED CAN WE DO BETTER
	}
	return "";
}



void kmer_Set_Light::get_monocolor_minitigs(const  vector<string>& minitigs, const vector<int64_t>& color, const  vector<uint16_t>& coverage, zstr::ofstream* out, const string& mini,uint64_t number_color){
	unordered_map<kmer,kmer,KmerHasher> next_kmer;
	unordered_map<kmer,kmer,KmerHasher> previous_kmer;
	unordered_map<kmer,vector<uint16_t>,KmerHasher> kmer_color;
	vector<uint16_t> bit_vector(number_color,0);
	//ASSOCIATE INFO TO KMERS
	for(uint64_t i_mini(0);i_mini<minitigs.size();++i_mini){
		kmer seq(str2num(minitigs[i_mini].substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq)),prev(-1);
		canon=(min_k(seq, rcSeq));
		if(kmer_color.count(canon)==0){
			kmer_color[canon]=bit_vector;
		}
		//~ print_kmer(canon);
		kmer_color[canon][color[i_mini]-1]=coverage[i_mini];

		for(uint i(0);i+k<minitigs[i_mini].size();++i){
			prev=canon;
			updateK(seq,minitigs[i_mini][i+k]);
			updateRCK(rcSeq,minitigs[i_mini][i+k]);
			canon=(min_k(seq, rcSeq));
			//~ print_kmer(canon);
			if(kmer_color.count(canon)==0){
				kmer_color[canon]=bit_vector;
			}
			kmer_color[canon][color[i_mini]-1]=coverage[i_mini];
			if(next_kmer.count(prev)==0){
				next_kmer[prev]=canon;
			}else{
				if(next_kmer[prev]!=canon and next_kmer[prev]!=-1){
					next_kmer[prev]=-1;
				}
			}
			if(previous_kmer.count(canon)==0){
				previous_kmer[canon]=prev;
			}else{
				if(previous_kmer[canon]!=prev and previous_kmer[canon]!=-1){
					previous_kmer[canon]=-1;
				}
			}
		}
	}
	string monocolor,compact;
	//~ cout<<"go max minitig"<<endl;
	uint kmer_number_check(0);
	//FOREACH KMER COMPUTE MAXIMAL MONOCOLOR MINITIG
	for (auto& it: kmer_color) {
		kmer_number_check++;
		if(not it.second.empty()){
			auto dumpcolor(it.second);
			kmer canon=it.first;
			monocolor=kmer2str(canon);
			//~ cout<<"go "<<monocolor<<endl;
			//GO FORWARD
			while(true){
				if(next_kmer.count(canon)==0){break;}
				if(next_kmer[canon]==-1){break;}
				if(kmer_color[next_kmer[canon]]!=it.second){break;}
				kmer_color[next_kmer[canon]].clear();
				compact=compaction(monocolor,kmer2str(next_kmer[canon]));
				if(compact.empty()){break;}
				monocolor=compact;
				canon=next_kmer[canon];
			}
			//GO BACKWARD
			canon=it.first;
			while(true){
				if(previous_kmer.count(canon)==0){break;}
				if(previous_kmer[canon]==-1){break;}
				if(kmer_color[previous_kmer[canon]]!=it.second){break;}
				kmer_color[previous_kmer[canon]].clear();
				compact=compaction(monocolor,kmer2str(previous_kmer[canon]));
				if(compact.empty()){break;}
				monocolor=compact;
				canon=previous_kmer[canon];
			}
			#pragma omp critical (monocolorFile)
			{
				*out<<">"+mini<<color_coverage2str(dumpcolor)<<"\n"<<monocolor<<"\n";
			}
		}
	}
}



void kmer_Set_Light::merge_super_buckets(const string& input_file, uint64_t number_color, zstr::ofstream* out){
	string line;
	vector<string> splitted;
	vector<string> minitigs;
	vector<int64_t> color;
	vector<uint16_t> coverage;
	zstr::ifstream in(input_file);
	int64_t minimizer=-1;
	while(not in.eof() and in.good()){
		getline(in,line);
		if(not line.empty()){
			splitted=split(line,':');
			if(stoi(splitted[0])>minimizer and not minitigs.empty()){
				get_monocolor_minitigs(minitigs,color,coverage,out,to_string(minimizer),number_color);
				minitigs.clear();
				color.clear();
				coverage.clear();
			}
			minimizer=stoi(splitted[0]);
			color.push_back(stoi(splitted[1]));
			minitigs.push_back(splitted[2]);
			coverage.push_back(stoi(splitted[3]));
		}
	}
	get_monocolor_minitigs(minitigs,color,coverage,out,to_string(minimizer),number_color);
}



void kmer_Set_Light::read_super_buckets(const string& input_file){
	uint64_t total_size(0);
	//~ #pragma omp parallel num_threads(coreNumber)
	{
		string useless,line;
		//~ #pragma omp for
		for(uint64_t SBC=0;SBC<number_superbuckets.value();++SBC){
			vector<uint64_t> number_kmer_accu(bucket_per_superBuckets.value(),0);
			uint64_t BC(SBC*bucket_per_superBuckets);
			//~ cout<<"go"<<BC<<endl;
			zstr::ifstream in((input_file+to_string(SBC)));
			while(not in.eof() and in.good()){
				useless=line="";
				getline(in,useless);
				getline(in,line);
				//~ if(line.size()>=k){cout<<line<<" "<<number_kmer<<" "<<BC<<endl;}
				if(not line.empty()){
					useless=useless.substr(1);
					uint64_t minimizer(stoi(useless));
					str2bool(line,minimizer);
					//~ #pragma omp critical(PSK)
					{
						position_super_kmers[number_kmer_accu[minimizer%bucket_per_superBuckets]+all_mphf[minimizer].mphf_size]=true;
					}
					number_kmer+=line.size()-k+1;
					number_kmer_accu[minimizer%bucket_per_superBuckets]+=line.size()-k+1;
					line.clear();
					number_super_kmer++;
				}
			}
			remove((input_file+to_string(SBC)).c_str());
			//~ create_mphf_mem(BC,BC+bucket_per_superBuckets);
			create_mphf_disk(BC,BC+bucket_per_superBuckets);
			position_super_kmers.optimize();
			position_super_kmers.optimize_gap_size();
			fill_positions(BC,BC+bucket_per_superBuckets);
			BC+=bucket_per_superBuckets.value();
			//~ cout<<"-"<<flush;
		}
	}
	position_super_kmers[number_kmer]=true;
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



vector<bool> kmer_Set_Light::get_seq(kmer mini,uint64_t pos,uint64_t n){
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



void kmer_Set_Light::print_kmer(kmer num,uint64_t n){
	Pow2<kmer> anc(2*(k-1));
	for(uint64_t i(0);i<k and i<n;++i){
		uint64_t nuc=num/anc;
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
	for(uint64_t i(0);i<k;++i){
		uint64_t nuc=num/anc;
		num=num%anc;
		if(nuc==3){
			res+="G";
		}
		if(nuc==2){
			res+="T";
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



void kmer_Set_Light::create_mphf_mem(uint64_t begin_BC,uint64_t end_BC){
	#pragma omp parallel  num_threads(coreNumber)
		{
		vector<kmer> anchors;
		uint32_t largest_bucket_anchor(0);
		uint32_t largest_bucket_nuc(0);
		#pragma omp for schedule(dynamic, number_bucket_per_mphf.value())
		for(uint64_t BC=(begin_BC);BC<end_BC;++BC){
			if(all_buckets[BC].nuc_minimizer!=0){
				largest_bucket_nuc=max(largest_bucket_nuc,all_buckets[BC].nuc_minimizer);
				largest_bucket_nuc_all=max(largest_bucket_nuc_all,all_buckets[BC].nuc_minimizer);
				uint32_t bucketSize(1);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				anchors.push_back(canon);
				for(uint64_t j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
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



void kmer_Set_Light::create_mphf_disk(uint64_t begin_BC,uint64_t end_BC){
	#pragma omp parallel  num_threads(coreNumber)
		{
		uint32_t largest_bucket_anchor(0);
		uint32_t largest_bucket_nuc(0);
		#pragma omp for schedule(dynamic, number_bucket_per_mphf.value())
		for(uint64_t BC=(begin_BC);BC<end_BC;++BC){
			uint64_t mphfSize(0);
			string name("_blkmers"+to_string(BC));
			if(all_buckets[BC].nuc_minimizer!=0){
				ofstream out(name,ofstream::binary | ofstream::trunc);
				largest_bucket_nuc=max(largest_bucket_nuc,all_buckets[BC].nuc_minimizer);
				largest_bucket_nuc_all=max(largest_bucket_nuc_all,all_buckets[BC].nuc_minimizer);
				uint32_t bucketSize(1);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				out.write(reinterpret_cast<char*>(&canon),sizeof(canon));
				mphfSize++;
				for(uint64_t j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
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



void kmer_Set_Light::int_to_bool(uint64_t n_bits_to_encode,uint64_t X, uint64_t pos,uint64_t start){
	for(uint64_t i(0);i<n_bits_to_encode;++i){
		uint64_t pos_check(i+pos*n_bits_to_encode+start);
		//~ positions_mutex[0].lock();
		positions_mutex[(pos_check/64)*4096/(positions_int)].lock();
		positions[pos_check]=X%2;
		//~ positions_mutex[0].unlock();
		positions_mutex[(pos_check/64)*4096/(positions_int)].unlock();
		X>>=1;
	}
}


uint64_t kmer_Set_Light::bool_to_int(uint64_t n_bits_to_encode,uint64_t pos,uint64_t start){
	uint64_t res(0);
	uint64_t acc(1);
	for(uint64_t i(0);i<n_bits_to_encode;++i, acc<<=1){
		if(positions[i+pos*n_bits_to_encode+start]){
			res |= acc;
		}
	}
	return res;
}



void kmer_Set_Light::fill_positions(uint64_t begin_BC,uint64_t end_BC){
	#pragma omp parallel for num_threads(coreNumber)
	for(uint64_t BC=(begin_BC);BC<end_BC;++BC){
		uint64_t super_kmer_id(0);
		if(all_buckets[BC].nuc_minimizer>0){
			uint64_t kmer_id(1);
			int n_bits_to_encode(all_mphf[BC/number_bucket_per_mphf].bit_to_encode);
			kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
			int_to_bool(n_bits_to_encode,super_kmer_id/positions_to_check.value(),all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
			for(uint64_t j(0);(j+k)<all_buckets[BC].nuc_minimizer;j++){
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



int64_t kmer_Set_Light::kmer_to_hash(const kmer canon,kmer minimizer){
	if(unlikely(all_mphf[minimizer/number_bucket_per_mphf].empty)){return -1;}
	uint64_t hash=(all_mphf[minimizer/number_bucket_per_mphf].kmer_MPHF->lookup(canon));
	if(unlikely(hash == ULLONG_MAX)){
		return -1;
	}else{
		return hash;
	}
}



int64_t kmer_Set_Light::hash_to_rank(const int64_t hash,kmer minimizer){
	int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf].bit_to_encode);
	int64_t rank(all_buckets[minimizer].skmer_number+bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf].start)*positions_to_check.value());
	return rank;
}



vector<kmer> kmer_Set_Light::kmer_to_superkmer(const kmer canon,kmer minimizer, int64_t& rank, int64_t& hash){
	hash=kmer_to_hash(canon,minimizer);
	if(hash<0){return {};}
	rank=(hash_to_rank(hash,minimizer));
	if(rank<0){return {};}
	vector<kmer> result;
	bool found(false);
	bm::id64_t pos,next_position,stop_position ;

	bool osef = position_super_kmers.select(rank+1, pos, *(position_super_kmers_RS));

	for(uint64_t check_super_kmer(0);check_super_kmer<positions_to_check.value() and not found;++check_super_kmer){
		next_position=(position_super_kmers.get_next(pos));
		if(next_position==0){
			cout<<"no succesor"<<endl;
			stop_position=bucketSeq.size()-k;
		}else{
			stop_position=next_position+(rank+check_super_kmer)*(k-1);
		}
		pos+=(rank+check_super_kmer)*(k-1);

		if(likely(((uint64_t)pos+k-1)<bucketSeq.size())){
			kmer seqR=get_kmer(0,pos);
			kmer rcSeqR, canonR;
			for(uint64_t j=(pos);j<stop_position;++j){
				rcSeqR=(rcb(seqR,k));
				canonR=(min_k(seqR, rcSeqR));
				result.push_back(canonR);
				if(canon==canonR){
					found=true;
				}
				if(likely(((j+k)*2<bucketSeq.size()))){
					seqR=update_kmer(j+k,0,seqR);//can be avoided
				}
			}
		}
		pos=next_position;
	}
	if(found){return result;}
	return{};
}



bool kmer_in_superkmer(const kmer canon,const vector<kmer>& V){
	for(uint64_t i(0);i<V.size();i++){
		if(canon==V[i]){return true;}
	}
	return false;
}



vector<bool> kmer_Set_Light::get_presence_query(const string& query){
	if(query.size()<k){return{};}
	vector<bool> result;
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
	int64_t i(0),rank,hash;
	vector<kmer> superKmers;
	canon=(min_k(seq, rcSeq));

	for(;i+k<=query.size();++i){
		if(superKmers.empty()){
			superKmers=kmer_to_superkmer(canon,regular_minimizer(canon),rank,hash);
			//~ if(superKmers.empty()){cout<<"WEIRD";cin.get();}
			result.push_back(not superKmers.empty());
		}else{
			if(kmer_in_superkmer(canon,superKmers)){
				//~ cout<<"in SK"<<endl;
				result.push_back(true);
			}else{
				superKmers=kmer_to_superkmer(canon,regular_minimizer(canon),rank,hash);
				//~ if(superKmers.empty()){cout<<"WEIRD";}
				result.push_back(not superKmers.empty());
			}
		}
		if(i+k<query.size()){
			updateK(seq,query[i+k]);
			updateRCK(rcSeq,query[i+k]);
			canon=(min_k(seq, rcSeq));
		}
	}
	return result;
}



vector<int64_t> kmer_Set_Light::get_rank_query(const string& query){
	if(query.size()<k){return{};}
	vector<int64_t> result;
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
	int64_t i(0),rank,hash;
	vector<kmer> superKmers;
	canon=(min_k(seq, rcSeq));
	kmer minimizer(regular_minimizer(canon));

	for(;i+k<=query.size();++i){
		if(superKmers.empty()){
			superKmers=kmer_to_superkmer(canon,minimizer,rank,hash);
			result.push_back(superKmers.empty()? -1 : rank);
		}else{
			if(kmer_in_superkmer(canon,superKmers)){
				result.push_back(rank);
			}else{
				superKmers=kmer_to_superkmer(canon,minimizer,rank,hash);
				result.push_back(superKmers.empty()? -1 : rank);
			}
		}
		if(i+k<query.size()){
			updateK(seq,query[i+k]);
			updateRCK(rcSeq,query[i+k]);
			canon=(min_k(seq, rcSeq));
			minimizer=regular_minimizer(canon);
		}
	}
	return result;
}



vector<int64_t> kmer_Set_Light::get_hashes_query(const string& query){
	if(query.size()<k){return{};}
	vector<int64_t> result;
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
	int64_t i(0),rank,hash;
	vector<kmer> superKmers;
	canon=(min_k(seq, rcSeq));
	kmer minimizer(regular_minimizer(canon));

	for(;i+k<=query.size();++i){
		if(superKmers.empty()){
			superKmers=kmer_to_superkmer(canon,minimizer,rank,hash);
			result.push_back(superKmers.empty() ? -1 : hash+all_mphf[minimizer/number_bucket_per_mphf].mphf_size);
		}else{
			if(kmer_in_superkmer(canon,superKmers)){
				result.push_back(kmer_to_hash(canon,minimizer)+all_mphf[minimizer/number_bucket_per_mphf].mphf_size);
			}else{
				superKmers=kmer_to_superkmer(canon,minimizer,rank,hash);
				result.push_back(superKmers.empty() ? -1 : hash+all_mphf[minimizer/number_bucket_per_mphf].mphf_size);
			}
		}
		if(i+k<query.size()){
			updateK(seq,query[i+k]);
			updateRCK(rcSeq,query[i+k]);
			canon=(min_k(seq, rcSeq));
			minimizer=(regular_minimizer(canon));
		}
	}
	return result;
}



void kmer_Set_Light::file_query_presence(const string& query_file){
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	auto in=new zstr::ifstream(query_file);
	atomic<uint64_t> TP(0),FP(0);
	#pragma omp parallel num_threads(coreNumber)
	{
		while(not in->eof() and in->good()){
			string query;
			vector<bool> presences;
			#pragma omp critical(dataupdate)
			{
				getline(*in,query);
				getline(*in,query);
			}
			if(query.size()>=k){
				presences=get_presence_query(query);
				uint64_t numberOne=count(presences.begin(), presences.end(), true);
				TP+=numberOne;
				FP+=presences.size()-numberOne;
				presences.clear();
			}
		}
	}
	cout<<endl<<"-----------------------QUERY PRESENCE RECAP ----------------------------"<<endl;
	cout<<"Good kmer: "<<intToString(TP)<<endl;
	cout<<"Erroneous kmers: "<<intToString(FP)<<endl;
	cout<<"Query performed: "<<intToString(FP+TP)<<endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole QUERY took me " << time_span.count() << " seconds."<< endl;
	delete in;
}



void kmer_Set_Light::file_query_hases(const string& query_file, bool check){
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	auto in=new zstr::ifstream(query_file);
	vector<int64_t> all_hashes;
	#pragma omp parallel num_threads(coreNumber)
	{
		vector<kmer> kmerV;
		while(not in->eof() and in->good()){
			string query;
			vector<int64_t> query_hashes;
			#pragma omp critical(query_file)
			{
				getline(*in,query);
				getline(*in,query);
			}
			if(query.size()>=k){
				query_hashes=get_hashes_query(query);
				if(check){
					#pragma omp critical(query_file)
					{
						all_hashes.insert( all_hashes.end(), query_hashes.begin(), query_hashes.end() );
					}
				}
			}
		}
	}
	cout<<endl<<"-----------------------QUERY HASHES RECAP----------------------------"<<endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole QUERY took me " << time_span.count() << " seconds."<< endl;
	if(not check){return;}
	sort(all_hashes.begin(),all_hashes.end());
	bool bijective(true);
	for(uint i(0);i<all_hashes.size();++i){
		if(all_hashes[i]!=i){
			//~ cout<<all_hashes[i]<<":"<<i<<"	";
			bijective=false;
		}
	}
	if(bijective){
		cout<<"HASHES ARE BIJECTIVE WITH [0::"<<all_hashes.size()-1<<"]"<<endl;
	}else{
		cout<<"HASHES ARE NOT BIJECTIVE"<<endl;
	}
	delete in;
}



void kmer_Set_Light::file_query_rank(const string& query_file){
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	auto in=new zstr::ifstream(query_file);
	vector<int64_t> all_hashes;
	#pragma omp parallel num_threads(coreNumber)
	{
		vector<kmer> kmerV;
		while(not in->eof() and in->good()){
			string query;
			vector<int64_t> query_hashes;
			#pragma omp critical(query_file)
			{
				getline(*in,query);
				getline(*in,query);
			}
			if(query.size()>=k){
				query_hashes=get_rank_query(query);
				#pragma omp critical(query_file)
				{
					all_hashes.insert( all_hashes.end(), query_hashes.begin(), query_hashes.end() );
				}
			}
		}
	}
	cout<<endl<<"-----------------------QUERY RANK RECAP----------------------------"<<endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole QUERY took me " << time_span.count() << " seconds."<< endl;
	sort(all_hashes.begin(),all_hashes.end());
	bool surjective(true);
	all_hashes	.erase( unique( all_hashes.begin(), all_hashes.end() ), all_hashes.end() );

	for(uint i(0);i<all_hashes.size();++i){
		if(all_hashes[i]!=i){
			//~ cout<<all_hashes[i]<<":"<<i<<"	";
			surjective=false;
		}
	}
	if(surjective){
		cout<<"RANKS ARE SURJECTIVE WITH [0::"<<all_hashes.size()-1<<"]"<<endl;
	}else{
		cout<<"RANKS ARE NOT BIJECTIVE"<<endl;
	}
	delete in;
}



void kmer_Set_Light::file_query_all_test(const string& query_file, bool full){
	file_query_presence(query_file);
	if(not full){return;}
	file_query_hases(query_file);

	file_query_rank(query_file);
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
	cout<<endl<<"-----------------------OLD QUERY RECAP----------------------------"<<endl;
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
	for(uint64_t i(0);i<V.size();++i){
		output=output | ((V[i] ? 1 : 0) << cmp);
		cmp++;
		if(cmp==8){
			buf.push_back(output);
			if(buf.size()==8000){
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



void read_vector_bool(vector<bool>& V, zstr::ifstream* out, uint64_t n_bits ){
	uint64_t size_buffer(8000);
	uint64_t n_bytes(n_bits/8+(n_bits%8==0 ? 0 :1));
	uint64_t position(0);
	vector<uint8_t> buf(size_buffer,0);
	while(position+size_buffer<n_bytes){
		out->read((char*)buf.data(),size_buffer);
		for(uint64_t i(0);i<buf.size();++i){
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
	for(uint64_t i(0);i<buf.size();++i){
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
	//OPEN
	filebuf fb;
	remove(output_file.c_str());
	fb.open (output_file, ios::out | ios::binary | ios::trunc);
	zstr::ostream out(&fb);


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


	for(uint64_t i(0);i<mphf_number;++i){
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


	for(uint64_t i(0);i<minimizer_number;++i){
		out.write(reinterpret_cast<const char*>(&all_buckets[i].skmer_number),sizeof(all_buckets[i].skmer_number));
		out.write(reinterpret_cast<const char*>(&all_buckets[i].nuc_minimizer),sizeof(all_buckets[i].nuc_minimizer));
		out.write(reinterpret_cast<const char*>(&all_buckets[i].current_pos),sizeof(all_buckets[i].current_pos));
		out.write(reinterpret_cast<const char*>(&all_buckets[i].start),sizeof(all_buckets[i].start));
	}
	out<<flush;
	fb.close();
	cout<<"Index dump"<<endl;
}



kmer_Set_Light::kmer_Set_Light(const string& index_file){

	zstr::ifstream out(index_file);

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
	for(uint64_t i(0);i<mphf_number;++i){
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

	for(uint64_t i(0);i<minimizer_number;++i){
		out.read(reinterpret_cast< char*>(&all_buckets[i].skmer_number),sizeof(all_buckets[i].skmer_number));
		out.read(reinterpret_cast< char*>(&all_buckets[i].nuc_minimizer),sizeof(all_buckets[i].nuc_minimizer));
		out.read(reinterpret_cast< char*>(&all_buckets[i].current_pos),sizeof(all_buckets[i].current_pos));
		out.read(reinterpret_cast< char*>(&all_buckets[i].start),sizeof(all_buckets[i].start));
	}

	position_super_kmers_RS=new bm::bvector<>::rs_index_type();
	position_super_kmers.build_rs_index(position_super_kmers_RS);

	cout<<"Index loaded"<<endl;
}





