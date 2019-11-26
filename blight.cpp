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
#include "utils.h"
#include "zstr.hpp"
#include "common.h"



using namespace std;
using namespace chrono;



uint kmer_number_check(0);
uint kmer_number_check2(0);
uint kmer_number_check3(0);
uint skmer_number_check(0);



inline __uint128_t kmer_Set_Light::rcb(const __uint128_t& in){

	assume(k <= 64, "k=%u > 64", k);

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

	res.m128i = mm_bitshift_right(res.m128i, 128 - 2*k);

	return res.k;
}



inline uint64_t  kmer_Set_Light::rcb(const uint64_t& in){
    assume(k <= 32, "k=%u > 32", k);
    // Complement, swap byte order
    uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
    // Swap nuc order in bytes
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
    const uint64_t c2 = 0x3333333333333333;
    res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
    res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

    // Realign to the right
    res >>= 64 - 2 * k;
    return res;
}



uint64_t kmer_Set_Light::canonize(uint64_t x,uint64_t n){
	return min(x,rcbc(x,n));
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
			kmer seq(str2num(ref.substr(0,m1))),rcSeq(rcbc(seq,m1)),canon(min_k(seq,rcSeq));
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

	nuc_minimizer=new uint32_t[minimizer_number.value()];
	current_pos=new uint64_t[minimizer_number.value()];
	start_bucket=new uint64_t[minimizer_number.value()];
	create_super_buckets(input_file);

	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Super bucket created: "<< time_span12.count() << " seconds."<<endl;

	read_super_buckets("_blout");

	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);
	cout<<"Indexes created: "<< time_span13.count() << " seconds."<<endl;
	duration<double> time_spant = duration_cast<duration<double>>(t13 - t1);
	cout << "The whole indexing took me " << time_spant.count() << " seconds."<< endl;
	delete [] nuc_minimizer;
	delete [] start_bucket;
	delete [] current_pos;
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
		//~ all_buckets[i].start=0;
		nuc_minimizer[i]=0;
		start_bucket[i]=0;
		current_pos[i]=0;
		all_buckets[i].skmer_number=0;
	}
}



void kmer_Set_Light::construct_index_fof(const string& input_file, bool countcolor, double max_divergence){
	count_color=countcolor;
	if(count_color){
		max_divergence_count=(max_divergence/100)+1;
	}
	cout<<max_divergence_count<<endl;
	if(m1<m2){
		cout<<"n should be inferior to m"<<endl;
		exit(0);
	}
	if(m2<m3){
		cout<<"s should be inferior to n"<<endl;
		exit(0);
	}
	nuc_minimizer=new uint32_t[minimizer_number.value()];
	current_pos=new uint64_t[minimizer_number.value()];
	start_bucket=new uint64_t[minimizer_number.value()];
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	ifstream infof(input_file);
	string file;
	uint64_t i_file(1);
	//STACK SUPERBUCKETS
	while(not infof.eof()){
		file.clear();
		getline(infof,file);
		if(exists_test(file)){
			create_super_buckets(file,i_file++);
		}
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	cout<<"Partition created	"<<read_kmer<<" kmers read "<<endl;
	duration<double> time_span12 = duration_cast<duration<double>>(t2 - t1);
	cout<<time_span12.count() << " seconds."<<endl;

	{
		zstr::ofstream out("_blmonocolor.fa.gz",ios::app);
		#pragma omp parallel for num_threads(coreNumber)
		for(uint i_superbuckets=0; i_superbuckets<number_superbuckets.value(); ++i_superbuckets){
			//SORT SUPERBUCKETS
			merge_super_buckets_mem("_blout"+to_string(i_superbuckets),i_file-1,&out);
			remove(("_blsout"+to_string(i_superbuckets)).c_str());
			cout<<"-"<<flush;
		}
	}
	cout<<endl;
	reset();
	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	cout<<"Monocolor minitig computed, now regular indexing start"<<endl;
	duration<double> time_span32 = duration_cast<duration<double>>(t3 - t2);
	cout<<time_span32.count() << " seconds."<<endl;
	create_super_buckets("_blmonocolor.fa.gz");

	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	duration<double> time_span43 = duration_cast<duration<double>>(t4 - t3);
	cout<<"Super buckets created: "<< time_span43.count() << " seconds."<<endl;

	read_super_buckets("_blout");
	delete [] nuc_minimizer ;
	delete [] start_bucket ;
	delete [] current_pos ;

	high_resolution_clock::time_point t5 = high_resolution_clock::now();
	duration<double> time_span53 = duration_cast<duration<double>>(t5 - t3);
	cout<<"Indexes created: "<< time_span53.count() << " seconds."<<endl;
	duration<double> time_spant = duration_cast<duration<double>>(t5 - t1);
	cout << "The whole indexing took me " << time_spant.count() << " seconds."<< endl;
}



void kmer_Set_Light::create_super_buckets(const string& input_file,int dbg_id){
	struct rlimit rl;
	getrlimit (RLIMIT_NOFILE, &rl);
	rl.rlim_cur = number_superbuckets.value()+10;
	setrlimit (RLIMIT_NOFILE, &rl);
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
			auto out =new  zstr::ofstream("_blout"+to_string(i),ofstream::app);
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
				read_kmer+=ref.size()-k+1;
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
						nuc_minimizer[old_minimizer]+=(i-last_position+k);
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
					nuc_minimizer[old_minimizer]+=(ref.substr(last_position)).size();
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
	uint64_t hash_base(0),old_hash_base(0),nb_skmer_before(0), last_skmer_number(0);
	for(uint64_t BC(0);BC<minimizer_number.value();++BC){
		start_bucket[BC]=i;
		current_pos[BC]=i;
		i+=nuc_minimizer[BC];
		max_bucket_mphf=max(all_buckets[BC].skmer_number,max_bucket_mphf);
		if (BC == 0){
			nb_skmer_before = 0; // I replace skmer_number by the total number of minitigs before this bucket
		}else{
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
			total_pos_size+=(n_bits_to_encode*all_mphf[BC/number_bucket_per_mphf].mphf_size);
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
				bucketSeq[(current_pos[mini]+i)*2]=(false);
				bucketSeq[(current_pos[mini]+i)*2+1]=(false);
				break;
			case 'C':
				bucketSeq[(current_pos[mini]+i)*2]=(false);
				bucketSeq[(current_pos[mini]+i)*2+1]=(true);
				break;
			case 'G':
				bucketSeq[(current_pos[mini]+i)*2]=(true);
				bucketSeq[(current_pos[mini]+i)*2+1]=(true);
				break;
			default:
				bucketSeq[(current_pos[mini]+i)*2]=(true);
				bucketSeq[(current_pos[mini]+i)*2+1]=(false);
				break;
			}
	}
	current_pos[mini]+=(str.size());
}

uint32_t misscompaction(0);

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
		cout<<++misscompaction<<endl;cin.get();
		//TODO SOME COMPACTION ARE MISSED CAN WE DO BETTER
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
		kmer seq(str2num(minitigs[i_mini].substr(0,k))),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq)),prev(-1);
		canon=(min_k(seq, rcSeq));
		if(kmer_color.count(canon)==0){
			kmer_color[canon]=bit_vector;
		}
		kmer_color[canon][color[i_mini]-1]=coverage[i_mini];

		for(uint i(0);i+k<minitigs[i_mini].size();++i){
			prev=canon;
			updateK(seq,minitigs[i_mini][i+k]);
			updateRCK(rcSeq,minitigs[i_mini][i+k]);
			canon=(min_k(seq, rcSeq));
			if(kmer_color.count(canon)==0){
				kmer_color[canon]=bit_vector;
			}
			kmer_color[canon][color[i_mini]-1]=coverage[i_mini];
			if(next_kmer.count(prev)==0){
				next_kmer[prev]=canon;
			}else{
				if(next_kmer[prev]!=canon and next_kmer[prev]!=(kmer)-1){
					next_kmer[prev]=(kmer)-1;
				}
			}
			if(previous_kmer.count(canon)==0){
				previous_kmer[canon]=prev;
			}else{
				if(previous_kmer[canon]!=prev and previous_kmer[canon]!=(kmer)-1){
					previous_kmer[canon]=(kmer)-1;
				}
			}
		}
	}
	string monocolor,compact;
	uint kmer_number_check(0);
	//FOREACH KMER COMPUTE MAXIMAL MONOCOLOR MINITIG
	for (auto& it: kmer_color) {
		kmer_number_check++;
		if(not it.second.empty()){
			auto dumpcolor(it.second);
			kmer canon=it.first;
			monocolor=kmer2str(canon);
			//GO FORWARD
			while(true){
				if(next_kmer.count(canon)==0){break;}
				if(next_kmer[canon]==(kmer)-1){break;}
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
				if(previous_kmer[canon]==(kmer)-1){break;}
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





bool equal_nonull(const vector<uint16_t>& V1,const vector<uint16_t>& V2){
	if(V1.empty()){return false;}
	if(V2.size()!=V1.size()){return false;}
	for(uint i(0);i<V1.size();++i){
		if((V1[i]==0)^(V2[i]==0)){
			return false;
		}
	}
	return true;
}



bool kmer_Set_Light::similar_count(const vector<uint16_t>& V1,const vector<uint16_t>& V2){
	if(V1.empty()){return false;}
	if(V2.size()!=V1.size()){return false;}
	for(uint i(0);i<V1.size();++i){
		if((double)max(V1[i],V2[i])/((double)min(V1[i],V2[i])>max_divergence_count)){
			return false;
		}
	}
	return true;
}





struct kmer_context{
	bool isdump;
	bool nextOK;
	bool prevOK;
	kmer next_kmer;
	kmer previous_kmer;
	vector<uint16_t> count;
};



void kmer_Set_Light::get_monocolor_minitigs_mem(const  vector<minitig>& minitigs , zstr::ofstream* out, const string& mini,uint64_t number_color){
	unordered_map<kmer,kmer_context> kmer2context;
	vector<uint16_t> bit_vector(number_color,0);
	string sequence;
	for(uint64_t i_mini(0);i_mini<minitigs.size();++i_mini){
		sequence=bool2strv(minitigs[i_mini].sequence);
		kmer seq(str2num(sequence.substr(0,k))),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq)),prev;
		canon=(min_k(seq, rcSeq));
		if(kmer2context.count(canon)==0){
			kmer2context[canon]={false,false,false,(kmer)-1,(kmer)-1,bit_vector};
		}
		kmer2context[canon].count[minitigs[i_mini].color-1]=minitigs[i_mini].coverage;
		for(uint i(0);i+k<sequence.size();++i){
			prev=canon;
			updateK(seq,sequence[i+k]);
			updateRCK(rcSeq,sequence[i+k]);
			canon=(min_k(seq, rcSeq));
			if(kmer2context.count(canon)==0){
				kmer2context[canon]={false,false,false,(kmer)-1,(kmer)-1,bit_vector};
			}
			kmer2context[canon].count[minitigs[i_mini].color-1]=minitigs[i_mini].coverage;
			if(not kmer2context[prev].nextOK){
				kmer2context[prev].next_kmer=canon;
				kmer2context[prev].nextOK=true;
			}
			if(kmer2context[canon].prevOK==0){
				kmer2context[canon].previous_kmer=prev;
				kmer2context[canon].prevOK=true;
			}
		}
	}

	string seq2dump,compact;
	for (auto&& it: kmer2context){
		if(not it.second.isdump){
			it.second.isdump=true;
			//~ kmer_number_check2++;
			auto colorV2dump(it.second.count);
			kmer canon=it.first;
			seq2dump=kmer2str(canon);

			while(true){
				if( not kmer2context[canon].nextOK){break;}
				kmer next(kmer2context[canon].next_kmer);
				if(kmer2context[next].isdump){break;}
				if(count_color){
					if(max_divergence_count==0){
						if(kmer2context[next].count!=colorV2dump){break;}
					}else{
						if(not similar_count(kmer2context[canon].count,colorV2dump)){break;}
					}
				}else{
					if(not equal_nonull(kmer2context[canon].count,colorV2dump)){break;}
				}
				compact=compaction(seq2dump,kmer2str(next));
				if(compact.empty()){break;}
				kmer2context[next].isdump=true;
				//~ kmer_number_check2++;
				seq2dump=compact;
				canon=next;
			}
			canon=it.first;
			while(true){
				if(not kmer2context[canon].prevOK){break;}
				kmer prev(kmer2context[canon].previous_kmer);
				if(kmer2context[prev].isdump){break;}
				if(count_color){
					if(kmer2context[prev].count!=colorV2dump){break;}
				}else{
					if(not equal_nonull(kmer2context[canon].count,colorV2dump)){break;}
				}
				compact=compaction(seq2dump,kmer2str(prev));
				if(compact.empty()){break;}
				kmer2context[prev].isdump=true;
				//~ kmer_number_check2++;
				seq2dump=compact;
				canon=prev;
			}

			#pragma omp critical (monocolorFile)
			{
				*out<<">"+mini<<color_coverage2str(colorV2dump)<<"\n"<<seq2dump<<"\n";
				//~ kmer_number_check+=seq2dump.size()-k+1;
				//~ skmer_number_check++;
				//~ if(kmer_number_check2!=kmer_number_check){
					//~ cout<<"STOP"<<endl;
					//~ cin.get();
				//~ }
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



void kmer_Set_Light::merge_super_buckets_mem(const string& input_file, uint64_t number_color, zstr::ofstream* out){
	string line;
	vector<string> splitted;
	vector<vector<minitig>> minimizer_to_minitigs(minimizer_number.value()/number_superbuckets.value());
	vector<int32_t> minimizers(minimizer_to_minitigs.size(),-1);
	zstr::ifstream in(input_file);
	while(not in.eof() and in.good()){
		getline(in,line);
		if(not line.empty()){
			splitted=split(line,':');
			minitig mini;
			mini.color=stoi(splitted[1]);
			mini.sequence=str2boolv(splitted[2]);
			mini.coverage=(stoi(splitted[3]));
			minimizer_to_minitigs[stoi(splitted[0])%minimizer_to_minitigs.size()].push_back(mini);
			minimizers[stoi(splitted[0])%minimizer_to_minitigs.size()]=(stoi(splitted[0]));
		}
	}
	#pragma omp parallel for num_threads(coreNumber)
	for(uint i=(0);i<minimizer_to_minitigs.size();i++){
		if(minimizers[i]!=-1){
			get_monocolor_minitigs_mem(minimizer_to_minitigs[i],out,to_string(minimizers[i]),number_color);
		}
	}
}



void kmer_Set_Light::read_super_buckets(const string& input_file){
	//~ #pragma omp parallel num_threads(coreNumber)
	{
		string useless,line;
		//~ #pragma omp for
		for(uint64_t SBC=0;SBC<number_superbuckets.value();++SBC){
			vector<uint64_t> number_kmer_accu(bucket_per_superBuckets.value(),0);
			uint64_t BC(SBC*bucket_per_superBuckets);
			zstr::ifstream in((input_file+to_string(SBC)));
			while(not in.eof() and in.good()){
				useless=line="";
				getline(in,useless);
				getline(in,line);
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
			cout<<"-"<<flush;
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
	uint64_t bit = (start_bucket[mini]+pos)*2;
	const uint64_t bitlast = bit + 2*k;
	for(;bit<bitlast;bit+=2){
		res<<=2;
		res |= bucketSeq[bit]*2 | bucketSeq[bit+1];
	}
	return res;
}



inline kmer kmer_Set_Light::get_kmer(uint64_t pos){
	kmer res(0);
	uint64_t bit = (pos)*2;
	const uint64_t bitlast = bit + 2*k;
	for(;bit<bitlast;bit+=2){
		res<<=2;
		res |= bucketSeq[bit]*2 | bucketSeq[bit+1];
	}
	return res;
}




inline kmer kmer_Set_Light::update_kmer(uint64_t pos,kmer mini,kmer input){
	return update_kmer_local(start_bucket[mini]+pos, bucketSeq, input);
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
			//~ if(all_buckets[BC].nuc_minimizer!=0){
			if(nuc_minimizer[BC]!=0){
				largest_bucket_nuc=max(largest_bucket_nuc,nuc_minimizer[BC]);
				largest_bucket_nuc_all=max(largest_bucket_nuc_all,nuc_minimizer[BC]);
				uint32_t bucketSize(1);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq));
				anchors.push_back(canon);
				for(uint64_t j(0);(j+k)<nuc_minimizer[BC];j++){
					if(position_super_kmers[all_mphf[BC].mphf_size+bucketSize]){
						j+=k-1;
						if((j+k)<nuc_minimizer[BC]){
							seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq)),canon=(min_k(seq,rcSeq));
							anchors.push_back(canon);
							bucketSize++;
						}
					}else{
						seq=update_kmer(j+k,BC,seq);
						rcSeq=(rcb(seq));
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
			if(nuc_minimizer[BC]!=0){
				ofstream out(name,ofstream::binary | ofstream::trunc);
				largest_bucket_nuc=max(largest_bucket_nuc,nuc_minimizer[BC]);
				largest_bucket_nuc_all=max(largest_bucket_nuc_all,nuc_minimizer[BC]);
				uint32_t bucketSize(1);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq));
				out.write(reinterpret_cast<char*>(&canon),sizeof(canon));
				mphfSize++;
				for(uint64_t j(0);(j+k)<nuc_minimizer[BC];j++){
					if(position_super_kmers[all_mphf[BC].mphf_size+bucketSize]){
						j+=k-1;
						if((j+k)<nuc_minimizer[BC]){
							seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq)),canon=(min_k(seq,rcSeq));
							out.write(reinterpret_cast<char*>(&canon),sizeof(canon));
							bucketSize++;
							mphfSize++;
						}
					}else{
						seq=update_kmer(j+k,BC,seq);
						rcSeq=(rcb(seq));
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
		positions_mutex[(pos_check/64)*4096/(positions_int)].lock();
		positions[pos_check]=X%2;
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
		if(nuc_minimizer[BC]>0){
			uint64_t kmer_id(1);
			int n_bits_to_encode(all_mphf[BC/number_bucket_per_mphf].bit_to_encode);
			kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq));
			int_to_bool(n_bits_to_encode,super_kmer_id/positions_to_check.value(),all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
			for(uint64_t j(0);(j+k)<nuc_minimizer[BC];j++){
				if(position_super_kmers[all_mphf[BC].mphf_size+kmer_id]){
					j+=k-1;
					super_kmer_id++;
					kmer_id++;
					if((j+k)<nuc_minimizer[BC]){
						seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq)),canon=(min_k(seq,rcSeq));
						{
							int_to_bool(n_bits_to_encode,super_kmer_id/positions_to_check.value(),all_mphf[BC/number_bucket_per_mphf].kmer_MPHF->lookup(canon),all_mphf[BC/number_bucket_per_mphf].start);
						}
					}
				}else{
					seq=update_kmer(j+k,BC,seq);
					rcSeq=(rcb(seq));
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

	position_super_kmers.select(rank+1, pos, *(position_super_kmers_RS));

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
			kmer seqR=get_kmer(pos);
			kmer rcSeqR, canonR;
			for(uint64_t j=(pos);j<stop_position;++j){
				rcSeqR=(rcb(seqR));
				canonR=(min_k(seqR, rcSeqR));
				result.push_back(canonR);
				if(canon==canonR){
					found=true;
				}
				if(likely(((j+k)*2<bucketSeq.size()))){
					//~ seqR=update_kmer(j+k,0,seqR);//can be avoided
					seqR=update_kmer_local(j+k, bucketSeq, seqR);
				}
			}
		}
		pos=next_position;
	}
	if(found){return result;}
	return{};
}



vector<bool> kmer_Set_Light::get_presence_query(const string& query){
	if(query.size()<k){return{};}
	vector<bool> result;
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq));
	int64_t i(0),rank,hash;
	vector<kmer> superKmers;
	canon=(min_k(seq, rcSeq));

	for(;i+k<=query.size();++i){
		if(superKmers.empty()){
			superKmers=kmer_to_superkmer(canon,regular_minimizer(canon),rank,hash);
			result.push_back(not superKmers.empty());
		}else{
			if(kmer_in_superkmer(canon,superKmers)){
				result.push_back(true);
			}else{
				superKmers=kmer_to_superkmer(canon,regular_minimizer(canon),rank,hash);
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
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq));
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
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq));
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
			cout<<all_hashes[i]<<":"<<i<<"	";
			bijective=false;
			break;
		}
	}
	if(bijective){
		cout<<"HASHES ARE BIJECTIVE WITH [0::"<<all_hashes.size()-1<<"]"<<endl;
	}else{
		cout<<"HASHES ARE NOT BIJECTIVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
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
	all_hashes.erase( unique( all_hashes.begin(), all_hashes.end() ), all_hashes.end() );

	for(uint i(0);i<all_hashes.size();++i){
		if(all_hashes[i]!=i){
			//~ cout<<all_hashes[i]<<":"<<i<<"	";
			surjective=false;
			break;
		}
	}
	if(surjective){
		cout<<"RANKS ARE SURJECTIVE WITH [0::"<<all_hashes.size()-1<<"]"<<endl;
	}else{
		cout<<"RANKS ARE NOT SURJECTIVE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	}
	delete in;
}



void kmer_Set_Light::file_query_all_test(const string& query_file, bool full){
	file_query_presence(query_file);
	if(not full){return;}
	file_query_hases(query_file);
	file_query_rank(query_file);
}



void kmer_Set_Light::dump_disk(const string& output_file){
	//OPEN
	filebuf fb;
	remove(output_file.c_str());
	fb.open (output_file, ios::out | ios::binary | ios::trunc);
	zstr::ostream out(&fb);

	//VARIOUS INTEGERS
	out.write(reinterpret_cast<const char*>(&k),sizeof(k));
	out.write(reinterpret_cast<const char*>(&m1),sizeof(m1));
	out.write(reinterpret_cast<const char*>(&m3),sizeof(m3));
	out.write(reinterpret_cast<const char*>(&minimizer_size_graph),sizeof(minimizer_size_graph));
	out.write(reinterpret_cast<const char*>(&bit_saved_sub),sizeof(bit_saved_sub));
	uint64_t positions_size(positions.size()), bucketSeq_size(bucketSeq.size());
	out.write(reinterpret_cast<const char*>(&positions_size),sizeof(positions_size));
	out.write(reinterpret_cast<const char*>(&bucketSeq_size),sizeof(bucketSeq_size));

	//BOOL VECTOR
	dump_vector_bool(bucketSeq,&out);
	dump_vector_bool(positions,&out);


	//BUCKETS INFORMATION

	for(uint64_t i(0);i<mphf_number;++i){
	out.write(reinterpret_cast<const char*>(&all_mphf[i].mphf_size),sizeof(all_mphf[i].mphf_size));
	out.write(reinterpret_cast<const char*>(&all_mphf[i].empty),sizeof(all_mphf[i].empty));
	out.write(reinterpret_cast<const char*>(&all_mphf[i].start),sizeof(all_mphf[i].start));
	out.write(reinterpret_cast<const char*>(&all_mphf[i].bit_to_encode),sizeof(all_mphf[i].bit_to_encode));
	if(not all_mphf[i].empty){
			all_mphf[i].kmer_MPHF->save(out);
		}
	}
	for(uint64_t i(0);i<minimizer_number;++i){
		out.write(reinterpret_cast<const char*>(&all_buckets[i].skmer_number),sizeof(all_buckets[i].skmer_number));
	}

	//BM VECTOR
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

	out<<flush;
	fb.close();
}




void kmer_Set_Light::dump_and_destroy(const string& output_file){
	//OPEN
	filebuf fb;
	remove(output_file.c_str());
	fb.open (output_file, ios::out | ios::binary | ios::trunc);
	zstr::ostream out(&fb);

	//VARIOUS INTEGERS
	out.write(reinterpret_cast<const char*>(&k),sizeof(k));
	out.write(reinterpret_cast<const char*>(&m1),sizeof(m1));
	out.write(reinterpret_cast<const char*>(&m3),sizeof(m3));
	out.write(reinterpret_cast<const char*>(&minimizer_size_graph),sizeof(minimizer_size_graph));
	out.write(reinterpret_cast<const char*>(&bit_saved_sub),sizeof(bit_saved_sub));
	uint64_t positions_size(positions.size()), bucketSeq_size(bucketSeq.size());
	out.write(reinterpret_cast<const char*>(&positions_size),sizeof(positions_size));
	out.write(reinterpret_cast<const char*>(&bucketSeq_size),sizeof(bucketSeq_size));

	//BOOL VECTOR
	vector<bool> nadine;
	dump_vector_bool(bucketSeq,&out);
	bucketSeq.clear();
	bucketSeq.swap(nadine);
	dump_vector_bool(positions,&out);
	positions.clear();
	positions.swap(nadine);

	//BCUKETS INFORMATION
	for(uint64_t i(0);i<mphf_number;++i){
		out.write(reinterpret_cast<const char*>(&all_mphf[i].mphf_size),sizeof(all_mphf[i].mphf_size));
		out.write(reinterpret_cast<const char*>(&all_mphf[i].empty),sizeof(all_mphf[i].empty));
		out.write(reinterpret_cast<const char*>(&all_mphf[i].start),sizeof(all_mphf[i].start));
		out.write(reinterpret_cast<const char*>(&all_mphf[i].bit_to_encode),sizeof(all_mphf[i].bit_to_encode));
		if(not all_mphf[i].empty){
			all_mphf[i].kmer_MPHF->save(out);
			delete all_mphf[i].kmer_MPHF;
			all_mphf[i].empty=true;
		}
	}
	for(uint64_t i(0);i<minimizer_number;++i){
		out.write(reinterpret_cast<const char*>(&all_buckets[i].skmer_number),sizeof(all_buckets[i].skmer_number));
	}

	//BM VECTOR
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

	out<<flush;
	fb.close();
}



kmer_Set_Light::kmer_Set_Light(const string& index_file){
	zstr::ifstream out(index_file);

	//VARIOUS INTEGERS
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

	//BOOL VECTOR
	read_vector_bool(bucketSeq,&out,bucketSeq_size);
	read_vector_bool(positions,&out,positions_size);

	//BUCKETS INFORMATION
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
	for(uint64_t i(0);i<minimizer_number;++i){
		out.read(reinterpret_cast< char*>(&all_buckets[i].skmer_number),sizeof(all_buckets[i].skmer_number));
	}

	//BM VECTOR
	uint64_t sz;
	out.read(reinterpret_cast< char*>(&sz),sizeof(sz));
	uint8_t* buff= new uint8_t[sz];
	out.read((char*)buff,sz);
	bm::deserialize(position_super_kmers, buff);
	delete buff;


	position_super_kmers_RS=new bm::bvector<>::rs_index_type();
	position_super_kmers.build_rs_index(position_super_kmers_RS);

	cout<<"Index loaded"<<endl;
}
