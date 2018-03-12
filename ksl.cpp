
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
#include "ksl.h"
#include <omp.h>




using namespace std;


string SEQOFINTERESET("TCAGCGAGGACGGGTATCCGGTTTCCGTCTT");



typedef boomphf::SingleHashFunctor<kmer>  hasher;
typedef boomphf::mphf<  kmer, hasher  > MPHF;




uint64_t nuc2int(char c){
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
	}
	cout<<"bug"<<c<<"!"<<endl;
	exit(0);
	return 0;
}



uint64_t nuc2intrc(char c){
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
	}
	cout<<"bug"<<c<<"!"<<endl;
	exit(0);
	return 0;
}

string intToString(uint64_t n){
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



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str){
	return (min(str,revComp(str)));
}


kmer str2num(const string& str){
	kmer res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			case 'T':res+=3;break;
			default:cout<<"bug"<<"!"<<endl;
	exit(0);
		}
	}
	return res;
}


//~ uint64_t xs(uint64_t y){
	//~ y^=(y<<13); y^=(y>>17);y=(y^=(y<<15)); return y;
//~ }
uint64_t xs(uint64_t y){
	return y;
}


string min_accordingtoXS(const string& seq1,const string& seq2){
	uint32_t u1(str2num(seq1)),u2(str2num(seq2));
	if(xs(u1)<xs(u2)){
		return seq1;
	}
	return seq2;
}



kmer min_accordingtoXS(const kmer& u1,const kmer& u2){
	if(xs(u1)<xs(u2)){
		return u1;
	}
	return u2;
}



kmer rcb(kmer min,uint n){
	kmer res(0);
	kmer offset(1);
	offset<<=(2*n-2);
	for(uint i(0); i<n;++i){
		res+=(3-(min%4))*offset;
		min>>=2;
		offset>>=2;
	}
	return res;
}

void kmer_Set_Light::updateK(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchor;
}

void kmer_Set_Light::updateM(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=(1<<(2*m1));
}



kmer min_k (const kmer& k1,const kmer& k2){
	if(k1<=k2){
		return k1;
	}
	return k2;
}



void kmer_Set_Light::updateRCK(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-2));
}

void kmer_Set_Light::updateRCM(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*m1-2));
}



uint32_t knuth_hash (uint32_t x){
	return x*2654435761;
}



size_t hash2(int i1)
{
    size_t ret = i1;
    ret *= 2654435761U;
    return ret ^ 69;
}




uint32_t kmer_Set_Light::minimizer_according_xs(kmer seq){
	uint32_t mini,mmer;
	mini=seq%minimizer_number;
	mini=min(mini,(uint32_t)rcb(mini,m1));
	for(uint i(0);i<k-m1;i++){
		seq>>=2;
		mmer=seq%minimizer_number;
		mmer=min(mmer,(uint32_t)rcb(mmer,m1));
		if(abundance_minimizer[mini]>abundance_minimizer[mmer]){
			mini=mmer;
		}else{
			if(abundance_minimizer[mini]==abundance_minimizer[mmer]){
				mini=((xs(mini)<=xs(mmer))?mini:mmer);
			}
		}
		//~ mini=((xs(mini)<=xs(mmer))?mini:mmer);
	}
	return mini;
}



void kmer_Set_Light::abundance_minimizer_construct(const string& input_file){
	ifstream inUnitigs(input_file);
	if( not inUnitigs.good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	string ref,useless;
	kmer super_minimizer,minimizer;
	while(not inUnitigs.eof()){
		getline(inUnitigs,useless);
		getline(inUnitigs,ref);
		//FOREACH UNITIG
		if(not ref.empty() and not useless.empty()){
			uint last_position(0);
			//FOREACH KMER
			kmer seq(str2num(ref.substr(0,m1))),rcSeq(rcb(seq,m1)),canon(min_k(seq,rcSeq));
			abundance_minimizer[canon]++;
			uint i(0);
			for(;i+m1<ref.size();++i){
				updateM(seq,ref[i+m1]);
				updateRCM(rcSeq,ref[i+m1]);
				canon=(min_k(seq,rcSeq));
				abundance_minimizer[canon]++;
			}
		}
	}
}



void kmer_Set_Light::create_super_buckets(const string& input_file){
	ifstream inUnitigs(input_file);
	if( not inUnitigs.good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	vector<ofstream*> out_files;
	for(uint i(0);i<number_superbuckets;++i){
		auto out =new ofstream("_out"+to_string(i));
		out_files.push_back(out);
	}
	omp_lock_t lock[number_superbuckets];
	for (int i=0; i<number_superbuckets; i++)
        omp_init_lock(&(lock[i]));

	#pragma omp parallel num_threads(coreNumber)
	{
		string ref,useless;
		minimizer_type super_minimizer,minimizer;
		while(not inUnitigs.eof()){
			#pragma omp critical(dataupdate)
			{
				getline(inUnitigs,useless);
				getline(inUnitigs,ref);
			}
			//FOREACH UNITIG
			if(not ref.empty() and not useless.empty()){
				super_minimizer=minimizer=minimizer_number;
				uint last_position(0);
				//FOREACH KMER
				kmer seq(str2num(ref.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				minimizer=minimizer_according_xs(canon);
				super_minimizer=minimizer;
				uint i(0);
				for(;i+k<ref.size();++i){
					updateK(seq,ref[i+k]);
					updateRCK(rcSeq,ref[i+k]);
					canon=(min_k(seq, rcSeq));
					//COMPUTE KMER MINIMIZER
					minimizer=minimizer_according_xs(canon);
					if(super_minimizer!=minimizer){
						omp_set_lock(&(lock[((super_minimizer))/bucket_per_superBuckets]));
						*(out_files[((super_minimizer))/bucket_per_superBuckets])<<">"+to_string(super_minimizer)+"\n"<<ref.substr(last_position,i-last_position+k)<<"\n";
						last_position=i+1;
						omp_unset_lock(&(lock[((super_minimizer))/bucket_per_superBuckets]));
						super_minimizer=minimizer;
					}
				}
				if(ref.size()-last_position>k-1){
					omp_set_lock(&(lock[((super_minimizer))/bucket_per_superBuckets]));
					*(out_files[((super_minimizer))/bucket_per_superBuckets])<<">"+to_string(super_minimizer)+"\n"<<ref.substr(last_position)<<"\n";
					omp_unset_lock(&(lock[((super_minimizer))/bucket_per_superBuckets]));
				}
			}
		}
	}
	for(uint i(0);i<number_superbuckets;++i){
		out_files[i]->close();
		delete(out_files[i]);
	}
}


void kmer_Set_Light::str2bool(const string& str,uint mini){
	for(uint i(0);i<str.size();++i){
		if(not Valid_kmer.empty()){
			Valid_kmer[mini].push_back(true);
		}
		switch (str[i]){
			case 'A':
				bucketSeq[mini].push_back(false);
				bucketSeq[mini].push_back(false);
				break;
			case 'C':
				bucketSeq[mini].push_back(false);
				bucketSeq[mini].push_back(true);
				break;
			case 'G':
				bucketSeq[mini].push_back(true);
				bucketSeq[mini].push_back(false);
				break;
			case 'T':
				bucketSeq[mini].push_back(true);
				bucketSeq[mini].push_back(true);
				break;
			default:
				cout<<"nope"<<endl;
			}
	}
	if(not Valid_kmer.empty()){
		for(uint i(0);i<k-1;++i){
			Valid_kmer[mini][Valid_kmer[mini].size()-k+i+1]=(false);
		}
	}
}


void kmer_Set_Light::read_super_buckets(const string& input_file){
	uint64_t total_size(0);

	uint number_superbuckets_by_buckets(minimizer_number/number_superbuckets);
	#pragma omp parallel num_threads(coreNumber)
	{
		string useless,line;
		#pragma omp for
		for(uint SBC=0;SBC<number_superbuckets;++SBC){
			uint BC(SBC*bucket_per_superBuckets);
			ifstream in(input_file+to_string(SBC));
			while(not in.eof()){
				getline(in,useless);
				getline(in,line);
				if(not useless.empty()){
					useless=useless.substr(1);
					uint minimizer(stoi(useless));
					str2bool(line,minimizer);
					#pragma omp atomic
					number_kmer+=line.size()-k+1;
					#pragma omp atomic
					number_super_kmer++;
					#pragma omp atomic
					total_size+=line.size();
				}
			}
			remove((input_file+to_string(SBC)).c_str());
			//~ cout<<"gomphf"<<endl;
			create_mphf(BC,BC+bucket_per_superBuckets);
			//~ cout<<"gofill"<<endl;
			fill_positions(BC,BC+bucket_per_superBuckets);

			BC+=bucket_per_superBuckets;
			cout<<"-"<<flush;
		}
	}
	cout<<endl;
	cout<<"----------------------INDEX RECAP----------------------------"<<endl;
	cout<<"Kmer in graph: "<<intToString(number_kmer)<<endl;
	cout<<"Super Kmer in graph: "<<intToString(number_super_kmer)<<endl;
	cout<<"Average size of Super Kmer: "<<intToString(total_size/(number_super_kmer))<<endl;
	cout<<"Total size of the partitionned graph: "<<intToString(total_size)<<endl;
	cout<<"Largest MPHF: "<<intToString(largest_MPHF)<<endl;
	cout<<"Largest Bucket: "<<intToString(largest_bucket_nuc_all)<<endl;

	cout<<"Size of the partitionned graph (MBytes): "<<intToString(total_size/(4*1024*1024))<<endl;
	cout<<"Total Positions size (MBytes): "<<intToString(positions_total_size/(8*1024*1024))<<endl;
	if(not light_mode){
		cout<<"Space used for separators (MBytes): "<<intToString(total_size/(8*1024*1024))<<endl;
	}
	cout<<"Size of the partitionned graph (bit per kmer): "<<((double)(2*total_size)/(number_kmer))<<endl;
	bit_per_kmer+=((double)(2*total_size)/(number_kmer));
	if(not light_mode){
		cout<<"Space used for separators (bit per kmer): "<<((double)total_size/(number_kmer))<<endl;
		bit_per_kmer+=((double)total_size/(number_kmer));
	}
	cout<<"Total Positions size (bit per kmer): "<<((double)positions_total_size/number_kmer)<<endl;
	bit_per_kmer+=((double)positions_total_size/number_kmer);

	cout<<"TOTAL Bits per kmer (without bbhash): "<<bit_per_kmer<<endl;
}


kmer kmer_Set_Light::get_kmer(uint32_t mini,uint32_t pos){
	kmer res(0);
	for(uint i(0);i<k;++i){
		res<<=2;
		if(bucketSeq[mini][(pos+i)*2]){
			if(bucketSeq[mini][(pos+i)*2+1]){
				res+=3;
			}else{
				res+=2;
			}
		}else{
			if(bucketSeq[mini][(pos+i)*2+1]){
				res+=1;
			}else{
			}
		}
	}
	return res;
}


kmer kmer_Set_Light::update_kmer(uint32_t pos,uint32_t mini,kmer input){
	input<<=2;
	if(bucketSeq[mini][(pos)*2]){
		if(bucketSeq[mini][(pos)*2+1]){
			input+=3;
		}else{
			input+=2;
		}
	}else{
		if(bucketSeq[mini][(pos)*2+1]){
			input+=1;
		}else{
		}
	}
	return input%offsetUpdateAnchor;
}

void kmer_Set_Light::print_kmer(kmer num){
	kmer anc=1;
	anc<<=2*(k-1);
	for(uint i(0);i<k;++i){
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




void kmer_Set_Light::create_mphf(uint begin_BC,uint end_BC){
	#pragma omp parallel  num_threads(1)
	{
		uint64_t anchors_number(0);
		auto anchors=new vector<kmer>;
		uint largest_bucket_anchor(0);
		uint64_t largest_bucket_nuc(0);
		//~ #pragma omp for schedule(dynamic, 16*number_bucket_per_mphf)
		for(uint BC=(begin_BC);BC<end_BC;++BC){
			if(bucketSeq[BC].size()!=0){
				largest_bucket_nuc=max(largest_bucket_nuc,bucketSeq[BC].size());
				largest_bucket_nuc_all=max(largest_bucket_nuc_all,bucketSeq[BC].size());
				uint bucketSize(1);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				anchors->push_back(canon);
				for(uint j(0);2*(j+k)<bucketSeq[BC].size();j++){
					if(not Valid_kmer.empty() and not Valid_kmer[BC][j+1]){
						j+=k-1;
						if(2*(j+k)<bucketSeq[BC].size()){
							seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
							anchors->push_back(canon);
							bucketSize++;
						}
					}else{
						seq=update_kmer(j+k,BC,seq);
						rcSeq=(rcb(seq,k));
						canon=(min_k(seq, rcSeq));
						anchors->push_back(canon);
						bucketSize++;
					}
				}
				largest_bucket_anchor=max(largest_bucket_anchor,bucketSize);
			}
			if((BC+1)%number_bucket_per_mphf==0){
				//~ cout<<BC/number_bucket_per_mphf<<endl;
				//~ cout<<BC<<endl;
				//~ cout<<number_bucket_per_mphf<<endl;

				largest_MPHF=max(largest_MPHF,anchors->size());
				anchors_number=anchors->size();
				auto data_iterator3 = boomphf::range(static_cast<const kmer*>(&(*anchors)[0]), static_cast<const kmer*>((&(*anchors)[0])+anchors->size()));
				kmer_MPHF[BC/number_bucket_per_mphf]= boomphf::mphf<kmer,hasher>(anchors->size(),data_iterator3,coreNumber,gammaFactor,false);
				mphf_size[BC/number_bucket_per_mphf]=largest_bucket_anchor;
				anchors->clear();
				int n_bits_to_encode((ceil(log2(largest_bucket_nuc/2+1))-bit_saved_sub));
				if(n_bits_to_encode<1){n_bits_to_encode=1;}
				bit_to_encode[BC/number_bucket_per_mphf]=n_bits_to_encode;
				positions[BC/number_bucket_per_mphf].resize((anchors_number*n_bits_to_encode),0);
				//~ cout<<anchors_number<<endl;
				#pragma omp atomic
				positions_total_size+=(anchors_number*n_bits_to_encode);
				largest_bucket_anchor=0;
				largest_bucket_nuc=(0);
			}
		}
		delete anchors;
	}
}



void int_to_bool(uint n_bits_to_encode,uint32_t X, uint32_t pos,vector<bool>& res){
	//~ cout<<"itb"<<endl;
	for(uint i(0);i<n_bits_to_encode;++i){
		res[i+pos*n_bits_to_encode]=X%2;
		X>>=1;
	}
	//~ cout<<"itbE"<<endl;
}


uint32_t bool_to_int(uint n_bits_to_encode,uint pos,const vector<bool>& V){
	uint32_t res(0);
	uint32_t acc(1);
	for(uint i(0);i<n_bits_to_encode;++i){
		if(V[i+pos*n_bits_to_encode]){
			res+=acc;
		}else{
		}
		acc<<=1;
	}
		return res;
	}



void kmer_Set_Light::fill_positions(uint begin_BC,uint end_BC){
	//~ uint64_t total_size(0);
	#pragma omp parallel for num_threads(1)
	for(uint BC=(begin_BC);BC<end_BC;++BC){
		if(bucketSeq[BC].size()==0){
			continue;
		}
		//~ cout<<bit_to_encode.size()<<" "<<BC/number_bucket_per_mphf<<endl;
		int n_bits_to_encode(bit_to_encode[BC/number_bucket_per_mphf]);
		//~ cout<<n_bits_to_encode<<endl;
		kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
		for(uint j(0);2*(j+k)<bucketSeq[BC].size();j++){
			if(not Valid_kmer.empty() and not Valid_kmer[BC][j+1]){
				j+=k-1;
				if(2*(j+k)<bucketSeq[BC].size()){
					seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
					int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,kmer_MPHF[BC/number_bucket_per_mphf].lookup(canon),positions[BC/number_bucket_per_mphf]);
				}
			}else{
				seq=update_kmer(j+k,BC,seq);
				rcSeq=(rcb(seq,k));
				canon=(min_k(seq, rcSeq));
				int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,kmer_MPHF[BC/number_bucket_per_mphf].lookup(canon),positions[BC/number_bucket_per_mphf]);
			}
		}
		if(light_mode){
			Valid_kmer[BC].clear();
			vector<bool>().swap(Valid_kmer[BC]);
		}
	}
}



int64_t kmer_Set_Light::correct_pos(uint32_t mini, uint64_t p){
	if(Valid_kmer.empty() or Valid_kmer[mini].size()<p+k){
		return p;
	}
	for(uint i(0);i<k;i++){
		if(Valid_kmer[mini][p+i]){
			return (p+i);
		}else{
		}
	}
	return p;
}



void kmer_Set_Light::multiple_query(const string& query_file){
	if(light_mode){
		Valid_kmer.clear();
	}
	ifstream in(query_file);
	uint TP(0),FP(0);
	#pragma omp parallel num_threads(coreNumber)
	while(not in.eof()){
		string query,Q;
		#pragma omp critical(dataupdate)
		{
			getline(in,query);
			getline(in,query);
		}
		//~ cout<<"go"<<endl;
		if(query.size()>=k){
			//~ cout<<"QUERY"<<endl;
			kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq)),canonR,seqR,rcSeqR;
			uint i(0);
			canon=(min_k(seq, rcSeq));
			//~ cout<<0<<endl;
			//COMPUTE KMER MINIMIZER
			uint32_t minimizer=(minimizer_according_xs(canon));
			int32_t hash;
			if(mphf_size[minimizer/number_bucket_per_mphf]==0){
				hash=-1;
			}else{
				hash=(kmer_MPHF[minimizer/number_bucket_per_mphf].lookup(canon));
			}
			//~ cout<<02<<endl;
			if(hash<0){
				#pragma omp atomic
				FP++;
				//~ cout<<1<<endl;
					//~ cin.get();
			}else{
				//~ cout<<2<<endl;
				int n_bits_to_encode(bit_to_encode[minimizer/number_bucket_per_mphf]);
				//~ cout<<hash<<endl;
				//~ cout<<positions[minimizer/number_bucket_per_mphf].size()<<endl;
				uint pos(bool_to_int( n_bits_to_encode, hash, positions[minimizer/number_bucket_per_mphf]));
				//~ cout<<20<<endl;
				if(2*(pos+k-1)<bucketSeq[minimizer].size()){
					pos=correct_pos(minimizer,(positions_to_check)*pos);
					//~ cout<<21<<endl;
					seqR=get_kmer(minimizer,pos);
					rcSeqR=rcb(seqR,k);
					canonR=(min_k(seqR, rcSeqR));
					//~ cout<<22<<endl;
					if(canon==canonR){
						//~ cout<<pos<<endl;
						#pragma omp atomic
						TP++;
					}else{
						//~ cout<<3<<endl;
						uint j;
						bool found(false);
						for(j=(pos);j<pos+positions_to_check;++j){
							if(not light_mode and not Valid_kmer[minimizer][j+1]){
								j+=k-1;
								if(2*(j+k)<bucketSeq[minimizer].size()){
									seqR=(get_kmer(minimizer,j+1)),rcSeqR=(rcb(seqR,k)),canonR=(min_k(seqR,rcSeqR));
									//~ canonR=(min_k(seqR, rcSeqR));
								}
							}else{
								seqR=update_kmer(j+k,minimizer,seqR);
								rcSeqR=(rcb(seqR,k));
								canonR=(min_k(seqR, rcSeqR));
							}

							if(canon==canonR){
								#pragma omp atomic
								TP++;
								found=true;
								break;
							}
						}
						if(not found){
							#pragma omp atomic
							FP++;
							//~ cout<<2;;
							//~ print_kmer(seq);
							//~ print_kmer(seqR);
						}
					}
				}else{
					#pragma omp atomic
					FP++;
				}
			}
			//~ cout<<11<<endl;
			for(;i+k<query.size();++i){
				//~ cout<<4<<endl;
				updateK(seq,query[i+k]);
				updateRCK(rcSeq,query[i+k]);
				canon=(min_k(seq, rcSeq));
				//COMPUTE KMER MINIMIZER
				uint32_t minimizer=(minimizer_according_xs(canon));
				if(mphf_size[minimizer/number_bucket_per_mphf]==0){
					hash=-1;
				}else{
					hash=(kmer_MPHF[minimizer/number_bucket_per_mphf].lookup(canon));
				}
				if(hash<0){
					#pragma omp atomic
					FP++;
					//~ cout<<3<<endl;
					//~ cin.get();
					continue;
				}
				int n_bits_to_encode(bit_to_encode[minimizer/number_bucket_per_mphf]);
				uint pos(bool_to_int( n_bits_to_encode, hash, positions[minimizer/number_bucket_per_mphf]));
				if(2*(pos+k-1)>=bucketSeq[minimizer].size()){
					#pragma omp atomic
					FP++;
					continue;
				}
				//~ cout<<5<<endl;
				pos=correct_pos(minimizer,(positions_to_check)*pos);
				seqR=(get_kmer(minimizer,pos));
				rcSeqR=rcb(seqR,k);
				canonR=(min_k(seqR, rcSeqR));
				if(canon==canonR){
					#pragma omp atomic
					TP++;
					continue;
				}
				uint j;
				bool found(false);
				for(j=(pos);j<pos+positions_to_check;++j){
					if(not light_mode and not Valid_kmer[minimizer][j+1]){
						j+=k-1;
						if(2*(j+k)<bucketSeq[minimizer].size()){
							seqR=(seqR=(get_kmer(minimizer,j+1))),rcSeqR=(rcb(seqR,k)),canonR=(min_k(seqR,rcSeqR));
						}
					}else{
						seqR=update_kmer(j+k,minimizer,seqR);
						rcSeqR=(rcb(seqR,k));
						canonR=(min_k(seqR, rcSeqR));
					}
					if(canon==canonR){
						#pragma omp atomic
						TP++;
						found=true;
						break;
					}
				}
				if(not found){
					#pragma omp atomic
					FP++;

				}
			}
			//~ cout<<12<<endl;
			//~ cout<<query<<endl;
		}
	}
	cout<<"-----------------------QUERY RECAP----------------------------"<<endl;
	cout<<"Good kmer: "<<intToString(TP)<<endl;
	cout<<"Erroneous kmers: "<<intToString(FP)<<endl;
}










