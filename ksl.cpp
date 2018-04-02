
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
		if(all_buckets[mini]==NULL or all_buckets[mmer]==NULL){
			return (3);
		}
		uint value_min(all_buckets[mini]->abundance_minimizer),value_mmer(all_buckets[mmer]->abundance_minimizer);

		if(value_min>value_mmer){
			mini=mmer;
		}else{
			if(value_min==value_mmer){
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
			//~ abundance_minimizer[canon]++;
			if(all_buckets[canon]==NULL){
				all_buckets[canon]=new bucket_minimizer;
				all_buckets[canon]->abundance_minimizer=0;
				all_buckets[canon]->nuc_minimizer=0;
				all_buckets[canon]->current_pos=0;
			}else{
				all_buckets[canon]->abundance_minimizer++;
			}
			uint i(0);
			for(;i+m1<ref.size();++i){
				updateM(seq,ref[i+m1]);
				updateRCM(rcSeq,ref[i+m1]);
				canon=(min_k(seq,rcSeq));
				if(all_buckets[canon]==NULL){
					all_buckets[canon]=new bucket_minimizer;
					all_buckets[canon]->abundance_minimizer=0;
					all_buckets[canon]->nuc_minimizer=0;
					all_buckets[canon]->current_pos=0;
				}else{
					all_buckets[canon]->abundance_minimizer++;
				}
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
	for (int i=0; i<number_superbuckets; i++){
        omp_init_lock(&(lock[i]));
	}

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

						omp_unset_lock(&(lock[((super_minimizer))/bucket_per_superBuckets]));
						all_buckets[super_minimizer]->nuc_minimizer+=(i-last_position+k);
						last_position=i+1;
						super_minimizer=minimizer;
					}
				}
				if(ref.size()-last_position>k-1){
					omp_set_lock(&(lock[((super_minimizer))/bucket_per_superBuckets]));
					*(out_files[((super_minimizer))/bucket_per_superBuckets])<<">"+to_string(super_minimizer)+"\n"<<ref.substr(last_position)<<"\n";
					omp_unset_lock(&(lock[((super_minimizer))/bucket_per_superBuckets]));
					all_buckets[super_minimizer]->nuc_minimizer+=(ref.substr(last_position)).size();
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
		Valid_kmer[mini%bucket_per_superBuckets].push_back(true);
		switch (str[i]){
			case 'A':
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2]=(false);
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2+1]=(false);
				break;
			case 'C':
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2]=(false);
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2+1]=(true);
				break;
			case 'G':
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2]=(true);
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2+1]=(false);
				break;
			case 'T':
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2]=(true);
				all_buckets[mini]->bucketSeq[all_buckets[mini]->current_pos+i*2+1]=(true);
				break;
			default:
				cout<<"nope"<<endl;
			}
	}
	all_buckets[mini]->current_pos+=(str.size()*2);
	for(uint i(0);i<k-1;++i){
		Valid_kmer[mini%bucket_per_superBuckets][Valid_kmer[mini%bucket_per_superBuckets].size()-k+i+1]=(false);
	}
}


void kmer_Set_Light::read_super_buckets(const string& input_file){
	uint64_t total_size(0);

	uint number_superbuckets_by_buckets(minimizer_number/number_superbuckets);
	#pragma omp parallel num_threads(1)
	{
		string useless,line;
		#pragma omp for
		for(uint SBC=0;SBC<number_superbuckets;++SBC){
			uint BC(SBC*bucket_per_superBuckets);
			for(uint ii(BC);ii<BC+bucket_per_superBuckets;++ii){
				if(all_buckets[ii]!=NULL){
					all_buckets[ii]->bucketSeq.resize(all_buckets[ii]->nuc_minimizer*2);
					all_buckets[ii]->bucketSeq.shrink_to_fit();
					total_size+=all_buckets[ii]->bucketSeq.capacity()/2;
				}
			}
			ifstream in(input_file+to_string(SBC));
			while(not in.eof()){
				useless="";
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
					//~ #pragma omp atomic
					//~ total_size+=line.size();
				}
			}
			remove((input_file+to_string(SBC)).c_str());
			create_mphf(BC,BC+bucket_per_superBuckets);
			fill_positions(BC,BC+bucket_per_superBuckets);
			BC+=bucket_per_superBuckets;
			cout<<"-"<<flush;
		}
	}
	//~ for(uint ii(0);ii<bucket_per_superBuckets;++ii){
		//~ Valid_kmer[ii].clear();
	//~ }
	delete[] Valid_kmer;
	cout<<endl;
	cout<<"----------------------INDEX RECAP----------------------------"<<endl;
	cout<<"Kmer in graph: "<<intToString(number_kmer)<<endl;
	cout<<"Super Kmer in graph: "<<intToString(number_super_kmer)<<endl;
	cout<<"Average size of Super Kmer: "<<intToString(number_kmer/(number_super_kmer))<<endl;
	cout<<"Total size of the partitionned graph: "<<intToString(total_size)<<endl;
	cout<<"Largest MPHF: "<<intToString(largest_MPHF)<<endl;
	cout<<"Largest Bucket: "<<intToString(largest_bucket_nuc_all)<<endl;

	cout<<"Size of the partitionned graph (MBytes): "<<intToString(total_size/(4*1024*1024))<<endl;
	if(not light_mode){
		cout<<"Space used for separators (MBytes): "<<intToString(total_size/(8*1024*1024))<<endl;
	}
	cout<<"Total Positions size (MBytes): "<<intToString(positions_total_size/(8*1024*1024))<<endl;
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
		if(all_buckets[mini]->bucketSeq[(pos+i)*2]){
			if(all_buckets[mini]->bucketSeq[(pos+i)*2+1]){
				res+=3;
			}else{
				res+=2;
			}
		}else{
			if(all_buckets[mini]->bucketSeq[(pos+i)*2+1]){
				res+=1;
			}else{
			}
		}
	}
	return res;
}


kmer kmer_Set_Light::update_kmer(uint32_t pos,uint32_t mini,kmer input){
	input<<=2;
	if(all_buckets[mini]->bucketSeq[(pos)*2]){
		if(all_buckets[mini]->bucketSeq[(pos)*2+1]){
			input+=3;
		}else{
			input+=2;
		}
	}else{
		if(all_buckets[mini]->bucketSeq[(pos)*2+1]){
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
		#pragma omp for schedule(dynamic, 1)
		for(uint BC=(begin_BC);BC<end_BC;++BC){

			if(all_buckets[BC]!=NULL){
				if(all_buckets[BC]->bucketSeq.size()!=0){
					largest_bucket_nuc=max(largest_bucket_nuc,all_buckets[BC]->bucketSeq.size());
					largest_bucket_nuc_all=max(largest_bucket_nuc_all,all_buckets[BC]->bucketSeq.size());
					uint bucketSize(1);
					kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
					anchors->push_back(canon);
					for(uint j(0);2*(j+k)<all_buckets[BC]->bucketSeq.size();j++){
						if(not Valid_kmer[BC%bucket_per_superBuckets][j+1]){
						//~ if(false){
							j+=k-1;
							if(2*(j+k)<all_buckets[BC]->bucketSeq.size()){
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
			}
			if((BC+1)%number_bucket_per_mphf==0){
				largest_MPHF=max(largest_MPHF,anchors->size());
				anchors_number=anchors->size();
				auto data_iterator3 = boomphf::range(static_cast<const kmer*>(&(*anchors)[0]), static_cast<const kmer*>((&(*anchors)[0])+anchors->size()));
				all_mphf[BC/number_bucket_per_mphf]=new info_mphf;
				all_mphf[BC/number_bucket_per_mphf]->kmer_MPHF= boomphf::mphf<kmer,hasher>(anchors->size(),data_iterator3,coreNumber,gammaFactor,false);
				all_mphf[BC/number_bucket_per_mphf]->mphf_size=largest_bucket_anchor;
				anchors->clear();
				int n_bits_to_encode((ceil(log2(largest_bucket_nuc/2+1))-bit_saved_sub));
				if(n_bits_to_encode<1){n_bits_to_encode=1;}
				all_mphf[BC/number_bucket_per_mphf]->bit_to_encode=n_bits_to_encode;
				all_mphf[BC/number_bucket_per_mphf]->positions.resize((anchors_number*n_bits_to_encode),0);
				all_mphf[BC/number_bucket_per_mphf]->positions.shrink_to_fit();
				#pragma omp atomic
				//~ positions_total_size+=(anchors_number*n_bits_to_encode);
				positions_total_size+=all_mphf[BC/number_bucket_per_mphf]->positions.capacity();
				largest_bucket_anchor=0;
				largest_bucket_nuc=(0);
			}
		}
		delete anchors;
	}
}



void int_to_bool(uint n_bits_to_encode,uint32_t X, uint32_t pos,vector<bool>& res){
	for(uint i(0);i<n_bits_to_encode;++i){
		res[i+pos*n_bits_to_encode]=X%2;
		X>>=1;
	}
}



uint32_t kmer_Set_Light::bool_to_int(uint n_bits_to_encode,uint pos,const vector<bool>& V){
	uint32_t res(0);
	uint32_t acc(1);
	for(uint i(0);i<n_bits_to_encode;++i){
		if(V[i+pos*n_bits_to_encode]){
			res+=acc;
		}else{
		}
		acc<<=1;
	}
	return res*positions_to_check;
}



void kmer_Set_Light::fill_positions(uint begin_BC,uint end_BC){
	//~ uint64_t total_size(0);
	#pragma omp parallel for num_threads(coreNumber)
	for(uint BC=(begin_BC);BC<end_BC;++BC){
		if(not all_buckets[BC]==NULL){
			if(all_buckets[BC]->bucketSeq.size()>0){
				int n_bits_to_encode(all_mphf[BC/number_bucket_per_mphf]->bit_to_encode);
				kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				for(uint j(0);2*(j+k)<all_buckets[BC]->bucketSeq.size();j++){
					if(not Valid_kmer[BC%bucket_per_superBuckets][j+1]){
						j+=k-1;
						if(2*(j+k)<all_buckets[BC]->bucketSeq.size()){
							seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));

							#pragma omp critical(dataupdate)
							{
								int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,all_mphf[BC/number_bucket_per_mphf]->kmer_MPHF.lookup(canon),all_mphf[BC/number_bucket_per_mphf]->positions);
							}
						}
					}else{
						seq=update_kmer(j+k,BC,seq);
						rcSeq=(rcb(seq,k));
						canon=(min_k(seq, rcSeq));
						#pragma omp critical(dataupdate)
						{
							int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,all_mphf[BC/number_bucket_per_mphf]->kmer_MPHF.lookup(canon),all_mphf[BC/number_bucket_per_mphf]->positions);
						}
					}
				}
			}
		}
		if(light_mode){
			//~ all_buckets[BC]->valid_kmers.clear();
			//~ vector<bool>().swap(all_buckets[BC]->valid_kmers);

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
		}else{
		}
	}
	return p;
}




















void kmer_Set_Light::get_anchors(const string& query,vector<uint>& minimizerV, vector<kmer>& kmerV){
	kmerV.clear();
	minimizerV.clear();
	kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq)),canonR,seqR,rcSeqR;
	uint i(0);
	canon=(min_k(seq, rcSeq));
	uint32_t minimizer=(minimizer_according_xs(canon));
	kmerV.push_back(canon);
	minimizerV.push_back(minimizer);
	for(;i+k<query.size();++i){
		updateK(seq,query[i+k]);
		updateRCK(rcSeq,query[i+k]);
		canon=(min_k(seq, rcSeq));
		minimizer=(minimizer_according_xs(canon));
		kmerV.push_back(canon);
		minimizerV.push_back(minimizer);
	}
}



uint kmer_Set_Light::multiple_query_serial(const vector<uint>& minimizerV, const vector<kmer>& kmerV){
	uint res(0);
	for(uint i(0);i<kmerV.size();++i){
		int32_t pos=query_get_pos_unitig(kmerV[i],minimizerV[i]);
		if(pos>=0){
			++res;
		}
	}
	return res;
}



uint next_different_value(const vector<uint>& minimizerV,uint start, uint m){
	uint i(0);
	for(;i+start<minimizerV.size();++i){
		if(minimizerV[i+start]!=m){
			return start+i-1;
		}
	}
	return start+i-1;
}



uint kmer_Set_Light::multiple_query_optimized(const vector<uint>& minimizerV, const vector<kmer>& kmerV){
	uint res(0);
	for(uint i(0);i<kmerV.size();++i){
		int32_t pos=query_get_pos_unitig(kmerV[i],minimizerV[i]);
		uint next(next_different_value(minimizerV,i,minimizerV[i]));
		if(next!=i){
			int32_t pos2=query_get_pos_unitig(kmerV[next],minimizerV[next]);
			if(pos2-pos==next-i){
				res+=next-i+1;
				i=next;
			}else{
				if(pos>=0){
					++res;
				}
			}
		}else{
			if(pos>=0){
				++res;
			}
		}
	}
	return res;
}



int32_t kmer_Set_Light::query_get_pos_unitig(const kmer canon,uint minimizer){
	int64_t hash(-1);
	if(all_mphf[minimizer/number_bucket_per_mphf]==NULL){
		return -1;
	}else{
		if(all_mphf[minimizer/number_bucket_per_mphf]->mphf_size==0){
			return -1;
		}else{
			hash=(all_mphf[minimizer/number_bucket_per_mphf]->kmer_MPHF.lookup(canon));
		}
	}
	if(hash<0){
		//~ cout<<"FAILHASH"<<endl;
		return -1;
	}else{
		int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf]->bit_to_encode);
		uint pos(bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf]->positions));
		if(2*(pos+k-1)<all_buckets[minimizer]->bucketSeq.size()){
			kmer seqR=get_kmer(minimizer,pos);
			kmer rcSeqR=rcb(seqR,k);
			kmer canonR=(min_k(seqR, rcSeqR));
			if(canon==canonR){
				return pos;
			}else{
				uint j;
				bool found(false);
				for(j=(pos);j<pos+positions_to_check;++j){
					if(not light_mode and not Valid_kmer[minimizer/number_superbuckets][j+1]){
						j+=k-1;
						if(2*(j+k)<all_buckets[minimizer]->bucketSeq.size()){
							seqR=(get_kmer(minimizer,j+1)),rcSeqR=(rcb(seqR,k)),canonR=(min_k(seqR,rcSeqR));
							canonR=(min_k(seqR, rcSeqR));
						}
					}else{
						seqR=update_kmer(j+k,minimizer,seqR);
						rcSeqR=(rcb(seqR,k));
						canonR=(min_k(seqR, rcSeqR));
					}

					if(canon==canonR){
						return j+1;
					}
				}
			}
		}
	}
	//~ cout<<"FAILPOS"<<endl;
	return -1;
}



void kmer_Set_Light::file_query(const string& query_file){
	ifstream in(query_file);
	uint TP(0),FP(0);
	//~ cout<<"go"<<endl;
	#pragma omp parallel num_threads(coreNumber)
	while(not in.eof()){
		string query;
		vector<kmer> kmerV;
		vector<uint> minimizerV;
		#pragma omp critical(dataupdate)
		{
			getline(in,query);
			getline(in,query);
		}
		if(query.size()>=k){
			get_anchors(query,minimizerV,kmerV);
			//~ uint32_t found=multiple_query_serial(minimizerV,kmerV);
			uint32_t found=multiple_query_optimized(minimizerV,kmerV);
			#pragma omp critical(dataupdate)
			{
				//~ #pragma atomic
				TP+=found;
				//~ #pragma atomic
				FP+=minimizerV.size()-found;
			}
			//~ cout<<minimizerV.size()<<" "<<found<<endl;
		}else{
			//~ cout<<"?"<<endl;
		}
	}
	cout<<"-----------------------QUERY RECAP 2----------------------------"<<endl;
	cout<<"Good kmer: "<<intToString(TP)<<endl;
	cout<<"Erroneous kmers: "<<intToString(FP)<<endl;
}




void kmer_Set_Light::multiple_query(const string& query_file){
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
			if(all_mphf[minimizer/number_bucket_per_mphf]==NULL){
				hash=-1;
			}else{
				if(all_mphf[minimizer/number_bucket_per_mphf]->mphf_size==0){
					hash=-1;
				}else{
					hash=(all_mphf[minimizer/number_bucket_per_mphf]->kmer_MPHF.lookup(canon));
				}
			}
			//~ cout<<02<<endl;
			if(hash<0){
				#pragma omp atomic
				FP++;
				//~ cout<<1<<endl;
				//~ cout<<"MEGA LOL"<<endl;
				//~ print_kmer(canon);
					//~ cin.get();
			}else{
				//~ cout<<2<<endl;
				int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf]->bit_to_encode);
				//~ cout<<hash<<endl;
				//~ cout<<positions[minimizer/number_bucket_per_mphf].size()<<endl;
				uint pos(bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf]->positions));
				//~ cout<<20<<endl;
				//~ cout<<pos<<endl;
				//~ if(all_buckets[minimizer]==NULL){
					//~ cout<<"lol"<<endl;
				//~ }
				if(2*(pos+k-1)<all_buckets[minimizer]->bucketSeq.size()){
					//~ cout<<201<<endl;
					//~ pos=correct_pos(minimizer,(positions_to_check)*pos);
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
							if(not light_mode and not Valid_kmer[minimizer/number_superbuckets][j+1]){
								j+=k-1;
								if(2*(j+k)<all_buckets[minimizer]->bucketSeq.size()){
									seqR=(get_kmer(minimizer,j+1)),rcSeqR=(rcb(seqR,k)),canonR=(min_k(seqR,rcSeqR));
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
							//~ print_kmer(canon);
							//~ print_kmer(canonR);
							//~ cout<<2<<endl;
							//~ cout<<pos<<endl;
							//~ cin.get();
							//~ print_kmer(seq);
							//~ print_kmer(seqR);
						}
					}
				}else{
					//~ cout<<"HERE"<<endl;
					#pragma omp atomic
					FP++;
					//~ cout<<7<<endl;
					//~ cin.get();
				}
			}
			//~ cout<<11<<endl;
			for(;i+k<query.size();++i){
				updateK(seq,query[i+k]);
				updateRCK(rcSeq,query[i+k]);
				canon=(min_k(seq, rcSeq));
				//COMPUTE KMER MINIMIZER
				uint32_t minimizer=(minimizer_according_xs(canon));
				if(all_mphf[minimizer/number_bucket_per_mphf]==NULL){
					hash=-1;
				}else{
					if(all_mphf[minimizer/number_bucket_per_mphf]->mphf_size==0){
						hash=-1;
					}else{
						hash=(all_mphf[minimizer/number_bucket_per_mphf]->kmer_MPHF.lookup(canon));
					}
				}
				if(hash<0){
					#pragma omp atomic
					FP++;
					//~ cout<<3<<endl;
					//~ cout<<"lol"<<endl;
					//~ cin.get();
					continue;
				}
				int n_bits_to_encode(all_mphf[minimizer/number_bucket_per_mphf]->bit_to_encode);
				uint pos(bool_to_int( n_bits_to_encode, hash, all_mphf[minimizer/number_bucket_per_mphf]->positions));
				if(2*(pos+k-1)>=all_buckets[minimizer]->bucketSeq.size()){
					#pragma omp atomic
					FP++;
					//~ cout<<"4"<<endl;
					continue;
				}
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
					if(not light_mode and not Valid_kmer[minimizer/number_superbuckets][j+1]){
						j+=k-1;
						if(2*(j+k)<all_buckets[minimizer]->bucketSeq.size()){
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
					//~ cout<<"5"<<endl;
					//~ cin.get();
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










