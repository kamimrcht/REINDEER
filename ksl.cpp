
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



using namespace std;



#define kmer __uint128_t



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


uint64_t str2num(const string& str){
	uint64_t res(0);
	for(uint i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}


uint32_t xs(uint32_t y){
	y^=(y<<13); y^=(y>>17);y=(y^=(y<<15)); return y;
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






void kmer_Set_Light::updateRCK(kmer& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-2));
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
	string ref,useless;
	string super_minimizer,kmer,minimizer,mmer;
	while(not inUnitigs.eof()){
		getline(inUnitigs,useless);
		getline(inUnitigs,ref);
		//FOREACH UNITIG
		if(not ref.empty() and not useless.empty()){
			super_minimizer=kmer=minimizer=mmer="";
			uint last_position(0);
			//FOREACH KMER
			uint i(0);
			for(;i+k<=ref.size();++i){
				kmer=getCanonical(ref.substr(i,k));
				minimizer=(getCanonical(kmer.substr(0,m)));
				//COMPUTE KMER MINIMIZER
				for(uint j(1);j+m<=kmer.size();++j){
					mmer=getCanonical(kmer.substr(j,m));
					minimizer=min_accordingtoXS(minimizer,mmer);
				}
				if(super_minimizer==""){
					super_minimizer=minimizer;
				}else{
					if(super_minimizer!=minimizer){
						*(out_files[xs(str2num(minimizer))%number_superbuckets])<<">"+(minimizer)+"\n"<<ref.substr(last_position,i-last_position+k-1)<<"\n";
						last_position=i;
						super_minimizer=minimizer;
					}
				}
			}
			if(i-last_position+k-1!=0){
				*(out_files[xs(str2num(minimizer))%number_superbuckets])<<">"+(minimizer)+"\n"<<ref.substr(last_position,i-last_position+k-1)<<"\n";
			}
		}
	}
}



void kmer_Set_Light::read_super_buckets(const string& input_file){
	string useless,line;
	uint number_superbuckets_by_buckets(minimizer_number/number_superbuckets);
	for(uint SBC(0);SBC<number_superbuckets;++SBC){
		ifstream in(input_file+to_string(SBC));
		while(not in.eof()){
			getline(in,useless);
			getline(in,line);
			useless=useless.substr(1);
			//~ cout<<useless<<endl;
			uint minimizer(str2num(useless));
			//~ cout<<minimizer<<" "<<buckets.size()<<endl;;
			buckets[minimizer]+=line;
		}
	}
}



void kmer_Set_Light::create_mphf(){
	//~ cout<<buckets.size()<<endl;
	//~ cin.get();
	for(uint BC(0);BC<buckets.size();++BC){
		//~ cout<<"bucket go"<<endl;
		auto anchors=new vector<kmer>;
		kmer seq(str2num(buckets[BC].substr(0,k))),rcSeq(rcb(seq,k)),canon(min(seq,rcSeq));
		anchors->push_back(canon);
		for(uint j(1);j+k<buckets[BC].size();++j){
			if(buckets[BC][j+k]=='\n'){
				cout<<"endl"<<endl;
				j+=k;
				seq=(str2num(buckets[BC].substr(j,k))),rcSeq=(rcb(seq,k)),canon=(min(seq,rcSeq));
				anchors->push_back(canon);
			}else{
				updateK(seq,buckets[BC][j+k]);
				updateRCK(rcSeq,buckets[BC][j+k]);
				canon=(min(seq, rcSeq));
				anchors->push_back(canon);
			}

		}
		//~ cout<<"gomphf"<<endl;
		//~ cout<<anchors->size()<<endl;;
		auto data_iterator3 = boomphf::range(static_cast<const kmer*>(&(*anchors)[0]), static_cast<const kmer*>((&(*anchors)[0])+anchors->size()));
		kmer_MPHF[BC]= boomphf::mphf<kmer,hasher>(anchors->size(),data_iterator3,1,gammaFactor,false);
		//~ cout<<"sucess"<<endl;
		buckets_size.push_back(anchors->size());

		delete anchors;
	}
	//~ cout<<"the end"<<endl;
}



void kmer_Set_Light::fill_positions(){
	for(uint BC(0);BC<buckets.size();++BC){
		//~ cout<<1<<endl;
		positions[BC].resize(buckets_size[BC]);
		kmer seq(str2num(buckets[BC].substr(0,k))),rcSeq(rcb(seq,k)),canon(min(seq,rcSeq));
		//~ cout<<(uint64_t)seq<<endl;
		//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
		positions[BC][kmer_MPHF[BC].lookup(canon)]=0;
		for(uint j(1);j+k<buckets[BC].size();++j){
			//~ cout<<2<<endl;
			updateK(seq,buckets[BC][j+k]);
			updateRCK(rcSeq,buckets[BC][j+k]);
			canon=(min(seq, rcSeq));
			positions[BC][kmer_MPHF[BC].lookup(canon)]=j;
		}
	}
}



bool kmer_Set_Light::exists(const string& query){
	uint32_t mini;
	for(uint j(0);j<query.size();++j){
		kmer seq(str2num(query.substr(0,m))),rcSeq(rcb(seq,m)),canon(min(seq,rcSeq));
		mini=canon;
		for(uint ii(0);ii+k<query.size();++ii){
			updateK(seq,query[ii+k]);
			updateRCK(rcSeq,query[ii+k]);
			canon=(min(seq, rcSeq));
			mini=min(mini,(uint32_t)canon);
		}
	}
	uint32_t pos( positions[mini][kmer_MPHF[mini].lookup(str2num(query))]);
	if(query.substr(pos,k)==query){
		return true;
	}
	if(query.substr(pos,k)==revComp(query)){
		return true;
	}
	return false;

}
