
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
#define kmer uint64_t



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


uint64_t xs(uint64_t y){
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

void kmer_Set_Light::updateM(kmer& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=(1<<(2*m));
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
	min+=(nuc2intrc(nuc)<<(2*m-2));
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


 hash<kmer> k_hash;


uint32_t kmer_Set_Light::minimizer_according_xs(kmer seq){
	uint32_t mini,mmer;
	mini=seq%minimizer_number;
	mini=min(mini,(uint32_t)rcb(mini,m));
	for(uint i(0);i<k-m;i++){
		seq>>=2;
		mmer=seq%minimizer_number;
		mmer=min(mmer,(uint32_t)rcb(mmer,m));
		if(abundance_minimizer[mini]>abundance_minimizer[mmer]){
			mini=mmer;
		}else{
			if(abundance_minimizer[mini]==abundance_minimizer[mmer]){
				mini=((xs(mini)<=xs(mmer))?mini:mmer);
			}
		}
		//~ min=((xs(min)<=xs(mmer))?min:mmer);
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
		//~ cout<<"go"<<endl;
		getline(inUnitigs,useless);
		getline(inUnitigs,ref);
		//FOREACH UNITIG
		if(not ref.empty() and not useless.empty()){
			uint last_position(0);
			//FOREACH KMER
			kmer seq(str2num(ref.substr(0,m))),rcSeq(rcb(seq,m)),canon(min_k(seq,rcSeq));
			//~ cout<<canon<<endl;
			//~ cout<<abundance_minimizer.size()<<endl;
			abundance_minimizer[canon]++;
			uint i(0);
			for(;i+m<ref.size();++i){
				updateM(seq,ref[i+m]);
				updateRCM(rcSeq,ref[i+m]);
				canon=(min_k(seq,rcSeq));
				//~ cout<<canon<<endl;
				//~ cout<<abundance_minimizer.size()<<endl;
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
		kmer super_minimizer,minimizer;
		while(not inUnitigs.eof()){
			#pragma omp critical(dataupdate)
			{
				getline(inUnitigs,useless);
				getline(inUnitigs,ref);
			}
			//~ cout<<ref<<endl;
			//FOREACH UNITIG
			if(not ref.empty() and not useless.empty()){
				super_minimizer=minimizer=minimizer_number;
				uint last_position(0);
				//FOREACH KMER
				kmer seq(str2num(ref.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
				//~ minimizer=min(minimizer_according_xs(seq),minimizer_according_xs(rcSeq));
				minimizer=minimizer_according_xs(canon);
				super_minimizer=minimizer;
				uint i(0);
				for(;i+k<ref.size();++i){
					updateK(seq,ref[i+k]);
					updateRCK(rcSeq,ref[i+k]);
					canon=(min_k(seq, rcSeq));
					//COMPUTE KMER MINIMIZER
					//~ minimizer=min(minimizer_according_xs(seq),minimizer_according_xs(rcSeq));
					minimizer=minimizer_according_xs(canon);
					if(super_minimizer!=minimizer){
						omp_set_lock(&(lock[xs((super_minimizer))%number_superbuckets]));
						*(out_files[xs((super_minimizer))%number_superbuckets])<<">"+to_string(super_minimizer)+"\n"<<ref.substr(last_position,i-last_position+k)<<"\n";
						last_position=i+1;
						omp_unset_lock(&(lock[xs((super_minimizer))%number_superbuckets]));
						super_minimizer=minimizer;
					}
				}
				if(ref.size()-last_position>k-1){
					//~ cout<<"NADINE"<<endl;
					//~ cout<<ref.substr(last_position)<<endl;
					//~ cout<<minimizer<<endl;
					omp_set_lock(&(lock[xs((super_minimizer))%number_superbuckets]));
					*(out_files[xs((super_minimizer))%number_superbuckets])<<">"+to_string(super_minimizer)+"\n"<<ref.substr(last_position)<<"\n";
					omp_unset_lock(&(lock[xs((super_minimizer))%number_superbuckets]));
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
	//~ cout<<str<<mini<<endl;
	//~ cout<<str<<endl;
	//~ cout<<str<<endl;
	//~ cout<<bucketSeq[mini].size()<<endl;
	for(uint i(0);i<str.size();++i){
		//~ cout<<"go"<<endl;
		Valid_kmer[mini].push_back(true);
		switch (str[i]){
			case 'A':
				bucketSeq[mini].push_back(false);
				bucketSeq[mini].push_back(false);
				//~ cout<<"00";
				break;
			case 'C':
				bucketSeq[mini].push_back(false);
				bucketSeq[mini].push_back(true);
				//~ cout<<"01";
				break;
			case 'G':
				bucketSeq[mini].push_back(true);
				bucketSeq[mini].push_back(false);
				//~ cout<<"10";
				break;
			case 'T':
				bucketSeq[mini].push_back(true);
				bucketSeq[mini].push_back(true);
				//~ cout<<"11";
				break;
			default:
				cout<<"nope"<<endl;
			}
	}
	//~ cout<<"endl"<<endl;
	//~ cout<<bucketSeq[mini].size()<<endl;
	//~ cout<<endl;
	for(uint i(0);i<k-1;++i){
		Valid_kmer[mini][Valid_kmer[mini].size()-k+i+1]=(false);
	}
	//~ for(uint i(0);i<Valid_kmer[mini].size();++i){
		//~ cout<<Valid_kmer[mini][i];
	//~ }
	//~ cout<<endl;cin.get();
	//~ Valid_kmer[mini].insert(Valid_kmer[mini].end()-k+1,k-1,false);
}


void kmer_Set_Light::read_super_buckets(const string& input_file){
	uint64_t total_size(0);
	uint number_superbuckets_by_buckets(minimizer_number/number_superbuckets);
	#pragma omp parallel num_threads(coreNumber)
	{
		string useless,line;
		#pragma omp for
		for(uint SBC=0;SBC<number_superbuckets;++SBC){
			//~ cin.get();
			//~ cout<<SBC<<endl;
			ifstream in(input_file+to_string(SBC));
			while(not in.eof()){
				getline(in,useless);
				getline(in,line);
				if(not useless.empty()){
					useless=useless.substr(1);
					uint minimizer(stoi(useless));
					//~ buckets[minimizer]+=(line+'\n');
					//~ cout<<1<<endl;
					//~ if(line.substr(0,k+1)==SEQOFINTERESET){
						//~ cout<<"YES"<<endl;
						//~ cout<<str2num(SEQOFINTERESET)<<endl;
						//~ str2bool(line,minimizer);
						//~ cin.get();
					//~ }
					str2bool(line,minimizer);
					//~ cout<<2<<endl;
					#pragma omp atomic
					total_size+=line.size();
				}
			}
		}
	}
	cout<<"Total size of the partitionned graph: "<<intToString(total_size)<<endl;
}


kmer kmer_Set_Light::get_kmer(uint32_t mini,uint32_t pos){
	kmer res(0);
	//~ pos<<=1;
	//~ kmer acc(1);
	for(uint i(0);i<k;++i){
		res<<=2;
		//~ cout<<"nad"<<1<<endl;
		//~ cout<<pos<<" "<<mini<<" "<<bucketSeq[mini].size()<<endl;
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
				//~ res+=0*acc;
			}
		}

	}
	return res;
}


kmer kmer_Set_Light::update_kmer(uint32_t pos,uint32_t mini,kmer input){
	input<<=2;
	//~ pos<<=1;
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
	//~ input%=offsetUpdateAnchor;
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




void kmer_Set_Light::create_mphf(){
	uint64_t max_size(0);
	uint64_t anchors_number(0);
	//~ cout<<"go"<<endl;
	#pragma omp parallel for num_threads(1)
	for(uint BC=(0);BC<bucketSeq.size();++BC){
		if(bucketSeq[BC].size()==0){
			buckets_size[BC]=0;
			continue;
		}
				//~ <<"begB"<<endl;

		//~ cout<<buckecoutts[BC]<<endl;
		//~ cout<<"endB"<<endl;
		//~ cout<<bucketSeq[BC].size()/2<<endl;
		//~ cout<<Valid_kmer[BC].size()<<endl;cin.get();
		auto anchors=new vector<kmer>;
		//~ cout<<1<<endl;
		//~ kmer seq2(str2num(buckets[BC].substr(0,k))),rcSeq2(rcb(seq2,k)),canon2(min_k(seq2,rcSeq2));
		//~ kmer seq2(str2num(buckets[BC].substr(0,k)));
		kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
		//~ if(str2num(SEQOFINTERESET)==seq){
						//~ cout<<"stop"<<endl;
						//~ cout<<0<<endl;
						//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
						//~ cin.get();
					//~ }
		//~ if(str2num(SEQOFINTERESET)==seq2 or buckets[BC].substr(0,k)==SEQOFINTERESET){
						//~ cout<<SEQOFINTERESET<<endl;
						//~ cout<<0<<endl;
						//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
						//~ cin.get();
					//~ }
		//~ cout<<canon<<endl;
		//~ cout<<canon2<<endl;
		//~ cin.get();
		//~ cout<<2<<endl;
		//~ print_kmer(seq);
		anchors->push_back(canon);
		for(uint j(0);2*(j+k)<bucketSeq[BC].size();j++){
			//~ cout<<3<<endl;
			//~ cout<<j<<endl;
			if(not Valid_kmer[BC][j+1]){
				//~ cout<<4<<endl;
				j+=k-1;
				//~ cout<<j<<endl;
				if(2*(j+k)<bucketSeq[BC].size()){
					//~ cout<<
					//~ seq=(str2num(buckets[BC].substr(j,k))),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
					seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
					//~ if(str2num(SEQOFINTERESET)==seq){
						//~ cout<<"stop"<<endl;
						//~ cout<<0<<endl;
						//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
						//~ cin.get();
					//~ }
					//~ cout<<"kmer after jump"<<endl;
					//~ print_kmer(seq);
					//~ cin.get();
					anchors->push_back(canon);
					//~ --j;
				}
			}else{
				//~ cout<<5<<endl;
				//~ updateK(seq,bucketSeq[BC][j+k]);
				seq=update_kmer(j+k,BC,seq);
				//~ kmer seq2=(str2num(buckets[BC].substr(j+1,k)));
					//~ cout<<seq<<endl;
					//~ cout<<seq2<<endl;
					//~ cin.get();
				//~ updateRCK(rcSeq,buckets[BC][j+k]);
				rcSeq=(rcb(seq,k));
				canon=(min_k(seq, rcSeq));
				//~ print_kmer(seq);
				//~ if(str2num(SEQOFINTERESET)==seq){
						//~ cout<<"stop"<<endl;
						//~ cout<<0<<endl;
						//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
						//~ cin.get();
					//~ }
				anchors->push_back(canon);

			}
		}
		//~ cout<<"nadine"<<endl;
		//~ cout<<anchors->size()<<endl;
		max_size=max(max_size,anchors->size());
		anchors_number+=anchors->size();
		auto data_iterator3 = boomphf::range(static_cast<const kmer*>(&(*anchors)[0]), static_cast<const kmer*>((&(*anchors)[0])+anchors->size()));
		kmer_MPHF[BC]= boomphf::mphf<kmer,hasher>(anchors->size(),data_iterator3,coreNumber,gammaFactor,false);
		buckets_size[BC]=anchors->size();

		delete anchors;
	}
	cout<<"Largest bucket: "<<intToString(max_size)<<endl;
	cout<<"Number of kmer: "<<intToString(anchors_number)<<endl;
}



void int_to_bool(uint n_bits_to_encode,uint32_t X, uint32_t pos,vector<bool>& res){
	//~ cout<<"ITB"<<endl;
	for(uint i(0);i<n_bits_to_encode;++i){
		//~ cout<<X<<endl;
		res[i+pos*n_bits_to_encode]=X%2;
		//~ cout<<X%2;
		//~ cout<<"pos"<<i+pos*n_bits_to_encode<<endl;
		X>>=1;
	}
	//~ cout<<endl;
}


uint32_t bool_to_int(uint n_bits_to_encode,uint pos,const vector<bool>& V){
	//~ cout<<"BTI"<<pos<<endl;
	uint32_t res(0);
	uint32_t acc(1);
	for(uint i(0);i<n_bits_to_encode;++i){
		//~ cout<<i<<endl;
		//~ cout<<pos<<endl;
		//~ cout<<i+pos*n_bits_to_encode<<endl;
		//~ cout<<V.size()<<endl;
		//~ cout<<"print bit"<<endl;
		if(V[i+pos*n_bits_to_encode]){
			//~ cout<<1;
			res+=acc;
		}else{
			//~ cout<<0;
		}
		;
		acc<<=1;
	}
	//~ cout<<endl;
	return res;
}



void kmer_Set_Light::fill_positions(){
	//~ cout<<"FILL POSITION"<<endl;
	uint64_t total_size(0);
	#pragma omp parallel for num_threads(coreNumber)
	for(uint BC=(0);BC<bucketSeq.size();++BC){
		if(buckets_size[BC]==0){
			continue;
		}
		//~ cout<<1<<endl;
		int n_bits_to_encode((ceil(log2(bucketSeq[BC].size()+1))-bit_saved_sub));
		//~ cout<<BC<<endl;
		//~ cout<<bucketSeq[BC].size()<<endl;
		//~ cout<<Valid_kmer[BC].size()<<endl;
		//~ cout<<n_bits_to_encode<<endl;
		//~ cout<<"go?"<<endl;
		if(n_bits_to_encode<1){n_bits_to_encode=1;}
		positions[BC].resize(buckets_size[BC]*n_bits_to_encode,0);
		//~ cout<<2<<endl;
		#pragma omp atomic
		total_size+=buckets_size[BC]*n_bits_to_encode;
		//~ kmer seq(str2num(buckets[BC].substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
		kmer seq(get_kmer(BC,0)),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq));
		//~ cout<<"go"<<1<<endl;
		//~ cout<<3<<endl;
		//~ if(str2num(SEQOFINTERESET)==seq){
						//~ cout<<"stop"<<endl;
						//~ cout<<0<<endl;
						//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
						//~ cout<<n_bits_to_encode<<endl;
						//~ cin.get();
					//~ }
		for(uint j(0);2*(j+k)<bucketSeq[BC].size();j++){
			//~ cout<<4<<endl;
			//~ cout<<"go"<<2<<endl;
			if(not Valid_kmer[BC][j+1]){
				//~ cout<<"go"<<3<<endl;
				j+=k-1;
				if(2*(j+k)<bucketSeq[BC].size()){
					//~ cout<<"go"<<4<<endl;
					//~ if(str2num(SEQOFINTERESET)==seq){
						//~ cout<<"stop"<<endl;
						//~ cout<<j+1<<endl;
						//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
						//~ cout<<n_bits_to_encode<<endl;
						//~ cin.get();
					//~ }
					//~ seq=(str2num(buckets[BC].substr(j,k))),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
					seq=(get_kmer(BC,j+1)),rcSeq=(rcb(seq,k)),canon=(min_k(seq,rcSeq));
					//~ --j;
					int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,kmer_MPHF[BC].lookup(canon),positions[BC]);
				}
			}else{
				//~ if(str2num(SEQOFINTERESET)==seq){
						//~ cout<<"stop"<<endl;
						//~ cout<<j+1<<endl;
						//~ cout<<kmer_MPHF[BC].lookup(canon)<<endl;
						//~ cout<<n_bits_to_encode<<endl;
						//~ cin.get();
					//~ }
				//~ cout<<"go"<<5<<endl;
				//~ updateK(seq,buckets[BC][j+k]);
				//~ updateRCK(rcSeq,buckets[BC][j+k]);
				seq=update_kmer(j+k,BC,seq);
				rcSeq=(rcb(seq,k));
				canon=(min_k(seq, rcSeq));
				//~ print_kmer(seq);
				//~ cout<<j+1<<endl;
				int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,kmer_MPHF[BC].lookup(canon),positions[BC]);
				//~ cout<<"nadine"<<endl;
			}
		}
	}
	cout<<"Total Positions size: (MBytes) "<<intToString(total_size/(8*1024*1024))<<endl;
}



int64_t kmer_Set_Light::correct_pos(uint32_t mini, uint64_t p){
	//~ cout<<"correct"<<endl;
	//~ cout<<p<<endl;
	//~ cout<<Valid_kmer[mini].size()<<endl;
	for(uint i(0);i<k;i++){
		//~ cout<<p+i<<endl;
		if(Valid_kmer[mini][p+i]){
			return (p+i);
		}else{
			//~ cout<<p+i<<"FALSE"<<endl;
		}
	}
	return p;
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
		if(query.size()>=k){
			//~ cout<<"QUERY"<<endl;
			kmer seq(str2num(query.substr(0,k))),rcSeq(rcb(seq,k)),canon(min_k(seq,rcSeq)),canonR,seqR,rcSeqR;
			//~ print_kmer(seq);
			uint i(0);
			canon=(min_k(seq, rcSeq));
			//COMPUTE KMER MINIMIZER
			//~ uint32_t minimizer=min(minimizer_according_xs(seq),minimizer_according_xs(rcSeq));
			uint32_t minimizer=(minimizer_according_xs(canon));
			int32_t hash(kmer_MPHF[minimizer].lookup(canon));
			if(hash<0){
				#pragma omp atomic
				FP++;
				cout<<"fail MPHF"<<endl;
				cin.get();
			}else{
				//~ TP++;
				//~ continue;
				int n_bits_to_encode((ceil(log2(bucketSeq[minimizer].size()+1))-bit_saved_sub));
				if(n_bits_to_encode<1){n_bits_to_encode=1;}
				uint pos(bool_to_int( n_bits_to_encode, hash, positions[minimizer]));
				//~ cout<<pos<<endl;
				pos=correct_pos(minimizer,(positions_to_check)*pos);
				//~ cout<<"first1"<<endl;
				//~ seqR=str2num(buckets[minimizer].substr(pos,k));
				seqR=get_kmer(minimizer,pos);
				//~ cout<<minimizer<<" "<<pos<<endl;
				//~ print_kmer(seqR);
				rcSeqR=rcb(seqR,k);
				canonR=(min_k(seqR, rcSeqR));
				if(canon==canonR){
					#pragma omp atomic
					TP++;
					//~ cout<<"sucess"<<endl;
				}else{
					uint j;
					bool found(false);
					for(j=(pos);j<pos+positions_to_check;++j){
						//~ if(buckets[minimizer][j+k]=='\n'){
						if(not Valid_kmer[minimizer][j+1]){
							j+=k-1;
							if(2*(j+k)<bucketSeq[minimizer].size()){
								//~ cout<<"first2"<<endl;
								seqR=(get_kmer(minimizer,j)),rcSeqR=(rcb(seqR,k)),canon=(min_k(seqR,rcSeqR));
								//~ --j;
								canonR=(min_k(seqR, rcSeqR));
								//~ int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,kmer_MPHF[BC].lookup(canon),positions[BC]);
							}
						}else{
							//~ updateK(seqR,buckets[minimizer][j+k]);
							//~ updateRCK(rcSeqR,buckets[minimizer][j+k]);
							seqR=update_kmer(j+k,minimizer,seq);
							rcSeqR=(rcb(seqR,k));
						}
						canonR=(min_k(seqR, rcSeqR));
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
						//~ cout<<2<<endl;
					}
				}
			}
			for(;i+k<query.size();++i){
				//~ cout<<"the res"<<endl;
				cout<<"QUREY"<<endl;
				updateK(seq,query[i+k]);
				print_kmer(seq);
				updateRCK(rcSeq,query[i+k]);
				canon=(min_k(seq, rcSeq));
				//COMPUTE KMER MINIMIZER
				//~ uint32_t minimizer=min(minimizer_according_xs(seq),minimizer_according_xs(rcSeq));
				uint32_t minimizer=(minimizer_according_xs(canon));
				int32_t hash(kmer_MPHF[minimizer].lookup(canon));
				if(hash<0){
					#pragma omp atomic
					FP++;
					cout<<"fail MPHF"<<endl;
					continue;
				}
				//~ TP++;
				//~ continue;
				int n_bits_to_encode((ceil(log2(bucketSeq[minimizer].size()+1))-bit_saved_sub));
				if(n_bits_to_encode<1){n_bits_to_encode=1;}
				//~ cout<<hash<<endl;
				//~ cout<<n_bits_to_encode<<endl;
				uint pos(bool_to_int( n_bits_to_encode, hash, positions[minimizer]));
				//~ cout<<pos<<" "<<minimizer<<endl;
				pos=correct_pos(minimizer,(positions_to_check)*pos);
				//~ seqR=str2num(buckets[minimizer].substr(pos,k));
				//~ cout<<pos<<endl;
				seqR=(get_kmer(minimizer,pos));
				print_kmer(seqR);
				rcSeqR=rcb(seqR,k);
				canonR=(min_k(seqR, rcSeqR));
				if(canon==canonR){
					#pragma omp atomic
					TP++;
					//~ cout<<"sucess"<<endl;
					continue;
				}
				cout<<TP<<endl;
				cout<<"FAIL"<<endl;
				cin.get();
				uint j;
				bool found(false);
				for(j=(pos);j<pos+positions_to_check;++j){
					if(not Valid_kmer[minimizer][j+1]){
						//~ cout<<"endl"<<endl;
						j+=k-1;
						if(2*(j+k)<bucketSeq[minimizer].size()){
							//~ seqR=(str2num(buckets[minimizer].substr(j,k))),rcSeqR=(rcb(seqR,k)),canonR=(min_k(seqR,rcSeqR));
							seqR=(seqR=(get_kmer(minimizer,j))),rcSeqR=(rcb(seqR,k)),canonR=(min_k(seqR,rcSeqR));
							//~ --j;
							canonR=(min_k(seqR, rcSeqR));
							//~ int_to_bool(n_bits_to_encode,(j+1)/positions_to_check,kmer_MPHF[BC].lookup(canon),positions[BC]);
						}
					}else{
						//~ updateK(seqR,buckets[minimizer][j+k]);
						//~ updateRCK(rcSeqR,buckets[minimizer][j+k]);
						seq=update_kmer(j+k,minimizer,seq);
						rcSeq=(rcb(seq,k));
						canonR=(min_k(seqR, rcSeqR));
					}
					if(canon==canonR){
						#pragma omp atomic
						TP++;
						found=true;
						break;
					}
				}
				//~ cout<<j<<" "<<epos+positions_to_check-1<<endl;
				if(not found){
					#pragma omp atomic
					FP++;
					//~ cout<<5<<endl;
					//~ cout<<pos<<endl;
				}
			}
		}
	}

	cout<<"Good kmer: "<<intToString(TP)<<endl;
	cout<<"Erroneous kmers: "<<intToString(FP)<<endl;
	cout<<"I am glad you are here with me. Here at the end of all things."<<endl;
}










