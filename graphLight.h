#ifndef BLIGHT
#define BLIGHT



#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include "bbhash.h"



using namespace std;



typedef __uint128_t kmer;
typedef boomphf::SingleHashFunctor<kmer>  hasher;
typedef boomphf::mphf<  kmer, hasher  > MPHF;



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string reverseComplements(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i]= revCompChar(s[i]);
		// rc[s.size()-1-i]=char2int[(uint)s[i]];
	}
	return rc;
}



string getRepresent(const string& s){
	return min(s,reverseComplements(s));
}



kmer seq2int(const string& seq){
	kmer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		switch(seq[i]){
			case 'A':
				break;
			case 'C':
				res+=1;
				break;
			case 'G':
				res+=2;
				break;
			case 'T':
				res+=3;
				break;
		}
	}
}



//TODO CAN BE IMPROVED
kmer seq2intCanon(const string& seqnoncanon){
	string seq(getRepresent(seqnoncanon));
	kmer res(0);
	for(uint i(0);i<seq.size();++i){
		res<<=2;
		switch(seq[i]){
			case 'A':
				break;
			case 'C':
				res+=1;
				break;
			case 'G':
				res+=2;
				break;
			case 'T':
				res+=3;
				break;
		}
	}
}



class graphLight{
public:
	//UNITIGS IN BINARY
	vector<bool> unitigs;
	//POSITIONS OF EACH UNIIG IN THE BITARRAY
	vector<uint> positions;
	//THE NEIBORS OF A UNITIG
	vector<uint8_t> neighbors;
	//INDEX FIRST AND LAST KMER OF EACH UNITIG
	MPHF kmer2unitig;
	//UNITIG NUMBER FROM HASH
	vector<uint> unitigNumber;
	uint kmerSize;



	kmer boolsToKmerbeg(const uint hashIndice){
		uint positionBeg(positions[hashIndice]);
		//~ uint positionEnd(positions[hashIndice+1]);
		kmer res(0);
		for(uint i(positionBeg);i<2*kmerSize;i+=2){
			res<<=2;
			res+=unitigs[i];
			res+=(unitigs[i+1]<<1);
		}
		return res;
	}



	kmer boolsToKmerend(const uint hashIndice){
		//~ uint positionBeg(positions[hashIndice]);
		uint positionEnd(positions[hashIndice+1]);
		kmer res(0);
		for(uint i(positionEnd-2*kmerSize);i<positionEnd;i+=2){
			res<<=2;
			res+=unitigs[i];
			res+=(unitigs[i+1]<<1);
		}
		return res;
	}



	graphLight(const string& unitigFileName, const uint kmerSize,const uint gamma=5,const uint coreNumber=4){
		string nuc("ACTG");
		vector<kmer> kmerToIndex;
		ifstream unitigFile(unitigFileName.c_str());
		string header,unitig;
		uint unitigCounter;
		//POSITIONS AND BIT ARRAY
		cout<<"first loop"<<endl;
		while(not unitigFile.eof()){
			getline(unitigFile,header);
			getline(unitigFile,unitig);
			cout<<unitig<<endl;
			if(unitig.size()>=kmerSize){
				//FILLING THE POSITIONS of the unitig
				positions.push_back(unitigs.size());
				//FILLING THE BITARRAY
				for(uint i(0);i<unitig.size();++i){
					switch(unitig[i]){
						case 'A':
							unitigs.push_back(false);
							unitigs.push_back(false);
							break;
						case 'C':
							unitigs.push_back(false);
							unitigs.push_back(true);
							break;
						case 'G':
							unitigs.push_back(true);
							unitigs.push_back(false);
							break;
						case 'T':
							unitigs.push_back(true);
							unitigs.push_back(true);
							break;
					}
				}
				cout<<"lol"<<endl;
				cout<<unitigs.size()<<endl;;
				kmer beg(seq2intCanon(unitig.substr(0,kmerSize)));
				kmer end(seq2intCanon(unitig.substr(unitig.size()-kmerSize,kmerSize)));
				if(beg==end){
					kmerToIndex.push_back(beg);
				}else{
					kmerToIndex.push_back(beg);
					kmerToIndex.push_back(end);
				}
			}
			unitig=header="";
		}
		positions.push_back(unitigs.size());

		//CREATING MPHF
		auto data_iterator = boomphf::range(static_cast<const kmer*>(&kmerToIndex[0]), static_cast<const kmer*>((&kmerToIndex[0])+kmerToIndex.size()));
		kmer2unitig= boomphf::mphf<kmer,hasher>(kmerToIndex.size(),data_iterator,coreNumber,gamma,false);
		kmerToIndex.clear();

		//FILLING UNITIG NUMBERS
		cout<<"second loop"<<endl;
		unitigFile.seekg(0, ios::beg);
		unitigFile.clear();
		uint count(0);
		while(not unitigFile.eof()){
			getline(unitigFile,header);
			getline(unitigFile,unitig);
			if(unitig.size()>=kmerSize){
				cout<<count++<<endl;
				kmer beg(seq2int(unitig.substr(0,kmerSize)));
				kmer end(seq2int(unitig.substr(unitig.size()-kmerSize,kmerSize)));
				if(beg==end){
					unitigNumber[kmer2unitig.lookup(beg)]=count++;
				}else{
					unitigNumber[kmer2unitig.lookup(beg)]=count;
					unitigNumber[kmer2unitig.lookup(end)]=count++;
				}
			}
			unitig=header="";
		}

		//FILLING NEIGHBORS
		cout<<"third loop"<<endl;
		unitigFile.seekg(0, ios::beg);
		unitigFile.clear();
		string subbeg,subend;
		while(not unitigFile.eof()){
			getline(unitigFile,header);
			getline(unitigFile,unitig);
			if(unitig.size()>=kmerSize){
				subbeg=(unitig.substr(1,kmerSize-1));
				subend=(unitig.substr(unitig.size()-kmerSize+1,kmerSize-1));
				for(uint i(0);i<4;++i){
					kmer nextkmer(seq2int(subbeg+nuc[i]));
					int hashN(kmer2unitig.lookup(nextkmer));
					if(hashN>=0){
						kmer kmerToCheckbeg(boolsToKmerbeg(hashN));
						if(kmerToCheckbeg==nextkmer){
							neighbors[count-1]+=(1<<i);
						}else{
							kmer kmerToCheckend(boolsToKmerend(hashN));
							if(kmerToCheckbeg==nextkmer){
								neighbors[count-1]+=(1<<i);
							}
						}
					}
					kmer prevkmer(seq2int(nuc[i]+subend));
					int hashP(kmer2unitig.lookup(prevkmer));
					if(hashP>=0){
						kmer kmerToCheckbeg(boolsToKmerbeg(hashP));
						if(kmerToCheckbeg==nextkmer){
							neighbors[count-1]+=(1<<(i+4));
						}else{
							kmer kmerToCheckend(boolsToKmerend(hashP));
							if(kmerToCheckbeg==nextkmer){
								neighbors[count-1]+=(1<<(i+4));
							}
						}
					}
				}
			}
			unitig=header="";
		}
		cout<<"THEPRINT"<<endl;
		cout<<unitigs.size()<<endl;;
		for(uint i(0);i<unitigs.size();++i){
			cout<<(uint)unitigs[i];
		}
		cout<<endl;
	}

string getUnitig(const uint n){
	string res;
	uint positionBeg(positions[n]);
	uint positionEnd(positions[n+1]);
	cout<<positionBeg<<endl;
	for(uint i(positionBeg);i<positionEnd;i+=2){
		char nuc(0);
		nuc+=unitigs[i+1];
		nuc<<=1;
		nuc+=unitigs[i];
		switch(nuc){
			case 0:
				res.push_back('A');
				break;
			case 1:
				res.push_back('C');
				break;
			case 2:
				res.push_back('G');
				break;
			case 3:
				res.push_back('T');
				break;
		}
	}
	return res;
}
};



#endif
