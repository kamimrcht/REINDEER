#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include <unordered_map>
#include <map>
#include <pthread.h>
#include <chrono>
#include <omp.h>
#include <tmmintrin.h>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <unistd.h>
#include "../blight/bbhash.h"
#include "../blight/blight.h"
#include "../blight/utils.h"
#include "../blight/zstr.hpp"
#include "../blight/common.h"
#include "../blight/robin_hood.h"
#include "../trle/trle.h"


using namespace std;
using namespace chrono;





uint64_t MAX_ABUNDANCE_DISCRETE;
vector<uint8_t> abundance_discretization;


uint16_t kmer_Set_Light::abundance_at (uint8_t index)
{
	if (index < abundance_discretization.size())
	{
		return floorf((abundance_discretization[index]  +  abundance_discretization[index+1])/2.0);
	}
	else
	{
		return 0;
	}
}



void kmer_Set_Light::init_discretization_scheme()
{
	MAX_ABUNDANCE_DISCRETE=65536;
	abundance_discretization.resize(257);
	int total =0;
	abundance_discretization[0] = 0;
	int idx=1;
	for(int ii=1; ii<= 70; ii++,idx++ )
	{
		total += 1;
		abundance_discretization[idx] = total ;
	}

	for(int ii=1; ii<= 15; ii++,idx++ )
	{
		total += 2;
		abundance_discretization[idx] = total  ;
	}

	for(int ii=1; ii<= 40; ii++,idx++ )
	{
		total += 10;
		abundance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 25; ii++,idx++ )
	{
		total += 20;
		abundance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 40; ii++,idx++ )
	{
		total += 100;
		abundance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 25; ii++,idx++ )
	{
		total += 200;
		abundance_discretization[idx] = total  ;
	}
	for(int ii=1; ii<= 40; ii++,idx++ )
	{
		total += 1000;
		abundance_discretization[idx] = total  ;
	}
	abundance_discretization[256] = total;
}


uint8_t kmer_Set_Light::return_count_bin(uint16_t abundance)
{
	int idx ;
	if (abundance >= MAX_ABUNDANCE_DISCRETE)
	{
		//~ _nb_abundances_above_precision++;
		//std::cout << "found abundance larger than discrete: " << abundance << std::endl;
		idx = abundance_discretization.size() -2 ;
	}
	else
	{
		//get first cell strictly greater than abundance
		auto  up = upper_bound(abundance_discretization.begin(), abundance_discretization.end(), abundance);
		up--; // get previous cell
		idx = up- abundance_discretization.begin() ;
	}
	return idx;
}




uint16_t kmer_Set_Light::parseCoverage_bin(const string& str){
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
	return return_count_bin((uint16_t)stof(str.substr(pos+5,i)));
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
		if( ((double)max(V1[i],V2[i])/((double)min(V1[i],V2[i])) >max_divergence_count)){
			return false;
		}
	}
	return true;
}



string compress_vector(const vector<uint16_t>& V){
	string res;
    unsigned char comp[V.size()*2+1024];
    uint32_t rle_length=trlec((unsigned char*)&V[0],V.size()*2,comp);
    res.assign((const char*)comp,rle_length);
	return res;
}


uint64_t max_size_bucket(0);
uint64_t min_size_bucket(10000000000);
uint64_t max_size_superbucket(0);
uint64_t min_size_superbucket(10000000000);




void kmer_Set_Light::construct_index_fof(const string& input_file, const string& tmp_dir, int colormode, double max_divergence){
	omp_set_nested(2);
	if(not tmp_dir.empty()){
		working_dir=tmp_dir+"/";
	}
	color_mode=colormode;
	if(color_mode==1){
		max_divergence_count=(max_divergence/100)+1;
	}
	if(color_mode==2){
		init_discretization_scheme();
	}
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
	vector<string> fnames;
	//STACK SUPERBUCKETS
	while(not infof.eof()){
		string file;
		getline(infof,file);
		if(not file.empty() and exists_test(file)){
			fnames.push_back(file);
		}
	}

	create_super_buckets_list(fnames);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	cout<<"Partition created	"<<intToString(read_kmer)<<" kmers read "<<endl;
	duration<double> time_span12 = duration_cast<duration<double>>(t2 - t1);
	cout<<time_span12.count() << " seconds."<<endl;
	{
		ofstream out(working_dir+"_blmonocolor.fa");
		#pragma omp parallel for num_threads(coreNumber)
		for(uint i_superbuckets=0; i_superbuckets<number_superbuckets.value(); ++i_superbuckets){
			merge_super_buckets_mem(working_dir+"_blout"+to_string(i_superbuckets),fnames.size(),&out);
			remove((working_dir+"_blsout"+to_string(i_superbuckets)).c_str());
			cout<<"-"<<flush;
		}
	}
	// cout<<endl;
	// cout<<"max size bucket:	"<<max_size_bucket<<endl;
	// cout<<"min size bucket:	"<<min_size_bucket<<endl;
	// cout<<"max size superbucket:	"<<max_size_superbucket<<endl;
	// cout<<"min size superbucket:	"<<min_size_superbucket<<endl;
	reset();


	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	cout<<"Monocolor minitig computed, now regular indexing start"<<endl;
	duration<double> time_span32 = duration_cast<duration<double>>(t3 - t2);
	cout<<time_span32.count() << " seconds."<<endl;
	// exit(0);

	create_super_buckets(working_dir+"_blmonocolor.fa");

	high_resolution_clock::time_point t4 = high_resolution_clock::now();
	duration<double> time_span43 = duration_cast<duration<double>>(t4 - t3);
	cout<<"Super buckets created: "<< time_span43.count() << " seconds."<<endl;

	read_super_buckets(working_dir+"_blout");
	delete [] nuc_minimizer ;
	delete [] start_bucket ;
	delete [] current_pos ;

	high_resolution_clock::time_point t5 = high_resolution_clock::now();
	duration<double> time_span53 = duration_cast<duration<double>>(t5 - t3);
	cout<<"Indexes created: "<< time_span53.count() << " seconds."<<endl;
	duration<double> time_spant = duration_cast<duration<double>>(t5 - t1);
	cout << "The whole indexing took me " << time_spant.count() << " seconds."<< endl;
}



vector<uint16_t> getcolorvector(const vector< pair<uint16_t,uint16_t> >&V,uint64_t number_color){
	vector<uint16_t> res(number_color,0);
	for(uint64_t i(0);i< V.size();++i){
		res[V[i].first]=V[i].second;
	}
	return res;
}


void kmer_Set_Light::merge_super_buckets_mem(const string& input_file, uint64_t number_color, ofstream* out,uint64_t number_pass ){
	bool toobig(false);
	for(uint pass(0);pass<number_pass;++pass){
		vector<robin_hood::unordered_node_map<kmer,kmer_context>> min2kmer2context(minimizer_number.value()/number_superbuckets.value());
		uint64_t mms=min2kmer2context.size();
		vector<int32_t> minimizers(mms,-1);
		ifstream in(input_file);
		uint64_t inserted_elements(0);

		// #pragma omp parallel num_threads(coreNumber)
		{
			string line,sequence;
			uint32_t color;
			uint16_t coverage;
			uint32_t min_int;
			uint8_t sizel(0);
            bool stop(false);
			while(not in.eof() and in.good() and not toobig){
				// #pragma omp critical
				{
					in.read(reinterpret_cast<char *>(&min_int), 4);
					in.read(reinterpret_cast<char *>(&color), 4);
					in.read(reinterpret_cast<char *>(&coverage), 2);
					in.read(reinterpret_cast<char *>(&sizel), 1);
					sequence.resize(sizel);
					in.read(reinterpret_cast<char *>(&sequence[0]), sizel);
                    if(in.eof()){
                        stop=true;
                    }else{
                        minimizers[min_int%mms]=(min_int);
                    }

				}
                if(not stop){
					uint64_t indice(min_int%mms);
					if( (number_pass!=1) and (indice%number_pass)!=pass){continue;}
					kmer seq(str2num(sequence.substr(0,k))),rcSeq(rcb(seq)),canon(min_k(seq,rcSeq));
					canon=(min_k(seq, rcSeq));
                    positions_mutex[indice%4096].lock();
					if(min2kmer2context[indice].count(canon)==0){
						min2kmer2context[indice][canon]={false,{}};
						// #pragma omp atomic
						inserted_elements++;
					}
					min2kmer2context[indice][canon].count.push_back({color,coverage});
					uint64_t sks=sequence.size();
					for(uint i(0);i+k<sks;++i){
						updateK(seq,sequence[i+k]);
						updateRCK(rcSeq,sequence[i+k]);
						canon=(min_k(seq, rcSeq));
						if(min2kmer2context[indice].count(canon)==0){
							min2kmer2context[indice][canon]={false,{}};
							// #pragma omp atomic
							inserted_elements++;
						}
						min2kmer2context[indice][canon].count.push_back({color,coverage});
					}
					if(inserted_elements>100000000){
						// #pragma omp critical
						toobig=true;
					}
					positions_mutex[indice%4096].unlock();
                }

			}
		}
		if (not toobig){
			get_monocolor_minitigs_mem(min2kmer2context,out,minimizers,number_color);
		}
	}
	number_pass*=2;
	if(toobig){
		merge_super_buckets_mem(input_file, number_color, out,number_pass);
	}
}



string nucleotides("ACGT");



kmer kmer_Set_Light::select_good_successor(const  robin_hood::unordered_node_map<kmer,kmer_context>& kmer2context,const kmer& start){
	kmer canon=canonize(start,k);
	if(kmer2context.count(canon)==0){return -1;}
	kmer_context kc(kmer2context.at(canon));
	for(uint64_t i(0);i<4;++i){
		kmer target=start;
		updateK(target,nucleotides[i]);
		kmer targetc=canonize(target,k);
		if(kmer2context.count(targetc)!=0){
			if(kmer2context.at(targetc).isdump){continue;}
            if(kmer2context.at(targetc).count==kc.count){
                return target;
            }
		}
	}
	return -1;
}



void kmer_Set_Light::get_monocolor_minitigs_mem(vector<robin_hood::unordered_node_map<kmer,kmer_context>>&  min2kmer2context , ofstream* out, const vector<int32_t>& mini,uint64_t number_color){
    vector<uint16_t> bit_vector(number_color,0);
    uint64_t size_superbucket(0);
    // #pragma omp parallel num_threads(coreNumber)
    {
    	string sequence, buffer,seq2dump,compact;
    	uint64_t ms=min2kmer2context.size();
    	// #pragma omp for schedule(static,ms/coreNumber)
    	for(uint i_set=(0);i_set<ms;i_set++){
    		for (auto& it: min2kmer2context[i_set]){
    			if(not it.second.isdump){
    				it.second.isdump=true;
    				auto colorV2dump(it.second.count);
    				kmer start=it.first;
    				seq2dump=kmer2str(start);
    				for(uint64_t step(0);step<2;step++){
    					while(true){
    						kmer next=select_good_successor(min2kmer2context[i_set],start);
    						if(next==(kmer)-1){break;}
    						compact=compaction(seq2dump,kmer2str(next),false);
    						if(compact.empty()){break;}
    						min2kmer2context[i_set].at(canonize(next,k)).isdump=true;
    						seq2dump=compact;
    						start=next;
    					}
    					start=rcb(it.first);
    					seq2dump=revComp(seq2dump);
    				}
    				buffer+=">"+to_string(mini[i_set])+color_coverage2str(getcolorvector(colorV2dump,number_color))+"\n"+seq2dump+"\n";
    				if(buffer.size()>8000){
    					#pragma omp critical (monocolorFile)
    					{
    						*out<<buffer;
    					}
    					buffer.clear();
    				}
    			}
    		}
		// 	uint64_t size_bucket=min2kmer2context[i_set].size();
        //     #pragma omp atomic
        //     size_superbucket+=size_bucket;
		// 	if(size_bucket>max_size_bucket){
		// 		// #pragma omp atomic
		// 		max_size_bucket=size_bucket;
		// 	}
		// 	if(size_bucket<min_size_bucket){
		// 		// #pragma omp atomic
		// 		min_size_bucket=size_bucket;
		// 	}
    	// 	min2kmer2context[i_set].clear();
    	}
    	#pragma omp critical (monocolorFile)
    	{
    		*out<<buffer<<flush;
    	}
    	buffer.clear();
    }
	// if(size_superbucket>max_size_superbucket){
	// 	max_size_superbucket=size_superbucket;
	// }
	// if(size_superbucket<min_size_superbucket){
	// 	min_size_superbucket=size_superbucket;
	// }
	//     cout<<size_superbucket/1000  <<' ';
}



uint16_t kmer_Set_Light::parseCoverage(const string& str){
	if(color_mode==0){return parseCoverage_bool(str);}
	if(color_mode==1){return parseCoverage_exact(str);}
	if(color_mode==2){return parseCoverage_bin(str);}
	return parseCoverage_log2(str);
}



void kmer_Set_Light::create_super_buckets_list(const vector<string>& input_files){
	struct rlimit rl;
	getrlimit (RLIMIT_NOFILE, &rl);
	rl.rlim_cur = number_superbuckets.value()+10+coreNumber;
	setrlimit (RLIMIT_NOFILE, &rl);
	// atomic<uint64_t> total_nuc_number(0);

	vector<ostream*> out_files;
	for(uint64_t i(0);i<number_superbuckets;++i){
		auto out =new  ofstream(working_dir+"_blout"+to_string(i),ofstream::app);
		out_files.push_back(out);
	}
	omp_lock_t lock[number_superbuckets.value()];
	for (uint64_t i=0; i<number_superbuckets; i++){
		omp_init_lock(&(lock[i]));
	}

	#pragma omp parallel num_threads(coreNumber)
    {
        #pragma omp for
    	for(uint32_t i_file=0;i_file<input_files.size();++i_file){
    		auto inUnitigsread=new zstr::ifstream(input_files[i_file]);
    		if(not inUnitigsread->good()){
    			cout<<"Problem with files opening"<<endl;
                delete inUnitigsread;
    			continue;
    		}

    		string ref,useless;
    		vector<string> buffer(number_superbuckets.value());
    		minimizer_type old_minimizer,minimizer;
    		while(not inUnitigsread->eof()){
    			ref=useless="";
    			getline(*inUnitigsread,useless);
    			getline(*inUnitigsread,ref);
    			if(ref.size()<k){
    				ref="";
    			}else{
    				#pragma omp atomic
    				read_kmer+=ref.size()-k+1;
    			}
    			//FOREACH UNITIG
    			if(not ref.empty() and not useless.empty()){
    				old_minimizer=minimizer=minimizer_number.value();
    				uint64_t last_position(0);
    				//FOREACH KMER
    				kmer seq(str2num(ref.substr(0,k)));
    				uint64_t position_min;
    				uint64_t min_seq=(str2num(ref.substr(k-minimizer_size_graph,minimizer_size_graph))),min_rcseq(rcbc(min_seq,minimizer_size_graph)),min_canon(min(min_seq,min_rcseq));
    				minimizer=regular_minimizer_pos(seq,position_min);
    				old_minimizer=minimizer;
    				uint64_t hash_min=unrevhash(minimizer);
    				uint64_t i(0);
    				for(;i+k<ref.size();++i){
    					updateK(seq,ref[i+k]);
    					updateM(min_seq,ref[i+k]);
    					updateRCM(min_rcseq,ref[i+k]);
    					min_canon=(min(min_seq,min_rcseq));
    					uint64_t new_h=unrevhash(min_canon);
    						//THE NEW mmer is a MINIMIZER
    					if(new_h<hash_min){
    						minimizer=(min_canon);
    						hash_min=new_h;
    						position_min=i+k-minimizer_size_graph+1;
    					}else{
    						//the previous minimizer is outdated
    						if(i>=position_min){
    							minimizer=regular_minimizer_pos(seq,position_min);
    							hash_min=unrevhash(minimizer);
    							position_min+=(i+1);
    						}else{
    						}
    					}
    					if(old_minimizer!=minimizer){
    						old_minimizer=(revhash(old_minimizer)%minimizer_number);
    						uint16_t cov(parseCoverage(useless));
    						uint8_t size_sk(i-last_position+k);
    						buffer[old_minimizer/bucket_per_superBuckets.value()].append(reinterpret_cast<char *> (&old_minimizer), 4);
    						buffer[old_minimizer/bucket_per_superBuckets.value()].append(reinterpret_cast<char *>(&i_file), 4);
    						buffer[old_minimizer/bucket_per_superBuckets.value()].append(reinterpret_cast<char *> (&cov), 2);
    						buffer[old_minimizer/bucket_per_superBuckets.value()].append(reinterpret_cast<char *> (&size_sk), 1);
    						buffer[old_minimizer/bucket_per_superBuckets.value()]+=ref.substr(last_position,i-last_position+k);
    						if(buffer[old_minimizer/bucket_per_superBuckets.value()].size()>80000){
    							omp_set_lock(&(lock[((old_minimizer))/bucket_per_superBuckets.value()]));
    							*(out_files[((old_minimizer))/bucket_per_superBuckets.value()])<<buffer[old_minimizer/bucket_per_superBuckets.value()];
    							omp_unset_lock(&(lock[((old_minimizer))/bucket_per_superBuckets.value()]));
    							buffer[old_minimizer/bucket_per_superBuckets.value()].clear();
    						}
    						last_position=i+1;
    						old_minimizer=minimizer;
    					}
    				}
    				if(ref.size()-last_position>k-1){
    					old_minimizer=(revhash(old_minimizer)%minimizer_number);
    					uint16_t cov(parseCoverage(useless));
    					uint8_t size_sk(ref.size()-last_position);
    					buffer[old_minimizer/bucket_per_superBuckets.value()].append((char *) &old_minimizer, 4);
    					buffer[old_minimizer/bucket_per_superBuckets.value()].append((char *) &i_file, 4);
    					buffer[old_minimizer/bucket_per_superBuckets.value()].append((char *) &cov, 2);
    					buffer[old_minimizer/bucket_per_superBuckets.value()].append((char *) &size_sk, 1);
    					buffer[old_minimizer/bucket_per_superBuckets.value()]+=ref.substr(last_position);
    					if(buffer[old_minimizer/bucket_per_superBuckets.value()].size()>80000){
    						omp_set_lock(&(lock[((old_minimizer))/bucket_per_superBuckets.value()]));
    						*(out_files[((old_minimizer))/bucket_per_superBuckets.value()])<<buffer[old_minimizer/bucket_per_superBuckets.value()];
    						omp_unset_lock(&(lock[((old_minimizer))/bucket_per_superBuckets.value()]));
    						buffer[old_minimizer/bucket_per_superBuckets.value()].clear();
    					}
    				}
    			}
    		}
    		for(uint64_t i(0);i<number_superbuckets.value();++i){
    			if(not buffer[i].empty()){
    				omp_set_lock(&(lock[i]));
    				*(out_files[i])<<buffer[i];
    				omp_unset_lock(&(lock[i]));
    			}
    		}
    		delete inUnitigsread;
    	}
    }
	for(uint64_t i(0);i<number_superbuckets;++i){
		*out_files[i]<<flush;
		delete(out_files[i]);
	}
}
