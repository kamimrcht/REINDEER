
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



kmer nuc2int(char c){
	return (c/2)%4;
}



kmer nuc2intrc(char c){
	return ((c/2)%4)^2;
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
	for(uint64_t i(0);i<str.size();i++){
		res<<=2;
		res+=(str[i]/2)%4;
	}
	return res;
}


  uint32_t revhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



  uint32_t unrevhash ( uint32_t x ) {
				x+=1234567890;
	x = ( ( x >> 16 ) ^ x ) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ( ( x >> 16 ) ^ x ) * 0x64ea2d65;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



  uint64_t revhash ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



  uint64_t unrevhash ( uint64_t x ) {
				x+=1234567890;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



vector<bool> str2boolv(const string& str){
	vector<bool> res;
	for(uint64_t i(0);i<str.size();++i){
		if(str[i]=='G' or str[i]=='T'){
			res.push_back(true);
		}else{
			res.push_back(false);
		}
		if(str[i]=='C' or str[i]=='G'){
			res.push_back(true);
		}else{
			res.push_back(false);
		}
	}
	return res;
}


string bool2strv(const vector<bool>& v){
	string res;
	for(uint64_t i(0);i<v.size();i+=2){
		if(v[i]){
			if(v[i+1]){
				res+='G';
			}else{
				res+='T';
			}
		}else{
			if(v[i+1]){
				res+='C';
			}else{
				res+='A';
			}
		}
	}
	return res;
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



void cat_stream( istream& is,  ostream& os)
{
    const  streamsize buff_size = 1 << 16;
    char * buff = new char [buff_size];
    while (true)
    {
        is.read(buff, buff_size);
         streamsize cnt = is.gcount();
        if (cnt == 0) break;
        os.write(buff, cnt);
    }
    delete [] buff;
} // cat_stream



void decompress_file(const string& file, const string& output_file)
{
    //
    // Set up sink ostream
    //
     unique_ptr<  ofstream > ofs_p;
     ostream * os_p = & cout;
    if (not output_file.empty())
    {
        ofs_p =  unique_ptr<  ofstream >(new strict_fstream::ofstream(output_file));
        os_p = ofs_p.get();
    }
    //
    // Process files
    //

	unique_ptr<  istream > is_p(new zstr::ifstream(file));
	cat_stream(*is_p, *os_p);
} // decompress_files



// It's quite complex to bitshift mmx register without an immediate (constant) count
// See: https://stackoverflow.com/questions/34478328/the-best-way-to-shift-a-m128i
  __m128i mm_bitshift_left(__m128i x, unsigned count)
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



  __m128i mm_bitshift_right(__m128i x, unsigned count)
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



  __uint128_t rcb(const __uint128_t& in, uint64_t n){

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



bool exists_test (const  string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}



uint64_t rcbc(uint64_t in, uint64_t n){
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



kmer min_k (const kmer& k1,const kmer& k2){
	if(k1<=k2){
		return k1;
	}
	return k2;
}


uint16_t parseCoverage_exact(const string& str){
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
	// cout<<"exact"<< (uint16_t)stof(str.substr(pos+5,i))<<endl;
	return (uint16_t)stof(str.substr(pos+5,i));
}



uint16_t parseCoverage_bool(const string& str){
	return 1;
}





uint64_t asm_log2(const uint64_t x) {
  uint64_t y;
  asm ( "\tbsr %1, %0\n"
      : "=r"(y)
      : "r" (x)
  );
  return y;
}


uint64_t mylog2 (uint64_t val) {
    if (val == 0) return 0;
    if (val == 1) return 0;
    uint64_t ret = 0;
    while (val > 1) {
        val >>= 1;
        ret++;
    }
    return ret;
}


uint16_t parseCoverage_log2(const string& str){
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
	// cout<<"log"<<asm_log2((uint16_t)stof(str.substr(pos+5,i)))<<endl;
	return asm_log2((uint16_t)stof(str.substr(pos+5,i)));
}



bool kmer_in_superkmer(const kmer canon,const vector<kmer>& V){
	for(uint64_t i(0);i<V.size();i++){
		if(canon==V[i]){return true;}
	}
	return false;
}




void dump_vector_bool(const vector<bool>& V, ostream* out ){
	int cmp = 0;
	uint8_t output=0;
	vector<uint8_t> buf;
	for(uint64_t i(0);i<V.size();++i){
		output=output | ((V[i] ? 1 : 0) << cmp);
		cmp++;
		if(cmp==8){
			buf.push_back(output);
			if(buf.size()>=8000){
				out->write((char*)buf.data(),buf.size());
				//~ *out<<flush;
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





//~ vector<string> split(const string &s, char delim){
	//~ stringstream ss(s);
	//~ string item;
	//~ vector<string> elems;
	//~ while (getline(ss, item, delim)) {
		//~ elems.push_back(move(item));
	//~ }
	//~ return elems;
//~ }

vector<string> split(const string &s, char delim){
	vector<string> res;
	uint pred(0);
	for(uint i(0);i<s.size();++i){
		if(s[i]==delim){
			res.push_back(s.substr(pred,i-pred));
			pred=i+1;
		}
	}
	res.push_back(s.substr(pred));
	return res;
}


void split(const string &s, char delim,vector<string>& res){
	res.clear();
	uint pred(0);
	for(uint i(0);i<s.size();++i){
		if(s[i]==delim){
			res.push_back(s.substr(pred,i-pred));
			pred=i+1;
		}
	}
	res.push_back(s.substr(pred));
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
