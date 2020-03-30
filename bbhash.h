// BooPHF library
// intended to be a minimal perfect hash function with fast and low memory construction, at the cost of (slightly) higher bits/elem than other state of the art libraries once built.
// should work with arbitray large number of elements, based on a cascade of  "collision-free" bit arrays

#pragma once
#include <stdio.h>
#include <climits>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <memory>
#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>
#include <string>
#include "common.h"
#include <sys/time.h>
#include <string.h>
#include <unistd.h>

namespace boomphf {

	// iterator from disk file of u_int64_t with buffered read,   todo template
	template <typename basetype>
	class bfile_iterator : public std::iterator<std::forward_iterator_tag, basetype>{
	public:

		bfile_iterator()
		: _is(nullptr)
		, _pos(0) ,_inbuff (0), _cptread(0)
		{
			_buffsize = 10000;
			_buffer = (basetype *) malloc(_buffsize*sizeof(basetype));
		}

		bfile_iterator(const bfile_iterator& cr)
		{
			_buffsize = cr._buffsize;
			_pos = cr._pos;
			_is = cr._is;
			_buffer = (basetype *) malloc(_buffsize*sizeof(basetype));
			 memcpy(_buffer,cr._buffer,_buffsize*sizeof(basetype) );
			_inbuff = cr._inbuff;
			_cptread = cr._cptread;
			_elem = cr._elem;
		}

		bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
		{
			//printf("bf it %p\n",_is);
			_buffsize = 10000;
			_buffer = (basetype *) malloc(_buffsize*sizeof(basetype));
			int reso = fseek(_is,0,SEEK_SET);
			advance();
		}

		~bfile_iterator()
		{
			if(_buffer!=NULL)
				free(_buffer);
		}


		basetype const& operator*()  {  return _elem;  }

		bfile_iterator& operator++()
		{
			advance();
			return *this;
		}

		friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
		{
			if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
			assert(lhs._is == rhs._is, "different input stream, not comparable");
			return rhs._pos == lhs._pos;
		}

		friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
	private:
		void advance()
		{

			//printf("_cptread %i _inbuff %i \n",_cptread,_inbuff);

			_pos++;

			if(_cptread >= _inbuff)
			{

				int res = fread(_buffer,sizeof(basetype),_buffsize,_is);

				//printf("read %i new elem last %llu  %p\n",res,_buffer[res-1],_is);
				_inbuff = res; _cptread = 0;

				if(res == 0)
				{
					_is = nullptr;
					_pos = 0;
					return;
				}
			}

			_elem = _buffer[_cptread];
			_cptread ++;
		}
		basetype _elem;
		FILE * _is;
		unsigned long _pos;

		basetype * _buffer; // for buffered read
		int _inbuff, _cptread;
		int _buffsize;
	};


	template <typename type_elem>
	class file_binary{
	public:

		file_binary(const char* filename)
		{
			_is = fopen(filename, "rb");

			if (!_is) {
				throw std::invalid_argument("Error opening " + std::string(filename));
			}
		}

		~file_binary()
		{
			fclose(_is);
		}

		bfile_iterator<type_elem> begin() const
		{
			return bfile_iterator<type_elem>(_is);
		}

		bfile_iterator<type_elem> end() const {return bfile_iterator<type_elem>(); }

		size_t		size () const  {  return 0;  }//todo ?

	private:
		FILE * _is;
	};


	typedef std::array<uint64_t,2> hash_pair_t;

/* alternative hash functor based on xorshift, taking a single hash functor as input.
we need this 2-functors scheme because HashFunctors won't work with unordered_map.
(rayan)
*/
	template <typename Item> class SingleHashFunctor
	{
	public:
		uint64_t operator ()  (const Item& key, uint64_t seed=0xAAAAAAAA55555555ULL) const  {  return hash64(key, seed);  }

	private:
		static uint64_t hash64 (uint64_t key, uint64_t seed) {
			return hash_bis(key, seed);
		}

		static uint64_t hash64 (unsigned __int128& key, uint64_t seed)
		{
			return hash_bis ((uint64_t) (key >> 64), seed^0xAAAAAAAA55555555ULL)  ^  hash_bis ((uint64_t)key, seed^0xBBBBBBBB66666666ULL);
		}

		static uint64_t hash_bis (uint64_t key, uint64_t seed)
		{
			uint64_t hash = seed;
			hash ^= (hash <<  7) ^  key * (hash >> 3) ^ (~((hash << 11) + (key ^ (hash >> 5))));
			hash = (~hash) + (hash << 21);
			hash = hash ^ (hash >> 24);
			hash = (hash + (hash << 3)) + (hash << 8);
			hash = hash ^ (hash >> 14);
			hash = (hash + (hash << 2)) + (hash << 4);
			hash = hash ^ (hash >> 28);
			hash = hash + (hash << 31);

			return hash;
		}
	};


	template <typename Item, class SingleHasher_t> class XorshiftHashFunctors
	{
		/*  Xorshift128*
			Written in 2014 by Sebastiano Vigna (vigna@acm.org)

			To the extent possible under law, the author has dedicated all copyright
			and related and neighboring rights to this software to the public domain
			worldwide. This software is distributed without any warranty.

			See <http://creativecommons.org/publicdomain/zero/1.0/>. */
		/* This is the fastest generator passing BigCrush without
		   systematic failures, but due to the relatively short period it is
		   acceptable only for applications with a mild amount of parallelism;
		   otherwise, use a xorshift1024* generator.

		   The state must be seeded so that it is not everywhere zero. If you have
		   a nonzero 64-bit seed, we suggest to pass it twice through
		   MurmurHash3's avalanching function. */

		uint64_t next(uint64_t * s) {
			uint64_t s1 = s[ 0 ];
			const uint64_t s0 = s[ 1 ];
			s[ 0 ] = s0;
			s1 ^= s1 << 23; // a
			return ( s[ 1 ] = ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) ) + s0; // b, c
		}

		public:


		uint64_t h0(hash_pair_t  & s, const Item& key ) const
		{
			s[0] =  singleHasher (key, 0xAAAAAAAA55555555ULL);
			return s[0];
		}

		uint64_t h1(hash_pair_t  & s, const Item& key ) const
		{
			s[1] =  singleHasher (key, 0x33333333CCCCCCCCULL);
			return s[1];
		}


		//return next hash an update state s
		uint64_t next(hash_pair_t & s ) const {
			uint64_t s1 = s[ 0 ];
			const uint64_t s0 = s[ 1 ];
			s[ 0 ] = s0;
			s1 ^= s1 << 23; // a
			return ( s[ 1 ] = ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) ) + s0; // b, c
		}

		uint64_t hot_fun iter(const Item& key, hash_pair_t& s, size_t i) const {
			uint64_t h;
			if(i == 0)
				h = h0(s, key);
			else if(i == 1)
				h = h1(s, key);
			else
				h = next(s);
			return h;
		}

	private:
		SingleHasher_t singleHasher;
	};


	using memreport_t = std::unordered_map<std::string, size_t>;
	inline void print_memreport(const memreport_t& report) {
		using pair_t = std::pair<std::string, size_t>;
		std::vector<pair_t> vec(report.begin(), report.end());
		std::sort(vec.begin(), vec.end());
		size_t max_l = 0;
		for(const pair_t& p : report) max_l = std::max(max_l, p.first.length());
		std::cout << "--------------------------------------------------------------------------------\nMemory usage:\n";
		uint64_t total = 0;
		for(const pair_t& p : vec) {
			total += p.second;
			float size = p.second;
			unsigned unit = 0;
			for(; size > 1024 ; unit++, size /= 1024);

			std::cout << p.first;
			for(size_t i=0 ; i < 1 + max_l - p.first.size() ; i++)
				std::cout << " ";
			std::cout << " : " << size<< "BKMGT"[unit] << std::endl;
		}

		float size = total;
		unsigned unit = 0;
		for(; size > 1024 ; unit++, size /= 1024);
		std::cout << "Total : " << size<< "BKMGT"[unit] << std::endl;
	}

	class bitVector {

	public:

		bitVector() : _nwords(0), _nranks(0)
		{
		}

		bitVector(uint64_t n)
		{
			n = (n + 63) / 64;
			_nwords = n;
			_bitArray = std::unique_ptr<uint64_t[]>(new uint64_t[_nwords]());
		}

		//copy constructor
		bitVector(bitVector const &r)
		{
		 *this = r;
		}

		// Copy assignment operator
		bitVector &operator=(bitVector const &r)
		{
			if (&r != this)
			{
				_nwords = r._nwords;
				_bitArray = std::unique_ptr<uint64_t[]>(new uint64_t[_nwords]);
				memcpy(_bitArray.get(), r._bitArray.get(), _nwords*sizeof(uint64_t) );

				_nranks = r._nranks;
				_ranks = std::unique_ptr<uint64_t[]>(new uint64_t[_nranks]);
				memcpy(_ranks.get(), r._ranks.get(), _nranks*sizeof(uint64_t) );
			}
			return *this;
		}

		// Move assignment operator
		bitVector &operator=(bitVector &&r)
		{
			_nwords = r._nwords;
			_nranks = r._nranks;
			_ranks = std::move(r._ranks);
			_bitArray = std::move(r._bitArray);
			r._bitArray = nullptr;
			r._nwords = 0;
			r._nranks = 0;
			return *this;
		}
		// Move constructor
		bitVector(bitVector &&r) : _bitArray ( nullptr)
		{
			*this = std::move(r);
		}

		size_t size() const
		{
			return _nwords * sizeof(uint64_t) * CHAR_BIT;
		}

		void report_memusage(memreport_t& report, const std::string& prefix="bitVector", bool add_struct=true) const {
			if(add_struct)
				report[prefix+"::sizeof(struct)"] += sizeof(bitVector);
			report[prefix+"::array"] += sizeof(uint64_t) * _nwords;
			report[prefix+"::ranks"] += sizeof(uint64_t) * _nranks;
		}

		//clear whole array
		void clear()
		{
			memset(_bitArray.get(),0,_nwords*sizeof(uint64_t));
		}

		//clear collisions in interval, only works with start and size multiple of 64
		void clearCollisions(uint64_t start, size_t size, const bitVector& cc)
		{
			assume( (start & 63) ==0, "start must be a multiple of 64");
			assume( (size & 63) ==0, "size must be a multiple of 64");
			uint64_t ids = (start/64ULL);
			assume( ids + (size/64ULL) <= _nwords, "clearCollisions called after start");
			for(uint64_t ii =0;  ii< (size/64ULL); ii++ )
			{
				_bitArray[ids+ii] =  _bitArray[ids+ii] & (~ (cc.get64(ii)) );
			}
		}

		//clear interval, only works with start and size multiple of 64
		void clear(uint64_t start, size_t size)
		{
			assume( (start & 63) ==0, "start must be a multiple of 64");
			assume( (size & 63) ==0, "size must be a multiple of 64");
			memset(&_bitArray[start/64ULL],0,(size/64ULL)*sizeof(uint64_t));
		}

		void clear(uint64_t start=0) { clear(start, _nwords); }

#ifndef NDEBUG
		//for debug purposes
		void print() const
		{
			printf("bit array of size %lli: \n",size());
			for(uint64_t ii = 0; ii< size(); ii++)
			{
				if(ii%10==0)
					printf(" (%llu) ",ii);
				int val = (_bitArray[ii >> 6] >> (ii & 63 ) ) & 1;
				printf("%i",val);
			}
			printf("\n");

			printf("rank array : size %lu \n",_nranks);
			for (uint64_t ii = 0; ii< _nranks; ii++)
			{
				printf("%llu :  %lli,  ",ii,_ranks[ii]);
			}
			printf("\n");
		}
#endif

		//return value at pos
		uint64_t operator[](uint64_t pos) const
		{
			assume(pos < size(), "pos=%llu > size()=%llu", pos,  size());
			return (_bitArray[pos >> 6ULL] >> (pos & 63 ) ) & 1;
		}

		//return old val and set to 1
		uint64_t test_and_set(uint64_t pos)
		{
			assume(pos < size(), "pos=%llu > size()=%llu", pos,  size());
			uint64_t& val = _bitArray[pos >> 6];
			uint64_t oldval = val;
			val = oldval | uint64_t(1) << (pos & 63);
			return  ( oldval >> (pos & 63 ) ) & 1;
		}

		uint64_t get(uint64_t pos) const
		{
			return (*this)[pos];
		}

		uint64_t get64(uint64_t pos_64) const
		{
			assume(pos_64 < size(), "pos_64=%llu > size()/64=%llu", pos_64, _nwords);
			return _bitArray[pos_64];
		}

		//set bit pos to 1
		void set(uint64_t pos)
		{
			assume(pos < size(), "pos=%llu > size()=%llu", pos,  size());
			_bitArray [pos >> 6] |= uint64_t(1) << (pos & 63) ;
		}

		//set bit pos to 0
		void reset(uint64_t pos)
		{
			assume(pos < size(), "pos=%llu > size()=%llu", pos,  size());
			_bitArray [pos >> 6] &= ~(uint64_t(1) << (pos & 63)) ;
		}

		//return value of  last rank
		// add offset to  all ranks  computed
		uint64_t build_ranks(uint64_t upto, uint64_t offset =0)
		{
			assume(upto <= size(), "build_ranks: upto=%llu > size()=%llu", upto, size());
			const uint64_t max_idx = (upto + 63)/64;
			assume(max_idx <= _nwords, "build_ranks: max_idx=%llu > _nwords=%llu", max_idx, _nwords);
			_nranks = (max_idx + (_words_per_rankingblock-1))/_words_per_rankingblock;
			_ranks = std::unique_ptr<uint64_t[]>(new uint64_t[_nranks]());

			uint64_t curent_rank = offset;
			for (size_t ii = 0, ridx = 0; ii < max_idx; ii++) {
				if (ii % _words_per_rankingblock == 0) {
					_ranks[ridx++] = curent_rank;
				}
				curent_rank += __builtin_popcountll(_bitArray[ii]);
				assert(ridx <=  _nranks, "%llu > %llu", ridx, _nranks);
			}

			return curent_rank;
		}

		uint64_t rank(uint64_t pos) const
		{
			assume(pos < size(), "pos=%llu > size()=%llu", pos,  size());
			uint64_t word_idx = pos / 64ULL;
			uint64_t word_offset = pos % 64;
			uint64_t block = word_idx / _words_per_rankingblock;
			uint64_t r = _ranks[block];
			for (uint64_t w = block * _words_per_rankingblock; w < word_idx; ++w) {
				r += __builtin_popcountll( _bitArray[w] );
			}
			uint64_t mask = (uint64_t(1) << word_offset ) - 1;
			r += __builtin_popcountll( _bitArray[word_idx] & mask);
			return r;
		}

		void save(std::ostream& os) const
		{
			os.write(reinterpret_cast<char const*>(&_nwords), sizeof(_nwords));
			os.write(reinterpret_cast<char const*>(_bitArray.get()), (std::streamsize)(sizeof(uint64_t) * _nwords));
			os.write(reinterpret_cast<char const*>(&_nranks),  sizeof(size_t));
			os.write(reinterpret_cast<char const*>(_ranks.get()), (std::streamsize)(sizeof(_ranks[0]) * _nranks));
		}

		void load(std::istream& is)
		{
			is.read(reinterpret_cast<char*>(&_nwords), sizeof(_nwords));
			_bitArray = std::unique_ptr<uint64_t[]>(new uint64_t[_nwords]);
			is.read(reinterpret_cast<char *>(_bitArray.get()), (std::streamsize)(sizeof(uint64_t) * _nwords));

			is.read(reinterpret_cast<char *>(&_nranks),  sizeof(size_t));
			_ranks = std::unique_ptr<uint64_t[]>(new uint64_t[_nranks]);
			is.read(reinterpret_cast<char*>(_ranks.get()), (std::streamsize)(sizeof(_ranks[0]) * _nranks));
		}

	protected:
		std::unique_ptr<uint64_t[]>  _bitArray;
		uint64_t _nwords;

		static const uint64_t _words_per_rankingblock = 16; // Have to popcount at most two cache lines
		//FIXME: ...if 8 words boundaries are aligned on cache lines, can we do something about that ?
		std::unique_ptr<uint64_t[]> _ranks;
		uint64_t _nranks;
	};




	/* Hasher_t returns a single hash when operator()(elem_t key) is called.
	   if used with XorshiftHashFunctors, it must have the following operator: operator()(elem_t key, uint64_t seed) */
	template <typename elem_t, typename Hasher_t, unsigned _nb_levels=16>
	class mphf : private XorshiftHashFunctors<elem_t,Hasher_t> { //Empty base opt
		/* this mechanisms gets P hashes out of Hasher_t */
		typedef XorshiftHashFunctors<elem_t,Hasher_t> MultiHasher_t ;

		struct buildState {
			bitVector accommodated;
			bitVector collisions_map;
			uint64_t level_offset;
		};

	public:
		mphf()
		{}

		template <typename Range>
		noinline_fun flatten_fun mphf( size_t n, Range const& input_range, double gamma = 2.0) :
		_gamma(gamma), _nelem(n)
		{
			if(n == 0) return;

			bitset = bitVector(configureLevels());

			buildState state = {
				.accommodated = bitVector(_nelem),
				.collisions_map = bitVector(_hash_domains[0]),
				.level_offset = 0
			};

			for(unsigned level = 0; ; level++) {
				if(level < _nb_levels) {
					if(processLevel(input_range, level, state)) {
						continue; // Some keys collided, needs another level
					} else {
						bitset.build_ranks(state.level_offset);
						break;
					}
				} else {
					uint64_t next_rank = bitset.build_ranks(bitset.size());
					processFallbackLevel(input_range, state, next_rank);
					break;
				}
			}
		}

		hot_fun uint64_t lookup(elem_t elem) const
		{
			const_level_finder finder(*this, elem);
			if(finder) {
				return finder.rank();
			} else {
				auto in_final_map = _final_hash.find(elem);
				if (in_final_map != _final_hash.end())
				{
					return in_final_map->second;
				} else {
					// elem was not in orignal set of keys
					//assert(false, "not in final map"); // Rayan: removed as we definitely want blight to report when a kmer isnt found; not throw an assert false
					return ULLONG_MAX; // means elem not in set
				}
			}
		}

		uint64_t nbKeys() const { return _nelem; }

		void report_memusage(memreport_t& report, const std::string& prefix="mphf", bool add_struct=true) const {
			if(add_struct)
				report[prefix+"::sizeof(struct)"] += sizeof(mphf);
			report[prefix+"::fallback_map"] += 24*_final_hash.size() + 8*_final_hash.bucket_count();
			bitset.report_memusage(report, prefix+"::bitvector", false);
		}

	private:

		//Configure the levels and return the total size of the bitVector
		uint64_t configureLevels() {
			double proba_collision = 1.0 -  pow(((_gamma*(double)_nelem -1 ) / (_gamma*(double)_nelem)),_nelem-1);

			//double sum_geom =_gamma * ( 1.0 +  proba_collision / (1.0 - proba_collision));
			//printf("proba collision %f  sum_geom  %f   \n",proba_collision,sum_geom);

			uint64_t previous_idx =0;
			size_t hash_domain = (size_t)  (ceil(double(_nelem) * _gamma)) ;
			for(unsigned ii=0; ii<_nb_levels; ii++)
			{

				//_levels[ii].idx_begin = previous_idx;

				// round size to nearest superior multiple of 64, makes it easier to clear a level
				_hash_domains[ii] =  (( (uint64_t) (hash_domain * pow(proba_collision,ii)) + 63) / 64 ) * 64;
				if(_hash_domains[ii] == 0 ) _hash_domains[ii]  = 64 ;
				previous_idx += _hash_domains[ii];

// 				printf("build level %i bit array : start %12llu, size %12llu  ",ii,_levels[ii].idx_begin,_levels[ii].hash_domain );
// 				printf(" expected elems : %.2f %% total \n",100.0*pow(proba_collision,ii));
			}

			return previous_idx;
		}


		template<typename Map> // either mphf or const mphf
		struct level_finder_tmp {
			hot_fun forceinline_fun level_finder_tmp(Map& map, const elem_t& key, uint64_t max_level=_nb_levels) : _map(map) {
				assume(max_level <= _nb_levels, "max_level=%u > _nb_levels=%u", max_level <= _nb_levels);
				hash_pair_t bbhash;
				uint64_t bit_begin = 0; // First bit of the current level in the bitvector
				for(unsigned level = 0; level < _nb_levels; level++) {
					uint64_t hash = _map.get_hasher().iter(key, bbhash, level);
					uint64_t hash_domain = _map._hash_domains[level];
					_level_bit = fastmod64(hash, hash_domain);
					_bit = _level_bit + bit_begin;


					if(level >= max_level) {
						break;
					} else if( _map.bitset[_bit]) {
						_found = true;
						break;
					}

					bit_begin += hash_domain;
				}
			}

			operator bool() { return _found; }

			bool set(bitVector& collisions) const {
				if(_map.bitset.test_and_set(_bit)) {
					collisions.test_and_set(_level_bit);
					return true;
				} else {
					return false;
				}
			}

			uint64_t rank() const { return _map.bitset.rank(_bit); }


		private:
			Map& _map;
			uint64_t _bit, _level_bit;
			bool _found = false;

			uint64_t fastmod64(uint64_t word, uint64_t p) {
				return (uint64_t)(((__uint128_t)word * (__uint128_t)p) >> 64);
			}
		};

		using level_finder = level_finder_tmp<mphf>;
		using const_level_finder = level_finder_tmp<const mphf>;

		template <typename Range>
		flatten_fun bool processLevel(const Range& range, unsigned level, buildState& state) {
			switch(level) {
				case 0: return processLevel_(range, 0, state);
				case 1: return processLevel_(range, 1, state);
				default: return processLevel_(range, level, state);
			}
		}

		template <typename Range>
		bool processLevel_(const Range& range, unsigned level, buildState& state)
		{
			assume(level < _nb_levels, "level=%u >= _nb_levels", level, _nb_levels);

			uint64_t hash_domain = _hash_domains[level];
			if(level > 0)
				state.collisions_map.clear(0, hash_domain);

			bool did_collide = false; // True if we have at least one collision
			size_t key_idx = 0;
			for(auto it = std::begin(range) ; it != std::end(range) ; ++it, key_idx++) {
				// After level 2 we can skip keys that were found at stable positions in level 1
				if(level >= 2 && state.accommodated.get(key_idx)) {
					// This branch is not useless, it allows to inform the branch predictor
					continue;
				} else {
					level_finder finder(*this, *it, level);
					if(finder) {
						// Mark the key as accommodated in previous level, avoiding rehasing it afterwards
						state.accommodated.set(key_idx);
					} else {
						did_collide |= finder.set(state.collisions_map);
					}
				}
			}

			if(did_collide)
				bitset.clearCollisions(state.level_offset, hash_domain , state.collisions_map);

			state.level_offset += hash_domain;
			return did_collide;
		}

		template <typename Range>
		cold_fun void processFallbackLevel(const Range& range, buildState& state, uint64_t next_rank)
		{
			size_t key_idx = 0;
			for(auto it = std::begin(range) ; it != std::end(range) ; ++it, key_idx++) {
				if(state.accommodated.get(key_idx)) {
					// This branch is not useless, it allows to inform the branch predictor
					continue;
				} else {
					const elem_t& val = *it;
					level_finder finder(*this, val);
					if(!finder) // Not in any level
					{
						_final_hash[val] = next_rank++;
					}
				}
			}
		}

	public:
		void save(std::ostream& os) const
		{
			os.write(reinterpret_cast<char const*>(&_gamma), sizeof(_gamma));
			os.write(reinterpret_cast<char const*>(&_nelem), sizeof(_nelem));
			bitset.save(os);

			//save final hash
			size_t final_hash_size = _final_hash.size();

			os.write(reinterpret_cast<char const*>(&final_hash_size), sizeof(size_t));

			// typename std::unordered_map<elem_t,uint64_t,Hasher_t>::iterator
			for (auto it = _final_hash.begin(); it != _final_hash.end(); ++it )
			{
				os.write(reinterpret_cast<char const*>(&(it->first)), sizeof(elem_t));
				os.write(reinterpret_cast<char const*>(&(it->second)), sizeof(uint64_t));
			}
		}

		void load(std::istream& is)
		{
			is.read(reinterpret_cast<char*>(&_gamma), sizeof(_gamma));
			is.read(reinterpret_cast<char*>(&_nelem), sizeof(_nelem));
			bitset.load(is);

			//mini setup, recompute size of each level
			configureLevels();

			//restore final hash
			_final_hash.clear();
			size_t final_hash_size ;

			is.read(reinterpret_cast<char *>(&final_hash_size), sizeof(size_t));

			for(unsigned int ii=0; ii<final_hash_size; ii++)
			{
				elem_t key;
				uint64_t value;

				is.read(reinterpret_cast<char *>(&key), sizeof(elem_t));
				is.read(reinterpret_cast<char *>(&value), sizeof(uint64_t));

				_final_hash[key] = value;
			}
		}

		const MultiHasher_t& get_hasher() const { return *this; }

	private:
		uint64_t _hash_domains[_nb_levels];
		std::unordered_map<elem_t,uint64_t,Hasher_t> _final_hash;
		bitVector bitset;

		double _gamma;
		uint64_t _nelem;
	};
}

