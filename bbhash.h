// BooPHF library
// intended to be a minimal perfect hash function with fast and low memory construction, at the cost of (slightly) higher bits/elem than other state of the art libraries once built.
// should work with arbitray large number of elements, based on a cascade of  "collision-free" bit arrays

#pragma once
#include <stdio.h>
#include <climits>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include <array>
#include <unordered_map>
#include <vector>
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

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark hasher
////////////////////////////////////////////////////////////////

	typedef std::array<uint64_t,2> hash_pair_t;

/* alternative hash functor based on xorshift, taking a single hash functor as input.
we need this 2-functors scheme because HashFunctors won't work with unordered_map.
(rayan)
*/

	// wrapper around HashFunctors to return only one value instead of 7
	template <typename Item> class SingleHashFunctor
	{
	public:
		uint64_t operator ()  (const Item& key, uint64_t seed=0xAAAAAAAA55555555ULL) const  {  return hash64(key, seed);  }

	private:
		inline static uint64_t hash64 (uint64_t key, uint64_t seed) {
			return hash_bis(key, seed);
		}

		inline static uint64_t hash64 (unsigned __int128& key, uint64_t seed)
		{
			return hash_bis ((uint64_t) (key >> 64), seed^0xAAAAAAAA55555555ULL)  ^  hash_bis ((uint64_t)key, seed^0xBBBBBBBB66666666ULL);
		}

		inline static uint64_t hash_bis (uint64_t key, uint64_t seed)
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

	  //  uint64_t s[ 2 ];

		uint64_t next(uint64_t * s) {
			uint64_t s1 = s[ 0 ];
			const uint64_t s0 = s[ 1 ];
			s[ 0 ] = s0;
			s1 ^= s1 << 23; // a
			return ( s[ 1 ] = ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) ) + s0; // b, c
		}

		public:


		uint64_t h0(hash_pair_t  & s, const Item& key )
		{
			s[0] =  singleHasher (key, 0xAAAAAAAA55555555ULL);
			return s[0];
		}

		uint64_t h1(hash_pair_t  & s, const Item& key )
		{
			s[1] =  singleHasher (key, 0x33333333CCCCCCCCULL);
			return s[1];
		}


		//return next hash an update state s
		uint64_t next(hash_pair_t & s ) {
			uint64_t s1 = s[ 0 ];
			const uint64_t s0 = s[ 1 ];
			s[ 0 ] = s0;
			s1 ^= s1 << 23; // a
			return ( s[ 1 ] = ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) ) + s0; // b, c
		}

		uint64_t hot_fun iter(const Item& key, hash_pair_t& s, size_t i) {
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


	template <typename Iterator>
	struct iter_range
	{
		iter_range(Iterator b, Iterator e)
		: m_begin(b)
		, m_end(e)
		{}

		Iterator begin() const
		{ return m_begin; }

		Iterator end() const
		{ return m_end; }

		Iterator m_begin, m_end;
	};

	template <typename Iterator>
	iter_range<Iterator> range(Iterator begin, Iterator end)
	{
		return iter_range<Iterator>(begin, end);
	}


	class bitVector {

	public:

		bitVector() : _nwords(0)
		{
			_bitArray = nullptr;
		}

		bitVector(uint64_t n)
		{
			n = (n + 63) / 64;
			_nwords = n;
			_bitArray =  (uint64_t *) calloc (_nwords,sizeof(uint64_t));
		}

		~bitVector()
		{
			if(_bitArray != nullptr)
				free(_bitArray);
		}

		 //copy constructor
		 bitVector(bitVector const &r)
		 {
			 _nwords = r._nwords;
			 _ranks = r._ranks;
			 _bitArray = (uint64_t *) calloc (_nwords,sizeof(uint64_t));
			 memcpy(_bitArray, r._bitArray, _nwords*sizeof(uint64_t) );
		 }

		// Copy assignment operator
		bitVector &operator=(bitVector const &r)
		{
			if (&r != this)
			{
				_nwords = r._nwords;
				_ranks = r._ranks;
				if(_bitArray != nullptr)
					free(_bitArray);
				_bitArray = (uint64_t *) calloc (_nwords,sizeof(uint64_t));
				memcpy(_bitArray, r._bitArray, _nwords*sizeof(uint64_t) );
			}
			return *this;
		}

		// Move assignment operator
		bitVector &operator=(bitVector &&r)
		{
			//printf("bitVector move assignment \n");
			if (&r != this)
			{
				if(_bitArray != nullptr)
					free(_bitArray);

				_nwords = std::move (r._nwords);
				_ranks = std::move (r._ranks);
				_bitArray = r._bitArray;
				r._bitArray = nullptr;
			}
			return *this;
		}
		// Move constructor
		bitVector(bitVector &&r) : _bitArray ( nullptr)
		{
			*this = std::move(r);
		}


		void resize(uint64_t newsize)
		{
			uint64_t new_nwords  = (newsize + 63) / 64;
			if(new_nwords > _nwords) {
				_bitArray = (uint64_t *) realloc(_bitArray,_nwords*sizeof(uint64_t));
			}
			_nwords = new_nwords;
		}

		size_t size() const
		{
			return _nwords * sizeof(uint64_t) * CHAR_BIT;
		}

		uint64_t bitSize() const {return (_nwords*64ULL + _ranks.capacity()*64ULL );}

		//clear whole array
		void clear()
		{
			memset(_bitArray,0,_nwords*sizeof(uint64_t));
		}

		//clear collisions in interval, only works with start and size multiple of 64
		void clearCollisions(uint64_t start, size_t size, const bitVector& cc)
		{
			assume( (start & 63) ==0, "start must be a multple of 64");
			assume( (size & 63) ==0, "siaz must be a multple of 64");
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
			assume( (start & 63) ==0, "start must be a multple of 64");
			assume( (size & 63) ==0, "siaz must be a multple of 64");
			memset(_bitArray + (start/64ULL),0,(size/64ULL)*sizeof(uint64_t));
		}

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

			printf("rank array : size %lu \n",_ranks.size());
			for (uint64_t ii = 0; ii< _ranks.size(); ii++)
			{
				printf("%llu :  %lli,  ",ii,_ranks[ii]);
			}
			printf("\n");
		}

		//return value at pos
		uint64_t operator[](uint64_t pos) const
		{
			//unsigned char * _bitArray8 = (unsigned char *) _bitArray;
			//return (_bitArray8[pos >> 3ULL] >> (pos & 7 ) ) & 1;
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
			_bitArray [pos >> 6] |= (uint64_t(1) << (pos & 63) ) ;
		}

		//set bit pos to 0
		void reset(uint64_t pos)
		{
			assume(pos < size(), "pos=%llu > size()=%llu", pos,  size());
			_bitArray [pos >> 6] &= ~(uint64_t(1) << (pos & 63) ) ;
		}

		//return value of  last rank
		// add offset to  all ranks  computed
		uint64_t build_ranks(uint64_t upto, uint64_t offset =0)
		{
			assume(upto <= size(), "build_ranks: upto=%llu > size()=%llu", upto, size());
			const uint64_t max_idx = (upto + 63)/64;
			assume(max_idx <= _nwords, "build_ranks: max_idx=%llu > _nwords=%llu", max_idx, _nwords);
			_ranks.reserve(2+ upto/_nb_bits_per_rank_sample);

			uint64_t curent_rank = offset;
			for (size_t ii = 0; ii < max_idx; ii++) {
				if (((ii*64)  % _nb_bits_per_rank_sample) == 0) {
					_ranks.push_back(curent_rank);
				}
				curent_rank += __builtin_popcountll(_bitArray[ii]);
			}

			return curent_rank;
		}

		uint64_t rank(uint64_t pos) const
		{
			assume(pos < size(), "pos=%llu > size()=%llu", pos,  size());
			uint64_t word_idx = pos / 64ULL;
			uint64_t word_offset = pos % 64;
			uint64_t block = pos / _nb_bits_per_rank_sample;
			uint64_t r = _ranks[block];
			for (uint64_t w = block * _nb_bits_per_rank_sample / 64; w < word_idx; ++w) {
				r += __builtin_popcountll( _bitArray[w] );
			}
			uint64_t mask = (uint64_t(1) << word_offset ) - 1;
			r += __builtin_popcountll( _bitArray[word_idx] & mask);

			return r;
		}



		void save(std::ostream& os) const
		{
			os.write(reinterpret_cast<char const*>(&_nwords), sizeof(_nwords));
			os.write(reinterpret_cast<char const*>(_bitArray), (std::streamsize)(sizeof(uint64_t) * _nwords));
			size_t sizer = _ranks.size();
			os.write(reinterpret_cast<char const*>(&sizer),  sizeof(size_t));
			os.write(reinterpret_cast<char const*>(_ranks.data()), (std::streamsize)(sizeof(_ranks[0]) * _ranks.size()));
		}

		void load(std::istream& is)
		{
			is.read(reinterpret_cast<char*>(&_nwords), sizeof(_nwords));
			this->resize(_nwords << 6);
			is.read(reinterpret_cast<char *>(_bitArray), (std::streamsize)(sizeof(uint64_t) * _nwords));

			size_t sizer;
			is.read(reinterpret_cast<char *>(&sizer),  sizeof(size_t));
			_ranks.resize(sizer);
			is.read(reinterpret_cast<char*>(_ranks.data()), (std::streamsize)(sizeof(_ranks[0]) * _ranks.size()));
		}


	protected:
		uint64_t*  _bitArray;
		uint64_t _nwords;

		 // epsilon =  64 / _nb_bits_per_rank_sample   bits
		// additional size for rank is epsilon * size()
		static const uint64_t _nb_bits_per_rank_sample = 512; //512 seems ok
		std::vector<uint64_t> _ranks;
	};




	/* Hasher_t returns a single hash when operator()(elem_t key) is called.
	   if used with XorshiftHashFunctors, it must have the following operator: operator()(elem_t key, uint64_t seed) */
	template <typename elem_t, typename Hasher_t, unsigned _nb_levels=16>
	class mphf {
		/* this mechanisms gets P hashes out of Hasher_t */
		typedef XorshiftHashFunctors<elem_t,Hasher_t> MultiHasher_t ;

		struct buildState {
			bitVector accommodated;
			bitVector collisions_map;
		};


	struct levelIndexer {
		levelIndexer(){ }

		~levelIndexer() {
		}

		static uint64_t fastrange64(uint64_t word, uint64_t p) {
			//return word %  p;

			return (uint64_t)(((__uint128_t)word * (__uint128_t)p) >> 64);
		}

		uint64_t hash2levelidx(uint64_t hash_raw) const {
			return fastrange64( hash_raw, hash_domain);
		}


		uint64_t hash2idx(uint64_t hash_raw) const {
			return hash2levelidx(hash_raw) + idx_begin;
		}

		uint64_t endidx() const {
			return idx_begin + hash_domain;
		}

		uint64_t idx_begin;
		uint64_t hash_domain;
	};

	public:
		mphf()
		{}

		template <typename Range>
		mphf( size_t n, Range const& input_range, double gamma = 2.0) :
		_gamma(gamma), _nelem(n)
		{
			if(n == 0) return;

			bitset = bitVector(configureLevels());

			buildState state = {
				.accommodated = bitVector(_nelem),
				.collisions_map = bitVector(_levels[0].hash_domain),
			};

			unsigned level = 0;
			for(; ; level++) {
				if(level < _nb_levels) {
					if(processLevel(input_range, level, state)) {
						continue; // Some keys collided, needs another level
					} else {
						break;
					}
				} else {
					processFallbackLevel(input_range, state);
					level = _nb_levels - 1;
					break;
				}
			}
			_lastbitsetrank = bitset.build_ranks(_levels[level].endidx());
		}

		uint64_t lookup(elem_t elem) const
		{
			hash_pair_t bbhash;  unsigned level;
			uint64_t bit_idx = getLevel(elem, level, bbhash);

			if(level < _nb_levels) {
				return bitset.rank(bit_idx);
			} else {
				auto in_final_map = _final_hash.find(elem);
				if (in_final_map != _final_hash.end())
				{
					return in_final_map->second + _lastbitsetrank;
				} else {
					// elem was not in orignal set of keys
					assert(false, "not in final map");
					return ULLONG_MAX; // means elem not in set
				}
			}
			return bitset.rank(bit_idx);
		}

		uint64_t nbKeys() const
		{
			return _nelem;
		}

		uint64_t totalBitSize() const
		{
			// unordered map takes approx 42B per elem [personal test] (42B with uint64_t key, would be larger for other type of elem)
			return bitset.bitSize() + CHAR_BIT * (_final_hash.size()*42 + sizeof(mphf));
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

				_levels[ii].idx_begin = previous_idx;

				// round size to nearest superior multiple of 64, makes it easier to clear a level
				_levels[ii].hash_domain =  (( (uint64_t) (hash_domain * pow(proba_collision,ii)) + 63) / 64 ) * 64;
				if(_levels[ii].hash_domain == 0 ) _levels[ii].hash_domain  = 64 ;
				previous_idx += _levels[ii].hash_domain;

// 				printf("build level %i bit array : start %12llu, size %12llu  ",ii,_levels[ii].idx_begin,_levels[ii].hash_domain );
// 				printf(" expected elems : %.2f %% total \n",100.0*pow(proba_collision,ii));
			}

			return previous_idx;
		}

		//compute level and returns bit vector position in the last level reached
		uint64_t hot_fun getLevel(elem_t val, unsigned& level, hash_pair_t& bbhash, unsigned maxlevel = _nb_levels) const
		{
			uint64_t bit_idx = ULLONG_MAX;
			assume(maxlevel <= _nb_levels, "getLevel: maxlevel=%u > _nb_levels=%u", maxlevel, _nb_levels);
			for (level = 0; level < maxlevel ; level++ )
			{
				bit_idx = _levels[level].hash2idx(get_hasher().iter(val, bbhash, level));
				if(bitset.get(bit_idx))
					break;
				else
					continue;
			}
			return bit_idx;
		}

		template <typename Range>
		bool flatten_fun processLevel(const Range& range, unsigned level, buildState& state) {
			switch(level) {
				case 0: return processLevel_(range, 0, state);
				case 1: return processLevel_(range, 1, state);
				default: return processLevel_(range, level, state);
			}
		}

		template <typename Range>
		bool processLevel_(const Range& range, unsigned level, buildState& state)
		{
			const levelIndexer& level_indexer = _levels[level];
			if(level > 0)
				state.collisions_map.clear(0, level_indexer.hash_domain);

			bool did_collide = false;
			size_t key_idx = 0;
			for(auto it = std::begin(range) ; it != std::end(range) ; it++, key_idx++) {
				// After level 2 we can skip keys that were found at stable positions in level 1
				if(level >= 2 && state.accommodated.get(key_idx)) {
					// This branch is not useless, it allows to inform the branch predictor
					continue;
				} else {
					const elem_t& val = *it;
					hash_pair_t bbhash; unsigned found_level;
					getLevel(val, found_level, bbhash, level);
					if(found_level < level) {
						// Mark the key as accommodated in previous level, avoid rehasing it afterwards
						state.accommodated.set(key_idx);
					} else {
						// Not found in previous levels => insert into this one
						uint64_t level_hash = get_hasher().iter(val, bbhash, level);

						// Put the key in the bit vector, checking for collisions
						if(bitset.test_and_set(level_indexer.hash2idx(level_hash)))
						{
							state.collisions_map.test_and_set(level_indexer.hash2levelidx(level_hash));
							did_collide = true;
						}
					}
				}
			}

			if(did_collide)
				bitset.clearCollisions(level_indexer.idx_begin, level_indexer.hash_domain , state.collisions_map);

			return did_collide;
		}

		template <typename Range>
		inline void processFallbackLevel(const Range& range, buildState& state)
		{
			uint64_t hashidx = 0;
			size_t key_idx = 0;
			for(auto it = std::begin(range) ; it != std::end(range) ; it++, key_idx++) {
				if(state.accommodated.get(key_idx)) {
					// This branch is not useless, it allows to inform the branch predictor
					continue;
				} else {
					const elem_t& val = *it;
					hash_pair_t bbhash;  unsigned found_level;
					getLevel(val, found_level, bbhash);
					if(found_level >= _nb_levels) // Not in any level
					{
						_final_hash[val] = hashidx++;
					}
				}
			}
		}

	public:
		void save(std::ostream& os) const
		{
			os.write(reinterpret_cast<char const*>(&_gamma), sizeof(_gamma));
			os.write(reinterpret_cast<char const*>(&_lastbitsetrank), sizeof(_lastbitsetrank));
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
			is.read(reinterpret_cast<char*>(&_lastbitsetrank), sizeof(_lastbitsetrank));
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

	private:
		levelIndexer _levels[_nb_levels];
		std::unordered_map<elem_t,uint64_t,Hasher_t> _final_hash;
		bitVector bitset;

		double _gamma;
		uint64_t _nelem;
		uint64_t _lastbitsetrank;

		MultiHasher_t get_hasher() const { return {}; }
	};
}
