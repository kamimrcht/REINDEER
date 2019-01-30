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


	inline u_int64_t printPt( pthread_t pt) {
	  unsigned char *ptc = (unsigned char*)(void*)(&pt);
		u_int64_t res =0;
	  for (size_t i=0; i<sizeof(pt); i++) {
		  res+= (unsigned)(ptc[i]);
	  }
		return res;
	}


////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark utils
////////////////////////////////////////////////////////////////


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

		size_t        size () const  {  return 0;  }//todo ?

	private:
		FILE * _is;
	};

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark hasher
////////////////////////////////////////////////////////////////

	typedef std::array<uint64_t,10> hash_set_t;
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
		uint64_t next(hash_pair_t  & s ) {
			uint64_t s1 = s[ 0 ];
			const uint64_t s0 = s[ 1 ];
			s[ 0 ] = s0;
			s1 ^= s1 << 23; // a
			return ( s[ 1 ] = ( s1 ^ s0 ^ ( s1 >> 17 ) ^ ( s0 >> 26 ) ) ) + s0; // b, c
		}

        //this one returns all the  hashes
        hash_set_t operator ()  (const Item& key)
        {
			uint64_t s[ 2 ];

            hash_set_t   hset;

            hset[0] =  singleHasher (key, 0xAAAAAAAA55555555ULL);
            hset[1] =  singleHasher (key, 0x33333333CCCCCCCCULL);

            s[0] = hset[0];
            s[1] = hset[1];

            for(size_t ii=2;ii< 10 /* it's much better have a constant here, for inlining; this loop is super performance critical*/; ii++)
            {
                hset[ii] = next(s);
            }

            return hset;
        }
    private:
        SingleHasher_t singleHasher;
    };


////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark iterators
////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark BitVector
////////////////////////////////////////////////////////////////

	class bitVector {

	public:

		bitVector() : _nwords(0)
		{
			_bitArray = nullptr;
		}

		bitVector(uint64_t n) : _nwords(n >> 6)
		{
			n = (n + 63) / 64;
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
			_nwords  = (newsize + 63) / 64;
			_bitArray = (uint64_t *) realloc(_bitArray,_nwords*sizeof(uint64_t));
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
		uint64_t build_ranks(uint64_t offset =0)
		{
			_ranks.reserve(2+ size()/_nb_bits_per_rank_sample);

			uint64_t curent_rank = offset;
			for (size_t ii = 0; ii < _nwords; ii++) {
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

////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark level
////////////////////////////////////////////////////////////////



	class level{
	public:
		level(){ }

		~level() {
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

		uint64_t idx_begin;
		uint64_t hash_domain;
	};


////////////////////////////////////////////////////////////////
#pragma mark -
#pragma mark mphf
////////////////////////////////////////////////////////////////


    /* Hasher_t returns a single hash when operator()(elem_t key) is called.
       if used with XorshiftHashFunctors, it must have the following operator: operator()(elem_t key, uint64_t seed) */
    template <typename elem_t, typename Hasher_t, unsigned _nb_levels=16>
	class mphf {

        /* this mechanisms gets P hashes out of Hasher_t */
        typedef XorshiftHashFunctors<elem_t,Hasher_t> MultiHasher_t ;
       // typedef HashFunctors<elem_t> MultiHasher_t; // original code (but only works for int64 keys)  (seems to be as fast as the current xorshift)
		//typedef IndepHashFunctors<elem_t,Hasher_t> MultiHasher_t; //faster than xorshift


		struct buildState {
			std::vector< elem_t > setLevelFastmode;
			bitVector collisions_map;
			unsigned fastModeLevel;
			bool fast_mode;
		};

	public:
		mphf()
		{}


		~mphf()
		{

		}


		// allow perc_elem_loaded  elements to be loaded in ram for faster construction (default 3%), set to 0 to desactivate
		template <typename Range>
		mphf( size_t n, Range const& input_range, double gamma = 2.0, float perc_elem_loaded = 0.99) :
		_gamma(gamma), _nelem(n)
		{
			if(n ==0) return;

			buildState state = setup(perc_elem_loaded);

			for(unsigned ii = 0; ii< _nb_levels; ii++)
			{
				processLevel(input_range,ii, state);

				bitset.clearCollisions(_levels[ii].idx_begin , _levels[ii].hash_domain , state.collisions_map);
				state.collisions_map.clear();
			}

			_lastbitsetrank = bitset.build_ranks(0) ;

			//printf("used temp ram for construction : %lli MB \n",setLevelFastmode.capacity()* sizeof(elem_t) /1024ULL/1024ULL);
		}


		uint64_t lookup(elem_t elem)
		{


			hash_pair_t bbhash;  unsigned level;
			uint64_t level_hash = getLevel(bbhash,elem,&level);

			if( level == (_nb_levels-1))
			{
				auto in_final_map  = _final_hash.find (elem);
				if ( in_final_map == _final_hash.end() )
				{
					//elem was not in orignal set of keys
					return ULLONG_MAX; //  means elem not in set
				}
				else
				{
					uint64_t minimal_hp = in_final_map->second + _lastbitsetrank;
					//printf("lookup %llu  level %i   --> %llu \n",elem,level,minimal_hp);

					return minimal_hp;
				}
			}

			return bitset.rank(_levels[level].hash2idx(level_hash));
		}

		uint64_t nbKeys() const
		{
            return _nelem;
        }

		uint64_t totalBitSize()
		{
			// unordered map takes approx 42B per elem [personal test] (42B with uint64_t key, would be larger for other type of elem)
			return bitset.bitSize() + CHAR_BIT * (_final_hash.size()*42 + sizeof(mphf));
		}

		template <typename Range>  //typename Range,
        inline void inner_processLevel(const Range& range, unsigned i, buildState& state)
		{
			uint64_t hashidx = 0;
			uint64_t idxLevelsetLevelFastmode = 0;

			for(const elem_t& val: range)
			{
				hash_pair_t bbhash;  unsigned level;
				uint64_t level_hash;
				getLevel(bbhash,val,&level, i);

				if(level == i) //insert into lvl i
				{
					if(state.fast_mode && i == state.fastModeLevel)
					{

						uint64_t idxl2 = idxLevelsetLevelFastmode++;
						//si depasse taille attendue pour setLevelFastmode, fall back sur slow mode mais devrait pas arriver si hash ok et proba avec nous
						if(idxl2>= state.setLevelFastmode.size()) {
							state.fast_mode = false;
						} else {
							state.setLevelFastmode[idxl2] = val; // create set for fast mode
						}
					}

					//insert to level i+1 : either next level of the cascade or final hash if last level reached
					if(i == _nb_levels-1) //stop cascade here, insert into exact hash
					{
						// calc rank de fin  precedent level qq part, puis init hashidx avec ce rank, direct minimal, pas besoin inser ds bitset et rank
						_final_hash[val] = hashidx++;
					}
					else
					{

						if ( level == 0)
							level_hash = get_hasher().h0(bbhash,val);
						else if ( level == 1)
							level_hash = get_hasher().h1(bbhash,val);
						else
						{
							level_hash = get_hasher().next(bbhash);
						}
						insertIntoLevel(level_hash, i, state.collisions_map); //should be safe
					}
				}
			}

			if(state.fast_mode && i == state.fastModeLevel) //shrink to actual number of elements in set
			{
				//printf("\nresize setLevelFastmode to %lli \n",_idxLevelsetLevelFastmode);
				state.setLevelFastmode.resize(idxLevelsetLevelFastmode);
			}
		}


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


		private :
		// Collision probability for the first level
		double probaCollision() {
			return 1.0 -  pow(((_gamma*(double)_nelem -1 ) / (_gamma*(double)_nelem)),_nelem-1);
		}

		//Configure the levels and return the total size of the bitVector
		uint64_t configureLevels() {
			double proba_collision = probaCollision();

			//double sum_geom =_gamma * ( 1.0 +  proba_collision / (1.0 - proba_collision));
			//printf("proba collision %f  sum_geom  %f   \n",proba_collision,sum_geom);

			//build levels
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

		buildState setup(float percent_elem_loaded_for_fastMode)
		{
			bitset = bitVector(configureLevels());
			double proba_collision = probaCollision();
			unsigned fastModeLevel = 0;
			for(; fastModeLevel<_nb_levels ; fastModeLevel++)
			{
				 if((float)pow(proba_collision,fastModeLevel) < percent_elem_loaded_for_fastMode)
				 {
				 	break;
				 }
			}

			return buildState {
				.setLevelFastmode = std::vector<elem_t>(size_t(percent_elem_loaded_for_fastMode * (float)_nelem), elem_t{}),
				.collisions_map = { _levels[0].hash_domain },
				.fastModeLevel = fastModeLevel,
				.fast_mode = true
			};
		}


		//compute level and returns hash of last level reached
		uint64_t getLevel(hash_pair_t & bbhash, elem_t val,unsigned * res_level, unsigned maxlevel = 100, unsigned minlevel =0)
		//uint64_t getLevel(hash_pair_t & bbhash, elem_t val,int * res_level, int maxlevel = 100, int minlevel =0)
		{
			unsigned level = 0;
			uint64_t hash_raw=0;

			for (unsigned ii = 0; ii<(_nb_levels-1) &&  ii < maxlevel ; ii++ )
			{

				//calc le hash suivant
				 if ( ii == 0)
					hash_raw = get_hasher().h0(bbhash,val);
				else if ( ii == 1)
					hash_raw = get_hasher().h1(bbhash,val);
				else
				{
					hash_raw = get_hasher().next(bbhash);
				}

				if( ii >= minlevel && bitset.get(_levels[ii].hash2idx(hash_raw)) )
				{
					break;
				}

				level++;
			}

			*res_level = level;
			return hash_raw;
		}


		//insert into bitarray
		void insertIntoLevel(uint64_t level_hash, unsigned i, bitVector& collisions_map)
		{
			if(bitset.test_and_set(_levels[i].hash2idx(level_hash)))
			{
				collisions_map.test_and_set(_levels[i].hash2levelidx(level_hash));
			}
		}


		//loop to insert into level i
		template <typename Range>
		void processLevel(Range const& input_range, unsigned i, buildState& state)
		{
			if(state.fast_mode && i >= (state.fastModeLevel+1))
			{
				inner_processLevel(state.setLevelFastmode, i, state);
			}
			else
			{
				inner_processLevel(input_range, i, state);
			}
		}

	private:
		level _levels[_nb_levels];
		std::unordered_map<elem_t,uint64_t,Hasher_t> _final_hash;
		bitVector bitset;

		double _gamma;
		uint64_t _nelem;
		uint64_t _lastbitsetrank;

		MultiHasher_t get_hasher() { return {}; }
	};
}

