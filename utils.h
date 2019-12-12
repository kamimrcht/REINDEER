


kmer nuc2int(char c);
kmer nuc2intrc(char c);
void read_vector_bool(vector<bool>& V, zstr::ifstream* out, uint64_t n_bits );
void dump_vector_bool(const vector<bool>& V, ostream* out );
string intToString(uint64_t n);
bool kmer_in_superkmer(const kmer canon,const vector<kmer>& V);
kmer min_k (const kmer& k1,const kmer& k2);
kmer str2num(const string& str);
uint64_t revhash ( uint64_t x );
uint64_t unrevhash ( uint64_t x );
uint16_t parseCoverage(const string& str);
string color_coverage2str(const vector<uint16_t>& V);
vector<string> split(const string &s, char delim);
kmer hash64shift(kmer key);
bool exists_test (const string& name);
__m128i mm_bitshift_right(__m128i x, unsigned count);
uint64_t rcbc(uint64_t in, uint64_t n);
string revComp(const string& s);
void decompress_file(const string& file, const string& output_file);
vector<bool> str2boolv(const string& str);
string bool2strv(const vector<bool>& v);
void split(const string &s, char delim,vector<string>& res);
template<typename T>
inline T xs(const T& x) { return unrevhash(x); }



template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '[';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}


// iterator from disk file of u_int64_t with buffered read,   todo template
class bfile_iterator : public std::iterator<std::forward_iterator_tag, kmer>{
	public:

	bfile_iterator()
	: _is(nullptr)
	, _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (kmer *) malloc(_buffsize*sizeof(kmer));
	}

	bfile_iterator(const bfile_iterator& cr)
	{
		_buffsize = cr._buffsize;
		_pos = cr._pos;
		_is = cr._is;
		_buffer = (kmer *) malloc(_buffsize*sizeof(kmer));
		 memcpy(_buffer,cr._buffer,_buffsize*sizeof(kmer) );
		_inbuff = cr._inbuff;
		_cptread = cr._cptread;
		_elem = cr._elem;
	}

	bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (kmer *) malloc(_buffsize*sizeof(kmer));
		fseek(_is,0,SEEK_SET);
		advance();
	}

	~bfile_iterator()
	{
		if(_buffer!=NULL)
			free(_buffer);
	}


	kmer const& operator*()  {  return _elem;  }

	bfile_iterator& operator++()
	{
		advance();
		return *this;
	}

	friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
	{
		if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
		assert(lhs._is == rhs._is);
		return rhs._pos == lhs._pos;
	}

	friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
	private:
	void advance()
	{
		_pos++;

		if(_cptread >= _inbuff)
		{
			int res = fread(_buffer,sizeof(kmer),_buffsize,_is);
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
	kmer _elem;
	FILE * _is;
	unsigned long _pos;

	kmer * _buffer; // for buffered read
	int _inbuff, _cptread;
	int _buffsize;
};


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

	bfile_iterator begin() const
	{
		return bfile_iterator(_is);
	}

	bfile_iterator end() const {return bfile_iterator(); }

	size_t        size () const  {  return 0;  }//todo ?

	private:
	FILE * _is;
};
