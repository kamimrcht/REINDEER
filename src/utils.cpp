#include "utils.hpp"

using namespace std;



uint64_t xorshift ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}




bool is_empty_file(ifstream& file)
{
	return file.peek() == ifstream::traits_type::eof();
}

int dirExists(string& path)
{
    struct stat info;

    int statRC = stat( path.c_str(), &info );
    if( statRC != 0 )
    {
        if (errno == ENOENT)  { return 0; } // something along the path does not exist
        if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
        return -1;
    }

    return ( info.st_mode & S_IFDIR ) ? 1 : 0;
}

//~ inline bool exists_test(const string& name) {
	//~ return ( access( name.c_str(), F_OK ) != -1 );
//~ }

vector<string> split_utils(const string &s, char delim){
	stringstream ss(s);
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)) {
		elems.push_back(move(item)); 
	}
	return elems;
}

///// parse bcalm headers //////////////
double parseCoverage_utils(const string& str){
	size_t pos(str.find("km:f:"));
	if(pos==string::npos)
	{
		pos=(str.find("KM:f:"));
	}
	if(pos==string::npos){
		return 1;
	}
	uint i(1);
	while(str[i+pos+5]!=' '){
		++i;
	}
	// WARNING WE RETURN COUNTS BETWEEN 0 AND 255
	float result(stof(str.substr(pos+5,i)));
	if (result > 255){
		return 255;
	}
	else
	{
		return result;
	}
}



uint32_t unitig_toui32(const string& s)
{
	uint32_t val(0);
	int32_t valI(stoi(s));
	val = (valI >= 0) ? (uint32_t) valI : (uint32_t) -valI;
	return val;
}



void parse_bgreat_output(string& input, vector<vector<uint64_t>>& unitigs_to_nodes)
{
	ifstream bgreat_file(input);
	string line;
	vector<string> unitigs_per_read;
	uint64_t readID(0);
	
	while(not bgreat_file.eof())
	{
		getline(bgreat_file,line);
		if(line.empty()){readID++; continue;}
		unitigs_per_read = split_utils(line,';');
		for (auto && u : unitigs_per_read)
		{ 
			uint32_t unitig(unitig_toui32(u) - 1);
			if (unitig >= unitigs_to_nodes.size())
			{
				unitigs_to_nodes.resize(unitig + 1);
			} 
			unitigs_to_nodes[unitig].push_back(readID);
		}
		readID ++;
	}
}

//~ uint16_t revhash ( uint32_t x ) {
	//~ x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	//~ x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	//~ x = ( ( x >> 16 ) ^ x ) * 0x64ea2d65;
	//~ x = ( ( x >> 16 ) ^ x );
	//~ return x;
//~ }

//~ vector<uint16_t> RLE16C(const vector<uint16_t>&V){
	//~ vector<uint16_t> res;
	//~ if(V.empty()){return res;}
	//~ uint16_t pred(V[0]);
	//~ uint64_t count(1);
	//~ for(uint64_t i(1);i<V.size();++i){
		//~ if(V[i]==pred){
			//~ count++;
		//~ }else{
			//~ res.push_back(pred);
			//~ res.push_back(count);
			//~ count=1;
			//~ pred=V[i];
		//~ }
	//~ }
	//~ res.push_back(pred);
	//~ res.push_back(count);
	//~ return res;
//~ }


//~ vector<uint8_t> RLE16C(const vector<uint16_t>&V){
vector<unsigned char> RLE16C(const vector<uint16_t>&V){
	vector<unsigned char> res;
	if(V.empty()){return res;}
	uint16_t pred(V[0]);

	uint8_t count(1);
	for(uint64_t i(1);i<V.size();++i){
		if(V[i]==pred){
			count++;
			if(count==254){
				res.push_back(min(pred/256,254));
				res.push_back(min(pred%256,254));
				res.push_back(count);
				count=1;
			}
		}else{
			res.push_back(min(pred/256,254));
			res.push_back(min(pred%256,254));
			res.push_back(count);
			count=1;
			pred=V[i];
		}
	}
	res.push_back(min(pred/256,254));
	res.push_back(min(pred%256,254));
	res.push_back(count);
	//~ vector<unsigned char> res_ca;
	//~ transform(begin(res), end(res), begin(res_ca), [](uint8_t i) { return '0' + i; });

	return res;
}



vector<uint16_t> RLE16D(const vector<uint8_t>&V){
	vector<uint16_t> res;
	if(V.size()<2){return res;}
	for(uint64_t i(0);i<V.size();i+=3){
		res.resize(V[i+2],V[i]*256+V[i+1]);
	}
	for (auto && k : res)
	cout << k <<  " ";
	cout << endl;
	return res;
}



vector<uint8_t> RLE8C(const vector<uint8_t>&V){
	vector<uint8_t> res;
	if(V.empty()){return res;}
	uint8_t pred(V[0]);
	uint8_t count(1);
	for(uint64_t i(1);i<V.size();++i){
		if(V[i]==pred){
			count++;
			if(count==254){
			res.push_back(pred);
			res.push_back(count);
			count=1;
			}
		}else{
			res.push_back(pred);
			res.push_back(count);
			count=1;
			pred=V[i];
		}
	}
	res.push_back(pred);
	res.push_back(count);
	return res;
}

vector<uint8_t> RLE8D(const vector<uint8_t>&V){
	vector<uint8_t> res;
	if(V.size()<2){return res;}
	for(uint64_t i(0);i<V.size();i+=2){
		res.resize(V[i+1],V[i]);
	}
	return res;
}





void new_paired_end_file(string& input, string& input2, string& output_file, bool fastq)
{
	string line,head,head2, sequence,junk;
	istream* in;
	istream* in2;
	in=new zstr::ifstream(input);
	in2=new zstr::ifstream(input2);
	ofstream out(output_file);
	vector<uint> lengths;
	while(not in->eof()){
		getline(*in,head);
		getline(*in2,head2);
		getline(*in,sequence);
		if(fastq)
		{
			getline(*in,junk);
			getline(*in,junk);
			head[0] = '>';
		}
		out << head << "\n";
		transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		out << sequence << "\n";
		sequence="";
		getline(*in2,sequence);
		if(fastq)
		{
			getline(*in2,junk);
			getline(*in2,junk);
			head2[0] = '>';
		}
		out << head2 << "\n";
		transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		out << sequence << "\n";
		sequence="";
	}
	delete(in); delete(in2);
}


void interleave_paired_end(string& fof, string& output)
{
	bool tested(false), fastq(true);
	ifstream file(fof);
	string new_fof_name(output + "/pe_fof.lst");
	string sample, sample2, header, output_name;
	uint file_index(0);
	ofstream new_fof(new_fof_name);
	while(not file.eof())
	{
		getline(file, sample);
		if (sample.empty()){break;}
		getline(file, sample2);
		if (sample.empty()){break;}
		if (not tested)
		{
			ifstream samp(sample);
			while(not file.eof())
			{
				getline(samp, header);
				if (header.empty()){break;}
				if (header[0] == '>')
					fastq = false;
				tested = true;
				break;
			}
		}
		output_name = output + "/PE_" + to_string(file_index);
		new_paired_end_file(sample, sample2, output_name, fastq);
		new_fof << output_name << endl;
		++file_index;
	}
	fof = new_fof_name;

}

uint64_t harmonic_mean(vector<uint64_t>& counts)
{
	if (counts.size() > 0)
	{
		float harmonicMean(0);
		for (uint i(0); i < counts.size(); ++i)
		{
			if (counts[i] > 0)
			{
				harmonicMean += 1/(float)(counts[i]);
			}
		}
		if (harmonicMean > 0)
		{
			return (uint64_t)(counts.size()/harmonicMean);
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

string getRealPath(string file, string& dir){
	char *symlinkpath = &dir[0];
	char actualpath [PATH_MAX+1];
	char *ptr;
	ptr = realpath(symlinkpath, actualpath);
	string rp(ptr);
	return rp + "/" + file;
}

uint16_t get_color_number(string& fof)
{	
	uint color(0);
	string line;
	ifstream fof_file(fof);
	while (not fof_file.eof())
	{
		getline(fof_file, line);
		if (not line.empty())
		{
			color++;
		}
	}
	return color;
}
