#ifndef UTIL
#define UTIL


#include <sys/types.h>
#include <sys/stat.h>


using namespace std;


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

inline bool exists_test(const string& name) {
	return ( access( name.c_str(), F_OK ) != -1 );
}

vector<string> split(const string &s, char delim){
	stringstream ss(s);
	string item;
	vector<string> elems;
	while (getline(ss, item, delim)) {
		elems.push_back(move(item)); 
	}
	return elems;
}

///// parse bcalm headers //////////////
double parseCoverage(const string& str){
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
		unitigs_per_read = split(line,';');
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
	cout << rp + "/" + file << endl;
	return rp + "/" + file;
}


#endif
