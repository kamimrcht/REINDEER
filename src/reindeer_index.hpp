#include <string>

#ifndef REIN
#define REIN
using namespace std;

struct reindeer_index
{
	// MPHF options
	uint m2;
	uint c;
	uint bit;
	uint ex;
	uint m1;
	uint m3;
	
	uint k; //size of k-mers
	
	// options
	uint threads; // number of threads
	bool record_counts; // if true, record counts else presence/absence //todo replace by color mode
	uint color_mode;
	bool do_query_on_disk; // if true the full index is dumped on the disk, else it is rebuilt and in ram
	
	// other counts modes
	bool quantize; // record quantization
	bool do_log; // record log
	
	// color variables
	uint64_t nb_colors; // number of samples
	uint64_t nb_monotig;
	long nb_eq_class;
	
	//files
	string color_load_file;
	string color_dump_file;
	string fof;
	string output;
	string monotig_files;
	string dumped_index;
};

#endif
