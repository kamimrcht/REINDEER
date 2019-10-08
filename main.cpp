#include <getopt.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <mutex>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include <omp.h>
#include <sys/stat.h>
#include <unistd.h>
#include "src/reindeer.hpp"
#include "src/launch_bcalm.hpp"

using namespace std;
using namespace chrono;

char ch;
string query,fof, color_dump_file(""), color_load_file(""), bgreat_paths_fof(""), output_bcalm("bcalm_out"),output_union_bcalm("bcalm_union_out"),output(""), output_index("index_out") ;
uint k(31);
uint t(1);
bool record_counts(false);
bool record_reads(false);
bool bcalm(false);
uint threshold(30);

void PrintHelp()
{
    cout <<
			"\n******************* REINDEER **************************************\n"
			"******************* REad Index for abuNDancE quERy ******************\n"


			"\n                    INDEX BUILDING\n"
			"* Mandatory\n"
            "--index <file>          :     file of file for colors\n"
            "                              either you've already computed each DBG on your samples, in this case the fof is a list of unitig files\n"
            "                              OR you need Bcalm to be run tp obtain unitigs per sample, in this case use --bcalm option"
            "--dump <file>           :     write index on disk in file\n"
            "* Options\n"
            "-k                      :     k-mer size (default 31)\n"
            "--abundance             :     retain abundances instead of presence/absence\n"
            "--bcalm                 :     launch bcalm on each single read dataset\n\n"
            "--output <file>         :     directory to write output files\n"


            "                    QUERY\n"
			"* Mandatory\n"
            "--load <file>           :     load index from file on disk\n"
            "* Options\n"
            "--abundance             :     query abundances (index must have been built with --abundance)\n"
            "--query <FASTA>         :     FASTA query file with query sequences\n\n\n"
            "--output <file>         :     directory to write output files\n"

            "                    Performances\n"
            "-t <integer>            :     number of threads (default 1)\n\n\n"

            "--help                  :     Show help\n";
    exit(1);
}

void ProcessArgs(int argc, char** argv)
{
	
    const char* const short_opts = "k:S:t:";
    const option long_opts[] = {
            {"index", required_argument, nullptr, 'g'},
            {"dump", no_argument, nullptr, 'w'},
            {"load", required_argument, nullptr, 'l'},
            {"help", no_argument, nullptr, 'h'},
            {"count", no_argument, nullptr, 'c'},
            {"bcalm", no_argument, nullptr, 'b'},
            {"query", required_argument, nullptr, 'q'},
            {"output", required_argument, nullptr, 'o'},
            {nullptr, no_argument, nullptr, 0}
    };

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
			case 'g':
				fof=optarg;
				break;
			case 't':
				t=stoi(optarg);
				break;
			case 'q':
				query=optarg;
				break;
			case 'k':
				k=stoi(optarg);
				break;
			case 'c':
				record_counts=true;
				break;
			case 'S':
				threshold=stoi(optarg);
				break;
			case 'l':
				color_load_file=optarg;
				break;
			case 'w':
				color_dump_file=optarg;
				break;
			case 'o':
				output=optarg;
				break;
        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            PrintHelp();
            break;
        }
    }
 
}

int main(int argc, char **argv)
{
    ProcessArgs(argc, argv);
    cout << "" << endl;
    //~ system("myfile.sh"); // myfile.sh should be chmod +x
	if ( fof.empty() or k == 0 )
	{
		PrintHelp();
		return 0;
	}
	if ( not output.empty() )
	{
		//check if the dir exists, mk it
		output_bcalm = output + "/" + output_bcalm;
		output_union_bcalm = output + "/"  + output_union_bcalm;
		output_index = output + "/"  + output_index;
	}
	if ( (not query.empty()) and (not fof.empty()) ){
		cout << "You must choose: either indexing or querying\n" << endl;
		PrintHelp();
		return 0;
	}
	if (bcalm)
	{
		
	}
	if (not fof.empty()) //indexing only
	{
		bcalm_launcher_union( fof,  k,  t, output_bcalm);
		string cl("");
		reindeer_index( k, output_union_bcalm, output_bcalm, color_dump_file, record_counts,record_reads, cl);
	} else {
		//~ reindeer_query( k, output_bcalm, fof,  string& color_load_file, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query);
	}
    return 0;
}
