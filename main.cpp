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
string query,fof, color_dump_file("reindeer_matrix"), color_load_file(""), bgreat_paths_fof(""), output_bcalm("bcalm_out"),output_union_bcalm("bcalm_union_out"),output("output_reindeer"), output_index("index_out"), graph ;
uint k(31);
uint t(1);
bool record_counts(false);
bool record_reads(false);
bool bcalm(false), do_Index(false), do_Query(false);
uint threshold(40);

void PrintHelp()
{
    cout <<
			"\n******************* REINDEER **************************************\n"
			"******************* REad Index for abuNDancE quERy ******************\n"


			"\n                    INDEX BUILDING\n"
			"* Mandatory\n"
            "--index <file>          :     Indexing mode"
            "-f                      :     file of file for colors\n"
            "                              either you've already computed each DBG on your samples, in this case the fof is a list of unitig files\n"
            "                              OR you need Bcalm to be run tp obtain unitigs per sample, in this case use --bcalm option"
            "--dump <file>           :     choose a filename to write index on disk\n"
            "* Options\n"
            "-k                      :     k-mer size (default 31)\n"
            "--abundance             :     retain abundances instead of presence/absence\n"
            "--bcalm                 :     launch bcalm on each single read dataset\n\n"
            "-g                      :     provide union DBG of all datasets\n\n"
            "--output <file>         :     directory to write output files\n"


            "                    QUERY\n"
			"* Mandatory\n"
			"--query                 :     Query mode\n"
            "-l <file>               :     load index from file on disk\n"
            "* Options\n"
            "--abundance             :     query abundances (index must have been built with --abundance)\n"
            "-q <FASTA>              :     FASTA query file with query sequences\n\n\n"
            "--g                     :     provide union DBG of all datasets\n\n"
            "--output <file>         :     directory to write output files\n"

            "                    Performances\n"
            "-t <integer>            :     number of threads (default 1)\n\n\n"

            "--help                  :     Show help\n";
    exit(1);
}

void ProcessArgs(int argc, char** argv)
{
	
    const char* const short_opts = "k:S:t:f:l:g:q:";
    const option long_opts[] = {
            {"index", no_argument, nullptr, 'i'},
            {"dump", no_argument, nullptr, 'w'},
            //~ {"load", required_argument, nullptr, 'l'},
            {"help", no_argument, nullptr, 'h'},
            {"count", no_argument, nullptr, 'c'},
            {"bcalm", no_argument, nullptr, 'b'},
            {"query", no_argument, nullptr, 'Q'},
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
			case 'i':
				do_Index=true;
				break;
			case 'f':
				fof=optarg;
				break;
			case 'g':
				graph=optarg;
				break;
			case 't':
				t=stoi(optarg);
				break;
			case 'Q':
				do_Query=true;
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
			case 'b':
				bcalm=true;
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
	int systRet;
    ProcessArgs(argc, argv);
	if ( fof.empty() or k == 0 )
	{
		PrintHelp();
		return 0;
	}

	if (not dirExists(output)){
		systRet=system(("mkdir " + output).c_str());
	//~ } else {
		//~ cout << "[WARNING] You should first remove current " << output << " directory or provide an output directory name using --output" << endl;
		//~ return 0;
	}
	if ( (do_Index and do_Query) or not(do_Index or do_Query) ){
		cout << "You must choose: either indexing or querying\n" << endl;
		PrintHelp();
		return 0;
	}
	if (bcalm)
	{
		cout << "Computing De Bruijn graphs on each dataset using Bcalm2...\n\n" << endl;
		fof = bcalm_launcher_single(fof,  k,  t, output, output_bcalm); // from here fof is a fof of unitig files
	} else {
		fof = getRealPaths(fof, output);
	}
	if (do_Index) //indexing only
	{
		if (graph.empty()){//todo try when fof is directly passed
			cout << "Computing Union De Bruijn graph using Bcalm2...\n\n" << endl;
			graph = bcalm_launcher_union(fof,  k,  t, output, output_union_bcalm);
		}
		string cl("");
		bcalm_cleanup();
		cout << "Indexing k-mers...\n\n" << endl;
		reindeer_index(k, graph, fof, color_dump_file, record_counts,record_reads, output, cl);
	} else {
		if (graph.empty()){
			cout << "Computing Union De Bruijn graph using Bcalm2...\n\n" << endl;
			graph = bcalm_launcher_union(fof,  k,  t, output, output_union_bcalm);
		}
		cout << "Querying..." << endl;
		string cd("");
		reindeer_query(k, graph, fof, cd, color_load_file,  record_counts,  record_reads,  threshold,  bgreat_paths_fof,  query);
		// todo interactive mode
	}
    return 0;
}
