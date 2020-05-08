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


string version("v1.0.2");

char ch;
string query,fof(""), color_dump_file("reindeer_matrix"), color_load_file(""), output_bcalm("bcalm_out"),output_union_bcalm("bcalm_union_out"),output("output_reindeer"), output_index("index_out");
uint k(31), threads(1);
bool record_counts(true), quantize(false), do_query_on_disk(false), bcalm(false), do_Index(false), do_Query(false), PE(false), do_log(false);
uint threshold(40);


void PrintHelp()
{
    cout <<
            "******************* REINDEER "<< version << "**********************************\n"
            "******************* REad Index for abuNDancE quERy ******************\n\n"
            
            "                    INDEX BUILDING\n\n"

            "      * Mandatory parameters\n"
            "--index                 :     Indexing mode\n"
            "-f <file>               :     File of file for colors. Either:\n"
            "                                  i) you've already computed each DBG on your samples, in this case the fof is a list of unitig files\n"
            "                              OR ii) you need Bcalm to be run to obtain unitigs per sample, in this case use --bcalm option\n"
            "      * Optional parameters\n"
            "-k                      :     k-mer size (default 31)\n"
            "--nocount               :     Only compute presence/absence, not abundance\n"
            "--bcalm                 :     Launch bcalm on each single read dataset\n"
            "--paired-end            :     Index using paired-end files (provide pairs of files one after another in the fof). Works only with --bcalm.\n"
            "--disk-query            :     Index for on-disk query (default: in-memory). To be used for large indexes that won't fit in RAM.\n"
            "--quantization          :     Quantize the abundances in bins (to use only with --count).\n"
            "--log-count             :     Record the log of the counts, gives approximate counts that save space (to use only with --count).\n"
            "      * Output options\n"
            "-o <file>               :     Directory to write output files (default: output_reindeer)\n\n"

            "                    QUERY\n\n"

            "      * Mandatory parameters\n"
            "--query                 :     Query mode\n"
            "-l                      :     Reindeer index directory (should be output_reindeer if you've not used -o during indexing)\n"
            "-q <FASTA>              :     FASTA query file with query sequences\n"
            "      * Optional parameters\n"
            "-P                      :     Threshold: at least P% of the positions in the query must be covered by k-mers present in a dataset for the dataset to be reported (default: " << threshold << "%)\n"
            "-o <file>               :     Directory to write output files (default: output_reindeer/)\n"
            "--disk-query            :     On-disk query (default: in-memory). To be used for large indexes that won't fit in RAM, if the index was constructed with the same option.\n\n"


            "                    General\n\n"
            "-t <integer>            :     Number of threads (default 1)\n"
            "--help                  :     Show help\n";
    exit(1);
}

void ProcessArgs(int argc, char** argv)
{
    
    const char* const short_opts = "k:P:t:f:l:g:q:o:w:";
    const option long_opts[] = {
            {"index", no_argument, nullptr, 'i'},
            {"help", no_argument, nullptr, 'h'},
            {"nocount", no_argument, nullptr, 'c'},
            {"bcalm", no_argument, nullptr, 'b'},
            {"query", no_argument, nullptr, 'Q'},
            {"paired-end", no_argument, nullptr, 'p'},
            {"disk-query", no_argument, nullptr, 'd'},
            {"quantization", no_argument, nullptr, 'u'},
            {"log-count", no_argument, nullptr, 'L'},
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
            case 'Q':
                do_Query=true;
                break;
            case 'u':
                quantize=true;
                break;
            case 'L':
                do_log=true;
                break;
            case 'q':
                query=optarg;
                break;
            case 'k':
                k=stoi(optarg);
                break;
            case 't':
                threads=stoi(optarg);
                break;
            case 'c':
                record_counts=false;
                break;
            case 'p':
                PE=true;
                break;
            case 'P':
                threshold=stoi(optarg);
                break;
            case 'l':
                color_load_file=optarg;
                break;
            case 'o':
                output=optarg;
                break;
            case 'b':
                bcalm=true;
                break;
            case 'd':
                do_query_on_disk=true;
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
    cout << "############# REINDEER version " << version << " #############" << endl << "Command line was: ";
    for(int i = 0; i < argc; ++i)
        cout << argv[i] << ' ';
    cout << endl << endl;
    if (not dirExists(output)){
        systRet=system(("mkdir " + output).c_str());
    }
    if ( (do_Index and do_Query) or not(do_Index or do_Query) ){
        cout << "You must choose: either indexing (--index) or querying (--query)\n" << endl;
        PrintHelp();
        return 0;
    }
    if (do_Index) //indexing only
    {
        if ( fof.size() == 0 ){
            cout << "Please specify input datasets to index using using -f.\n" << endl;
            PrintHelp();
            return 0;
        }
        if (bcalm)
        {
            cout << "Computing De Bruijn graphs on each dataset using Bcalm2...\n\n" << endl;
            fof = bcalm_launcher_single(fof,  k,  threads, output, output_bcalm, PE);
        }
        
        if ( fof.empty() or k == 0 )
        {
            cout << "Missing argument " << endl;
            PrintHelp();
            return 0;
        }
        string cl("");
        bcalm_cleanup();
        cout << "Indexing k-mers...\n\n" << endl;
        color_dump_file = output + "/" + color_dump_file;
        reindeer_index(k, fof, color_dump_file, record_counts, output, cl, threads,  do_query_on_disk, quantize, do_log);
        cout << "INDEX BUILDING = THE END" <<endl;
    } else {
        if (color_load_file.empty())
        {
            cout << "Missing argument -l" << endl;
            PrintHelp();
            return 0;
        } else {
            cout << "Querying..." << endl;
            string output_query(output + "/query_results");
            if (not dirExists(output_query)){
                systRet=system(("mkdir " + output_query).c_str());
            }
            if (not(query.empty() or exists_test(query)))
            {   
                cout << "[ERROR] Invalid query file" << endl;
                return 0;
			}
            reindeer_query(color_load_file, output_query,  threshold,   query, threads,  do_query_on_disk);
        }
    }
    return 0;
}
