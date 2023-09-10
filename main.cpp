#include "src/launch_bcalm.hpp"
#include "src/reindeer.hpp"
#include "src/utils.hpp"
#include "version.h" // manage version from git
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <ctype.h>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <map>
#include <mutex>
#include <omp.h>
#include <pthread.h>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include <filesystem>

using namespace std;
using namespace chrono;

// MPHF options
uint m1(10);
uint m3(5);

char ch;
string query, fof(""), color_load_file(""), output_bcalm("bcalm_out"), output_union_bcalm("bcalm_union_out"), output("reindeer_index_files"), output_format("raw");
uint k(31), threads(1);
bool record_counts(true), quantize(false), do_query_on_disk(true), bcalm(false), do_Index(false), do_Query(false), PE(false), do_log(false), keep_tmp(false), output_monotigs(false);
uint threshold(40);

void PrintHelp()
{
    cout << "******************* REINDEER " << VERSION << "**********************************\n"
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
                                                          "--mem-query             :     Index for in-memory query (default: on-disk). To be used for very large query batches if the index fits in RAM.\n"
                                                          "--quantization          :     Quantize the abundances in bins (to use only with --count).\n"
                                                          "--log-count             :     Record the log of the counts, gives approximate counts that save space (to use only with --count).\n"
                                                          "      * Output options\n"
                                                          "-o <file>               :     Directory to write reindeer index files (default: reindeer_index_files)\n\n"
                                                          "      * Advanced parameters (we recommend not to change these values unless you are very aware of REINDEER's inner components)\n"
                                                          "--minimizer-size <integer>   :    MPHF option: minimizer size\n"
                                                          "--buckets <integer>          :    MPHF option: number of buckets (log)\n"
                                                          "--keep-tmp                   :    keep tmp files\n"
                                                          "--monotigs                   :    Save the monotigs into monotigs.tsv inside of the output directory\n"
                                                          "                    QUERY\n\n"

                                                          "      * Mandatory parameters\n"
                                                          "--query                 :     Query mode\n"
                                                          "-l                      :     Reindeer index directory (should be reindeer_index_files if you've not used -o during indexing)\n"
                                                          "-q <FASTA>              :     FASTA query file with query sequences\n"
                                                          "      * Optional parameters\n"
                                                          "-P                      :     Threshold: at least P% of the positions in the query must be covered by k-mers present in a dataset for the dataset to be reported (default: "
         << threshold << "%)\n"
                         "-o <dir|file>           :     Directory to write query output files (default: query_results/out_query_Reindeer), or a file\n"
                         "--disk-query            :     On-disk query (default: in-memory). To be used for large indexes that won't fit in RAM, if the index was constructed with the same option.\n\n"
                         "--format, -F <format>   :     Choose output format : raw | <sum> <s> | <average> <a> | <normalize> <n> (default : raw)\n"

                         "                    General\n\n"
                         "-t <integer>            :     Number of threads (default 1)\n"
                         "--help, -h              :     Show help\n"
                         "--version, -V           :     Show version\n";
    exit(1);
}

void ProcessArgs(int argc, char** argv)
{

    const char* const short_opts = "pk:P:t:f:l:g:q:o:w:F:rVh";
    const option long_opts[] = {
        { "index", no_argument, nullptr, 'i' },
        { "help", no_argument, nullptr, 'h' },
        { "nocount", no_argument, nullptr, 'c' },
        { "bcalm", no_argument, nullptr, 'b' },
        { "query", no_argument, nullptr, 'Q' },
        { "paired-end", no_argument, nullptr, 'p' },
        { "disk-query", no_argument, nullptr, 'd' },
        { "quantization", no_argument, nullptr, 'u' },
        { "log-count", no_argument, nullptr, 'L' },
        { "minimizer-size", required_argument, nullptr, 'm' },
        { "buckets", required_argument, nullptr, 'n' },
        { "keep-tmp", no_argument, nullptr, 'r' },
        { "version", no_argument, nullptr, 'V' },
        { "format", required_argument, nullptr, 'F' },
        { "monotigs", no_argument, nullptr, 'M' },
        { nullptr, no_argument, nullptr, 0 }
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt) {
        case 'i':
            do_Index = true;
            break;
        case 'f':
            fof = optarg;
            break;
        case 'Q':
            do_Query = true;
            break;
        case 'u':
            quantize = true;
            break;
        case 'L':
            do_log = true;
            break;
        case 'q':
            query = optarg;
            break;
        case 'm':
            m1 = stoi(optarg);
            break;
        case 'n':
            m3 = stoi(optarg);
            break;
        case 'k':
            k = stoi(optarg);
            break;
        case 't':
            threads = stoi(optarg);
            break;
        case 'c':
            record_counts = false;
            break;
        case 'p':
            PE = true;
            break;
        case 'P':
            threshold = stoi(optarg);
            break;
        case 'r':
            keep_tmp = true;
            break;
        case 'l':
            color_load_file = optarg;
            break;
        case 'o':
            output = optarg;
            break;
        case 'b':
            bcalm = true;
            break;
        case 'd':
            do_query_on_disk = false;
            break;
        case 'V':
            cout << VERSION << endl;
            exit(0);
            break;
        case 'F':
            if (strcmp("s",optarg) == 0) {
                output_format = "sum";
            } else if (strcmp("a",optarg) == 0) {
                output_format = "average";
            } else if (strcmp("n",optarg) == 0) {
                output_format = "normalize";
            } else {
                output_format = optarg;
            }
            break;
        case 'M':
            output_monotigs = true;
            break;
        case 'h': // -h or --help
        case '?': // Unrecognized option
        default:
            PrintHelp();
            break;
        }
    }
}

int main(int argc, char** argv)
{
    int systRet;
    string rno(get_run_tag());
    ProcessArgs(argc, argv);
    cout << "############# REINDEER version " << VERSION << " #############" << endl
         << "Command line was: ";
    for (int i = 0; i < argc; ++i)
        cout << argv[i] << ' ';
    cout << endl
         << endl;

    if ((do_Index and do_Query) or not(do_Index or do_Query)) {
        cout << strerror(E2BIG) << " | " << "You must choose: either indexing (--index) or querying (--query)\n"
             << endl;
        PrintHelp();
        return E2BIG;
    }
    if (do_Index) //indexing only
    {
        if (fof.size() == 0) {
            cout << strerror(EINVAL) << " | " << "Please specify input datasets to index using using -f.\n"
                 << endl;
            PrintHelp();
            return EINVAL;
        }
        string reindeer_index_files;
        if (output == "reindeer_index_files")
            reindeer_index_files = output + "_" + rno;
        else
            reindeer_index_files = output;
        if (!filesystem::exists(reindeer_index_files)){
            try {
                filesystem::create_directories(reindeer_index_files);
            } catch (const exception& e) {
                return EINVAL;
            }
        }
        if (bcalm) {
            cout << "Computing De Bruijn graphs on each dataset using Bcalm2...\n\n"
                 << endl;
            fof = bcalm_launcher_single(fof, k, threads, reindeer_index_files, output_bcalm, PE);
        }

        if (fof.empty() or k == 0) {
            cout << strerror(EINVAL) << " | " << "Missing argument " << endl;
            PrintHelp();
            return EINVAL;
        }
        bcalm_cleanup();
        cout << "Indexing k-mers...\n\n"
             << endl;
        Reindeer_Index<uint16_t> reindeer_index(k, fof, record_counts, reindeer_index_files, threads, do_query_on_disk, quantize, do_log, m1, m3, !(keep_tmp), output_monotigs);
        //~ reindeer_index(k, fof, color_dump_file, record_counts, output, cl, threads,  do_query_on_disk, quantize, do_log, m1, m3);
        cout << "Reindeer index files written in " << reindeer_index_files << endl;
        cout << "INDEX BUILDING = THE END." << endl;
    } else {
        if (color_load_file.empty()) {
            cout << strerror(EINVAL) << " | " << "Missing argument -l" << endl;
            PrintHelp();
            return EINVAL;
        } else {
            cout << "Querying..." << endl;
            // Test if we have write access to output directory, and avoid loading index
            // define output filename later (in ::query_by_file)
            filesystem::path testPath("");
            if (output == "reindeer_index_files") {
                // default value, no output arg given, test current path to write
                // save in output, the final directory
                testPath = filesystem::current_path();
                testPath /= "query_results";
                output = testPath.string();
            } else {
                // -o is given check if a directory or a file given
                filesystem::path outputPath(output);
                // test if has an extension = file, otherwise consider directory
                // can't test if exists, because the output file will exists after the query
                if (outputPath.extension() == "") {
                    // it's directory
                    if (! Is_Writeable(outputPath)) {
                        cerr << "Can't write in the given path: " << outputPath << ", change to a writeable path with -o option" << endl;
                        return EINVAL;
                    }
                    testPath = outputPath;
                    testPath /= "query_results";
                    output = testPath.string();
                } else {
                    // it's a file
                    testPath = outputPath.parent_path();
                    // stupid, return empty if current_path
                    if (testPath == "")
                      testPath = ".";
                }
            }
            cerr << "final output: " << output << " and testPath: " << testPath << endl;
            if (filesystem::exists(testPath)) {
                if (! Is_Writeable(testPath)) {
                    cerr << "Can't write in the given path: " << testPath << ", change to a writeable path with -o option" << endl;
                    return EINVAL;
                }
            } else {
                // try to create it
                try {
                    filesystem::create_directories(testPath);
                } catch (const exception& e) {
                    cerr << "error eival, createdir: final output: " << output << " and testPath: " << testPath << endl;
                    return EINVAL;
                }
            }

            if (not(query.empty() or exists_test(query))) {
                cout << "[ERROR] " << strerror(EINVAL) << " | " << "Invalid query file" << endl;
                return EINVAL;
            }

            if (dirExists(color_load_file)) {
                if (not(exists_test(color_load_file + "/reindeer_matrix_eqc") and exists_test(color_load_file + "/reindeer_matrix_eqc_info") and exists_test(color_load_file + "/reindeer_index.gz") and exists_test(color_load_file + "/reindeer_matrix_eqc_position"))) {
                    cerr << "[ERROR] " << strerror(ENOENT) << " | " << "REINDEER index files are missing. Stopped." << endl;
                    return ENOENT;
                }
            } else {
                cerr << "[ERROR] " << strerror(ENOENT) << " | " << "REINDEER index directory is missing (or path in -l is wrong). Stopped." << endl;
                return ENOENT;
            }
            Reindeer_Index<uint16_t> reindeer_index(color_load_file, output, threads, !(keep_tmp));
            reindeer_index.load_index();
            reindeer_index.querying(query, threshold, output_format);
        }
    }
    return 0;
}
