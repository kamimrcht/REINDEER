#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdint>
#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

using namespace std;
namespace fs = filesystem;

typedef struct info{
    uint64_t nb_monotig = 0;
    uint k = 0;
    uint record_option = 0;
    long nb_eq_class = 0;
    uint64_t nb_colors = 0;
    bool do_query_on_disk = 0;
    vector<pair<string,uint64_t>> kbf;
} Info;

void help() {
    cout << endl;
    cout <<
    "* Mandatory parameters *\n"
    "-f, --fof <file>                   :       fof (File of files)\n"
    "-i, --info <directory or file>     :       Either the directory where the index is stored or the reindeer_maxtrix_eqc_info file\n\n"
    "* General options *\n"
    "-h, --help                         :       Print help (what you are seeing)\n"
    << endl;
    exit(0);
}

vector<pair<string,string>> process_args(int argc, char** argv) {

    vector<pair<string,string>> paths;

    const char* const short_opts = "hf:i:";
    const option long_opts[] = {
        { "fof", required_argument, nullptr, 'f' },
        { "info", required_argument, nullptr, 'i' },
        { "help", no_argument, nullptr, 'h' },
        { nullptr, no_argument, nullptr, 0 }
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt) {
            break;
        }
        switch (opt) {
        case 'f':
            paths.push_back(make_pair("fof",optarg));
            break;
        case 'i':
            paths.push_back(make_pair("path_info",optarg));
            break;
        case 'h':
        case '?': // Unrecognized option
        default:
            help();
            break;
        }
    }
    return paths;
}

uint16_t parseCoverage(const string& str) {
	size_t pos(str.find("km:f:"));
	if (pos == string::npos) {
		pos = (str.find("KM:f:"));
	}
	if (pos == string::npos) {
		return 1;
	}
	uint i(1);
	while (str[i + pos + 5] != ' ') {
		++i;
	}
	// cout<<"exact"<< (uint16_t)stof(str.substr(pos+5,i))<<endl;
	return (uint16_t)stof(str.substr(pos + 5, i));
}

string parse_filename(string& filename) {
    string clean_filename = "";
    auto last_slash_pos = filename.find_last_of('/');
    if (last_slash_pos != string::npos) {
        auto first_period_pos = filename.find_first_of('.', last_slash_pos);
        if (first_period_pos != string::npos) {
            clean_filename = filename.substr(last_slash_pos + 1, first_period_pos - last_slash_pos - 1);
        } else {
            clean_filename = filename;
        }
    } else {
        auto first_period_pos = filename.find_first_of('.');
        if (first_period_pos != string::npos) {
            clean_filename = filename.substr(0, first_period_pos);
        } else {
            clean_filename = filename;
        }
    }
    return clean_filename;
}

Info read_binary(string& in_file) {
    // ------ reading binary info file ------ //
    Info data_read_binary;
    ifstream info_file_bin(in_file,ios::binary);
    info_file_bin.seekg(info_file_bin.beg);
    info_file_bin.read(reinterpret_cast<char*>(&data_read_binary.nb_monotig), sizeof(uint64_t));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_binary.k), sizeof(uint));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_binary.record_option), sizeof(uint));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_binary.nb_eq_class), sizeof(long));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_binary.nb_colors), sizeof(uint64_t));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_binary.do_query_on_disk), sizeof(bool));
    info_file_bin.close();
    return data_read_binary;
}

vector<pair<string,uint64_t>> calc_kbf(string& fof) {
    // ------ calculating kmers_by_file ------ //

    ifstream ifof(fof);
    vector<string> input_files;
    while (!ifof.eof()) {
        string line;
        getline(ifof,line);
        if (fs::exists(line) && !fs::is_empty(line)) {
            input_files.push_back(line);
        }
    }
    vector<pair<string,uint64_t>> kbf(input_files.size(),make_pair("",0));
    for (uint32_t file = 0; file < input_files.size(); file++) {
        kbf[file].first = parse_filename(input_files[file]);
        auto reading = new ifstream(input_files[file]);
        string ref = "", header = "";
        while (!reading->eof()) {
            uint64_t occurrence_kmer = 0;
            getline(*reading, header);
            getline(*reading, ref);
            if (ref.size() < 31){
                ref = "";
            } else {
                occurrence_kmer = ref.size() - 30;
            }
            if (!header.empty() && !ref.empty()) {
                occurrence_kmer *= parseCoverage(header);
            }
            kbf[file].second += occurrence_kmer;
        }
        delete reading;
    }
    ifof.close();
    return kbf;
}

void write_info_txt(string& in_file, Info& data) {
    string extension = "_txt";
    ofstream info_file_txt(in_file+extension);
    info_file_txt << "nb_monotig:" << to_string(data.nb_monotig) << endl;
    info_file_txt << "k:" << to_string(data.k) << endl;
    info_file_txt << "record_option:" << to_string(data.record_option) << endl;
    info_file_txt << "nb_eq_class:" << to_string(data.nb_eq_class) << endl;
    info_file_txt << "nb_colors:" << to_string(data.nb_colors) << endl;
    info_file_txt << "do_query_on_disk:" << to_string(data.do_query_on_disk) << endl;
    for (auto& v : data.kbf) {
        info_file_txt << v.first << ":" << to_string(v.second) << endl;
    }
    info_file_txt.close();
}

void clear(string& in_file) {
    string extension = "_txt";
    fs::path path1 = in_file;
    fs::remove(in_file);
    fs::rename(in_file+extension,path1);
}

int main (int argc, char* argv[]) {
    if (argc == 1) {
        cout << strerror(EINVAL) << " : Missing argument(s) | Try -h or --help" << endl;
        exit(EINVAL);
    }
    vector<pair<string,string>> paths = process_args(argc, argv);
    string fof = "", in_file = "";
    for (auto& p : paths) {
        if (p.first == "fof") {
            fof = p.second;
        } else {
            auto last_underscore_pos = p.second.find_last_of('_');
            if (last_underscore_pos != string::npos) {
                if (p.second.substr(last_underscore_pos + 1) == "info") {
                    in_file = p.second;
                } else {
                    in_file = p.second+"/reindeer_matrix_eqc_info";
                }
            } else {
                in_file = p.second+"/reindeer_matrix_eqc_info";
            }
        }
    }
    Info data_read_binary = read_binary(in_file);
    data_read_binary.kbf = calc_kbf(fof); 
    write_info_txt(in_file, data_read_binary);
    clear(in_file);
    return 0;
}
