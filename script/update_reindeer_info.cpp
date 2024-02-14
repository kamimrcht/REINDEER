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
#include <system_error>
#include "../version.h" // manage version from git

using namespace std;
namespace fs = filesystem;

typedef struct info{
typedef struct filespath {
    string fof = "";        // store fof file = list of unitigs files
    string info_file = "";  // reindeer_matrix_eqc_info file to convert
    string info_new = "";   // reindeer_matrix_eqc_info file converted
    string info_path = "";  // path given for index directory
} Files_Path;
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
    "  -f, --fof <file>                 :       fof (File of files), file containing paths to files used to construct index\n"
    "  -i, --info <directory or file>   :       The directory where the index is stored\n"
    "\n"
    "* Optional parameters *\n"
    "  --text                           :       If the info file is already a text file and you want to update it again\n"
    "\n"
    "* General parameters *\n"
    "  --version, -V                    :       Show version\n"
    "  -h, --help                       :       Print help (what you are currently seeing)\n"
    << endl;
    exit(0);
}

Files_Path process_args(int argc, char** argv, bool& is_text) {

    Files_Path paths;

    const char* const short_opts = "hf:i:V";
    const option long_opts[] = {
        { "fof", required_argument, nullptr, 'f' },
        { "info", required_argument, nullptr, 'i' },
        { "version", no_argument, nullptr, 'V' },
        { "help", no_argument, nullptr, 'h' },
        { "text", no_argument, nullptr, 't' },
        { nullptr, no_argument, nullptr, 0 }
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt) {
            break;
        }
        switch (opt) {
        case 'f':
            paths.fof = optarg;
            break;
        case 'i':
            paths.info_path = optarg;
            break;
        case 't':
            is_text = true;
            break;
        case 'V':
            cout << VERSION << endl;
            exit(0);
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

bool is_binary(const string& filename) {
    ifstream file(filename, ios::binary);
    bool isBinary = false;
    char ch;
    while (file.get(ch)) {
        if (static_cast<unsigned char>(ch) > 127) {
        // If any byte value is greater than 127, it is likely binary data
            isBinary = true;
            break;
        }
    }
    file.close();
    return isBinary;
}

string parse_filename(const string& filename) {
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

Info read_binary(const string& in_file) {
    // ------ reading binary info file ------ //
    if (!is_binary(in_file)) {
        cerr << "[ERROR] " << strerror(EINVAL) << " : Input info file is not binary, try --text if you still want to update it or see help (-h, --help)" << endl;
        exit(EINVAL);
    }
    Info data_read_bin;
    ifstream info_file_bin(in_file,ios::binary);
    info_file_bin.seekg(info_file_bin.beg);
    info_file_bin.read(reinterpret_cast<char*>(&data_read_bin.nb_monotig), sizeof(uint64_t));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_bin.k), sizeof(uint));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_bin.record_option), sizeof(uint));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_bin.nb_eq_class), sizeof(long));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_bin.nb_colors), sizeof(uint64_t));
    info_file_bin.read(reinterpret_cast<char*>(&data_read_bin.do_query_on_disk), sizeof(bool));
    info_file_bin.close();
    cout << "Binary file read ✓" << endl;
    return data_read_bin;
}

Info read_text(const string& in_file){
    if (is_binary(in_file)) {
        cerr << "[ERROR] " << strerror(EINVAL) << " : Input info file is binary, try without --text if you want to update it or see help (-h, --help)" << endl;
        exit(EINVAL);
    }
    Info data_read_txt;
    ifstream info_file_txt(in_file);

    string line;
    getline(info_file_txt, line);
    size_t colon_pos = line.find(':');
    if (line.substr(0,colon_pos) != "nb_monotig") {
        cerr << "[ERROR]" << strerror(EINTR) << " | Update stopped : Missing variable 'nb_monotig', try with original binary info file" << endl;
        exit(EINTR);
    } else {
        data_read_txt.nb_monotig = stoull(line.substr(colon_pos + 1));
    }

    getline(info_file_txt, line);
    colon_pos = line.find(':');
    if (line.substr(0,colon_pos) != "k") {
        cerr << "[ERROR]" << strerror(EINTR) << " | Update stopped : Missing variable 'k', try with original binary info file" << endl;
        exit(EINTR);
    } else {
        data_read_txt.k = stoul(line.substr(colon_pos + 1));
    }

    getline(info_file_txt, line);
    colon_pos = line.find(':');
    if (line.substr(0,colon_pos) != "record_option") {
        cerr << "[ERROR]" << strerror(EINTR) << " | Update stopped : Missing variable 'record_option', try with original binary info file" << endl;
        exit(EINTR);
    } else {
        data_read_txt.record_option = stoul(line.substr(colon_pos + 1));
    }

    getline(info_file_txt, line);
    colon_pos = line.find(':');
    if (line.substr(0,colon_pos) != "nb_eq_class") {
        cerr << "[ERROR]" << strerror(EINTR) << " | Update stopped : Missing variable 'nb_eq_class', try with original binary info file" << endl;
        exit(EINTR);
    } else {
        data_read_txt.nb_eq_class = stol(line.substr(colon_pos + 1));
    }

    getline(info_file_txt, line);
    colon_pos = line.find(':');
    if (line.substr(0,colon_pos) != "nb_colors") {
        cerr << "[ERROR]" << strerror(EINTR) << " | Update stopped : Missing variable 'nb_colors', try with original binary info file" << endl;
        exit(EINTR);
    } else {
        data_read_txt.nb_colors = stoull(line.substr(colon_pos + 1));
    }

    getline(info_file_txt, line);
    colon_pos = line.find(':');
    if (line.substr(0,colon_pos) != "do_query_on_disk") {
        cerr << "[ERROR]" << strerror(EINTR) << " | Update stopped : Missing variable 'do_query_on_disk', try with original binary info file" << endl;
        exit(EINTR);
    } else {
        string do_query_on_disk_string = line.substr(colon_pos + 1);
        data_read_txt.do_query_on_disk = (do_query_on_disk_string != "0");
    }
    cout << "Text file read ✓" << endl;
    return data_read_txt;
}

vector<pair<string,uint64_t>> calc_kbf(const string& fof) {
    // ------ calculating kmers_by_file ------ //
    cout << "Calculating kmers by file..." << endl;
    ifstream ifof(fof);
    vector<string> input_files;
    while (!ifof.eof()) {
        string line;
        getline(ifof,line);
        if (fs::exists(line) && !fs::is_empty(line)) {
            input_files.push_back(line);
        } else
          cout << "Error with file: " << line << " not present or empty" << endl;
    }
    vector<pair<string,uint64_t>> kbf(input_files.size(),make_pair("",0));
    for (uint32_t file = 0; file < input_files.size(); file++) {
        //TODO: kbf[file].first = fs::path(input_files[file]).stem();
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
    cout << "Done ✓" << endl;
    return kbf;
}

void write_info_txt(const string& in_file, Info& data) {
    // rename old binary info file
    if (fs::exists(in_file) & is_binary(in_file))
      fs::rename(in_file, in_file + ".reindeer1.1");
    ofstream info_file_txt(in_file);
    info_file_txt << "nb_monotig:" << to_string(data.nb_monotig) << endl;
    info_file_txt << "k:" << to_string(data.k) << endl;
    info_file_txt << "record_option:" << to_string(data.record_option) << endl;
    info_file_txt << "nb_eq_class:" << to_string(data.nb_eq_class) << endl;
    info_file_txt << "nb_colors:" << to_string(data.nb_colors) << endl;
    info_file_txt << "do_query_on_disk:" << to_string(data.do_query_on_disk) << endl;
    for (auto& v : data.kbf) {
        info_file_txt << v.first << ":" << to_string(v.second) << endl;
    }
    cout << "Results written in ";
    info_file_txt.close();
}

/*
void clear(string& in_file) {
    string extension = "_txt";
    fs::path path1 = in_file;
    fs::remove(in_file);
    fs::rename(in_file+extension,path1);
}
*/

int main (int argc, char* argv[]) {
    if (argc == 1 || argc == 3 || argc == 4) {
        cout << strerror(EINVAL) << " : Missing argument(s) | Try -h or --help" << endl;
        exit(EINVAL);
    }
    bool is_text = false;
    Files_Path paths = process_args(argc, argv, is_text);
    // check if give directory for index files
    if (fs::is_directory(paths.info_path)) {
        paths.info_file = paths.info_path + "/reindeer_matrix_eqc_info";
        paths.info_new = paths.info_file + ".txt";
    } else if (fs::is_regular_file(paths.info_path)) {
        // given the info file ?
        string extension = paths.info_path.substr(paths.info_path.size() -4, 4);
        if (extension == "info") {
           paths.info_file = paths.info_path;
           paths.info_path = fs::path(paths.info_path).parent_path();
        } else if (extension == ".txt") {
           paths.info_new = paths.info_path;
           paths.info_file = paths.info_path.substr(0, paths.info_path.size()-4);
           paths.info_path = fs::path(paths.info_path).parent_path();
           is_text = true;
        }
    } else {
        cerr << "Error, the directory given for the reindeer index (" << paths.info_path << ") is not right" << endl;
        exit(1);
    }
    Info data_read;
    if (!is_text) {
        if (! fs::is_regular_file(paths.info_file)) {
            cerr << "Error: reindeer info file (" << paths.info_file << ") is not present, check it" << endl;
            exit(1);
        }
        data_read = read_binary(paths.info_file);
    } else {
        if (! fs::is_regular_file(paths.info_new)) {
            cerr << "Error: reindeer new info file (text format) (" << paths.info_new << ") is not present, check it" << endl;
            exit(1);
        }
        data_read = read_text(paths.info_new);
    }
    data_read.kbf = calc_kbf(paths.fof, data_read);
    write_info_txt(paths.info_new, data_read);
    cout << paths.info_path << "reindeer_matrix_eqc_info" << endl;
    //clear(in_file);
    return 0;
}
