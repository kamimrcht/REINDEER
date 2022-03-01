#include "utils.hpp"

using namespace std;

string get_run_tag()
{
    srand(time(NULL));
    int run_no(rand() % 10000);
    stringstream ss_rno;
    ss_rno << run_no;
    string rno = ss_rno.str();
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d-%m-%Y %H-%M-%S");
    string timestr = oss.str();
    vector<string> splitdate(split_utils(timestr, ' '));
    return splitdate[0] + "_" + rno;
}

void get_all_blout(string& path, vector<string>& files)
{
    DIR* dirp = opendir(path.c_str());

    struct dirent* dp;
    while ((dp = readdir(dirp)) != NULL) {
        files.push_back(dp->d_name);
    }
    closedir(dirp);
}

string do_fof(string& path, string& output)
{
    vector<string> bloutv;
    get_all_blout(path, bloutv);
    string fname(output + "/fof_blout");
    ofstream out(fname);
    for (auto&& f : bloutv)
        if (f != "." and f != "..")
            out << path + f << endl;
    out.close();
    return fname;
}

// convert char [] counts/colors to uint
vector<uint16_t> count_string_to_count_vector(unsigned char* count_char_monotig, unsigned size)
{
    vector<uint16_t> counts;
    for (uint i(0); i < size - 1; i += 2) {
        counts.push_back((uint16_t)count_char_monotig[i] + (uint16_t)count_char_monotig[i + 1] * 256);
    }
    return counts;
}

// convert char [] counts/colors to uint
vector<uint8_t> count_string_to_count_vector8(unsigned char* count_char_monotig, unsigned size)
{
    vector<uint8_t> counts;
    for (uint i(0); i < size; i++) {
        counts.push_back((uint8_t)count_char_monotig[i]);
    }
    return counts;
}

uint64_t getMemorySelfUsed()
{
    u_int64_t result = 0;
    FILE* file = fopen("/proc/self/status", "r");
    if (file) {
        char line[128];

        while (fgets(line, 128, file) != NULL) {
            if (strncmp(line, "VmRSS:", 6) == 0) {
                char* loop = line;
                result = strlen(line);
                while (*loop < '0' || *loop > '9')
                    loop++;
                loop[result - 3] = '\0';
                result = atoi(loop);
                break;
            }
        }
        fclose(file);
    } else
        cout << "warning: could not fopen /proc/self/status" << std::endl;
    return result;
}

uint64_t getMemorySelfMaxUsed()
{
    u_int64_t result = 0;
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        result = usage.ru_maxrss;
    }
    return result;
}

uint64_t xorshift(uint64_t x)
{
    x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
    x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
    x = ((x >> 32) ^ x);
    return x;
}

void throw_character_issue()
{
    cerr << "[Warning] Sequence containing a 'N' character or invalid header is disregarded." << endl;
}

void throw_size_issue()
{
    cerr << "[Warning] Sequence smaller than k is disregarded." << endl;
}

void throw_empty_issue()
{
    cerr << "[Warning] Query file contains empty lines. Please amend it before querying REINDEER" << endl;
}

//~ void getLineFasta(ifstream* in, string& fasta, string& header) {
//~ void getLineFasta_buffer(ifstream* in, string& fasta, string& header, uint buffer) {
string getLineFasta_buffer(ifstream* in)
{
    string line, result;
    getline(*in, line);
    char c = static_cast<char>(in->peek());
    if (line[0] == '>') {
        return line;
    } else {
        if (not check_character(line))
            result += line;
        while (c != '>' and c != EOF) {
            getline(*in, line);
            if (not check_character(line))
                result += line;
            c = static_cast<char>(in->peek());
        }
        return result;
    }
}

vector<string> getLineFasta_buffer2(ifstream* in, uint stop, uint k)
{
    vector<string> lines;
    string line, result;
    uint i(0);
    bool error(false);
    char c;
    c = static_cast<char>(in->peek());
    while (not(i >= stop or in->eof())) {
        ++i;
        if (c == '>') // header
        {
            if (not result.empty()) {
                if (result.size() >= k) {
                    lines.push_back(result);
                } else {
                    lines.pop_back(); //remove last header corresponding to too small sequence
                    throw_size_issue();
                }
                result = "";
            }
            getline(*in, line);
            if (line.empty()) {
                throw_empty_issue();
                break;
            };
            c = static_cast<char>(in->peek());
            lines.push_back(line);
            error = false;
        } else // sequence
        {
            while (c != '>' and c != EOF) {
                getline(*in, line);
                if (line.empty()) {
                    throw_empty_issue();
                    break;
                };
                {

                    if (check_character(line)) //invalid character, do no include read
                    {
                        lines.pop_back(); //remove last header corresponding to an erroneous sequence
                        result = "";
                        throw_character_issue();
                        c = static_cast<char>(in->peek());
                        error = true;
                        break;
                    } else {
                        if (lines.back()[0] == '>' and (not error))
                            result += line;
                    }
                }

                c = static_cast<char>(in->peek());
            }
        }
    }
    if (not lines.empty()) {
        if (lines.back()[0] == '>') {
            if (not result.empty()) {
                if (result.size() >= k) {
                    lines.push_back(result);
                    result = "";
                } else {
                    lines.pop_back(); //remove last header corresponding to too small sequence
                    throw_size_issue();
                }
            } else {
                lines.pop_back();
            }
        }
        if (lines.back().size() < k) {
            lines.pop_back();
            lines.pop_back();
            throw_size_issue();
        }
    }
    return lines;
}

bool is_empty_file(ifstream& file)
{
    return file.peek() == ifstream::traits_type::eof();
}

bool is_empty_zfile(zstr::ifstream& file)
{
    return file.peek() == zstr::ifstream::traits_type::eof();
}

int dirExists(string& path)
{
    struct stat info;

    int statRC = stat(path.c_str(), &info);
    if (statRC != 0) {
        if (errno == ENOENT) {
            return 0;
        } // something along the path does not exist
        if (errno == ENOTDIR) {
            return 0;
        } // something in path prefix is not a dir
        return -1;
    }

    return (info.st_mode & S_IFDIR) ? 1 : 0;
}

bool check_character(string& s)
{
    if (s[0] != '>') {
        for (auto&& c : s) {
            if (c != 'A' and c != 'T' and c != 'G' and c != 'C') {
                return true;
            }
        }
    }
    return false;
}

//~ inline bool exists_test(const string& name) {
//~ return ( access( name.c_str(), F_OK ) != -1 );
//~ }

vector<string> split_utils(const string& s, char delim)
{
    stringstream ss(s);
    string item;
    vector<string> elems;
    while (getline(ss, item, delim)) {
        elems.push_back(move(item));
    }
    return elems;
}

string get_file_name(string& path)
{
    string result(""), tmp("");
    vector<string> split(split_utils(path, '/'));
    if (not split.empty()) {
        tmp = split.back();
        split = split_utils(tmp, '.');
        if (not split.empty())
            result = split_utils(tmp, '.')[0];
    }
    return result;
}

///// parse bcalm headers //////////////
double parseCoverage_utils(const string& str)
{
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
    // WARNING WE RETURN COUNTS BETWEEN 0 AND 255
    float result(stof(str.substr(pos + 5, i)));
    if (result > 255) {
        return 255;
    } else {
        return result;
    }
}

uint32_t unitig_toui32(const string& s)
{
    uint32_t val(0);
    int32_t valI(stoi(s));
    val = (valI >= 0) ? (uint32_t)valI : (uint32_t)-valI;
    return val;
}

void parse_bgreat_output(string& input, vector<vector<uint64_t>>& unitigs_to_nodes)
{
    ifstream bgreat_file(input);
    string line;
    vector<string> unitigs_per_read;
    uint64_t readID(0);

    while (not bgreat_file.eof()) {
        getline(bgreat_file, line);
        if (line.empty()) {
            readID++;
            continue;
        }
        unitigs_per_read = split_utils(line, ';');
        for (auto&& u : unitigs_per_read) {
            uint32_t unitig(unitig_toui32(u) - 1);
            if (unitig >= unitigs_to_nodes.size()) {
                unitigs_to_nodes.resize(unitig + 1);
            }
            unitigs_to_nodes[unitig].push_back(readID);
        }
        readID++;
    }
}

vector<unsigned char> RLE16C(const vector<uint16_t>& V)
{
    vector<unsigned char> res;
    if (V.empty()) {
        return res;
    }
    uint16_t pred(V[0]);

    uint8_t count(1);
    for (uint64_t i(1); i < V.size(); ++i) {
        if (V[i] == pred) {
            count++;
            if (count == 254) {
                res.push_back(min(pred / 256, 254));
                res.push_back(min(pred % 256, 254));
                res.push_back(count);
                count = 1;
            }
        } else {
            res.push_back(min(pred / 256, 254));
            res.push_back(min(pred % 256, 254));
            res.push_back(count);
            count = 1;
            pred = V[i];
        }
    }
    res.push_back(min(pred / 256, 254));
    res.push_back(min(pred % 256, 254));
    res.push_back(count);
    //~ vector<unsigned char> res_ca;
    //~ transform(begin(res), end(res), begin(res_ca), [](uint8_t i) { return '0' + i; });

    return res;
}

vector<uint16_t> RLE16D(const vector<uint8_t>& V)
{
    vector<uint16_t> res;
    if (V.size() < 2) {
        return res;
    }
    for (uint64_t i(0); i < V.size(); i += 3) {
        res.resize(V[i + 2], V[i] * 256 + V[i + 1]);
    }
    return res;
}

vector<uint8_t> RLE8C(const vector<uint8_t>& V)
{
    vector<uint8_t> res;
    if (V.empty()) {
        return res;
    }
    uint8_t pred(V[0]);
    uint8_t count(1);
    for (uint64_t i(1); i < V.size(); ++i) {
        if (V[i] == pred) {
            count++;
            if (count == 254) {
                res.push_back(pred);
                res.push_back(count);
                count = 1;
            }
        } else {
            res.push_back(pred);
            res.push_back(count);
            count = 1;
            pred = V[i];
        }
    }
    res.push_back(pred);
    res.push_back(count);
    return res;
}

vector<uint8_t> RLE8D(const vector<uint8_t>& V)
{
    vector<uint8_t> res;
    if (V.size() < 2) {
        return res;
    }
    for (uint64_t i(0); i < V.size(); i += 2) {
        res.resize(V[i + 1], V[i]);
    }
    return res;
}

void new_paired_end_file(string& input, string& input2, string& output_file, bool fastq)
{
    string line, head, head2, sequence, junk;
    istream* in;
    istream* in2;
    in = new zstr::ifstream(input);
    in2 = new zstr::ifstream(input2);
    ofstream out(output_file);
    vector<uint> lengths;
    while (not in->eof()) {
        getline(*in, head);
        getline(*in2, head2);
        getline(*in, sequence);
        if (fastq) {
            getline(*in, junk);
            getline(*in, junk);
            head[0] = '>';
        }
        out << head << "\n";
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        out << sequence << "\n";
        sequence = "";
        getline(*in2, sequence);
        if (fastq) {
            getline(*in2, junk);
            getline(*in2, junk);
            head2[0] = '>';
        }
        out << head2 << "\n";
        transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
        out << sequence << "\n";
        sequence = "";
    }
    delete (in);
    delete (in2);
}

void interleave_paired_end(string& fof, string& output)
{
    bool tested(false), fastq(true);
    ifstream file(fof);
    string new_fof_name(output + "/pe_fof.lst");
    string sample, sample2, header, output_name;
    uint file_index(0);
    ofstream new_fof(new_fof_name);
    while (not file.eof()) {
        getline(file, sample);
        if (sample.empty()) {
            break;
        }
        getline(file, sample2);
        if (sample.empty()) {
            break;
        }
        if (not tested) {
            ifstream samp(sample);
            while (not file.eof()) {
                getline(samp, header);
                if (header.empty()) {
                    break;
                }
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
    if (counts.size() > 0) {
        float harmonicMean(0);
        for (uint i(0); i < counts.size(); ++i) {
            if (counts[i] > 0) {
                harmonicMean += 1 / (float)(counts[i]);
            }
        }
        if (harmonicMean > 0) {
            return (uint64_t)(counts.size() / harmonicMean);
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}

string getRealPath(string file, string& dir)
{
    char* symlinkpath = &dir[0];
    char actualpath[PATH_MAX + 1];
    char* ptr;
    ptr = realpath(symlinkpath, actualpath);
    string rp(ptr);
    return rp + "/" + file;
}

uint64_t get_color_number(string& fof)
{
    uint color(0);
    string line;
    ifstream fof_file(fof);
    while (not fof_file.eof()) {
        getline(fof_file, line);
        if (not line.empty()) {
            color++;
        }
    }
    return color;
}
