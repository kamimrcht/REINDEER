#include "query.hpp"
using namespace std;
using namespace chrono;

// decode rl encoded vector (color/count of ONE monotig)
unsigned char* decode_vector(unsigned char* monotig_counts, unsigned vector_size, uint64_t color_number, bool record_counts)
{
    unsigned char* decoded_vector;
    unsigned sz;
    if (record_counts) {
        decoded_vector = new unsigned char[color_number * 2 + 4096];
        sz = trled(monotig_counts, vector_size, decoded_vector, color_number * 2);
    } else {
        decoded_vector = new unsigned char[color_number + 1024];
        sz = trled(monotig_counts, vector_size, decoded_vector, color_number);
    }
    return decoded_vector;
}

// for ONE monotig, get a vector of its counts/colors (uint) from the rle encoded matrix
vector<uint16_t> get_count_monotig(unsigned char* monotig_counts, unsigned vector_size, uint64_t color_number, bool record_counts)
{
    unsigned char* decoded(decode_vector(monotig_counts, vector_size, color_number, record_counts));
    vector<uint16_t> counts(count_string_to_count_vector(decoded, color_number * 2));
    delete[] decoded;
    return counts;
}

template <class T>
void Reindeer_Index<T>::get_position_vector_query_disk(vector<long>& position_in_file)
{
    ifstream in(matrix_eqc_position_file); //positions file
    long position;
    uint nb(0);
    while (nb < nb_monotig) {
        in.read(reinterpret_cast<char*>(&position), sizeof(long));
        position_in_file.push_back(position);
        ++nb;
    }
    in.close();
}

long get_matrix_line_query_disk(int64_t rank, unsigned char* color, unsigned& line_size, vector<long>& position_in_file, ifstream& in)
{
    int64_t rank2;
    long pos = position_in_file[rank];
    in.seekg(pos);
    in.read(reinterpret_cast<char*>(&line_size), sizeof(unsigned));
    in.read(reinterpret_cast<char*>(&rank2), sizeof(int64_t));
    in.read((char*)(color), line_size);
    return pos;
}

long get_matrix_line_query(int64_t rank, unsigned char* color, unsigned& line_size, vector<long>& position_in_file, vector<unsigned char*>& compr_monotig_color)
{
    long pos = position_in_file[rank];
    color = compr_monotig_color[pos];
    return pos;
}

void get_colors_counts_query_eq_classes(vector<int64_t>& kmer_ids, uint64_t color_number, vector<vector<uint16_t>>& query_counts, vector<unsigned char*>& compr_monotig_color, vector<unsigned>& compr_monotig_color_size, vector<long>& position_in_file, bool record_counts, vector<vector<uint8_t>>& query_colors, bool do_query_on_disk, string& rd_file)
{
    ifstream in(rd_file);
    unsigned char* color;
    color = new unsigned char[2 * color_number + 1204];
    unsigned size;
    int64_t lastId(-1);
    vector<uint16_t> qcounts, lastV;
    vector<uint8_t> qcolors;
    vector<uint16_t> vec(color_number, 0);

    for (uint64_t i(0); i < kmer_ids.size(); ++i) {
        if (kmer_ids[i] >= 0) {
            // no need to compute for a kmer if it comes from the same monotig than the previous, just copy the result
            if (kmer_ids[i] == lastId) {
                lastId = kmer_ids[i];
                query_counts.push_back(lastV);
            } else {
                long pos;
                unsigned char* lo;
                if (do_query_on_disk) {
                    pos = get_matrix_line_query_disk(kmer_ids[i], color, size, position_in_file, in);
                    lo = decode_vector((unsigned char*)&color[0], size, color_number, record_counts);
                } else {
                    pos = get_matrix_line_query(kmer_ids[i], color, size, position_in_file, compr_monotig_color);
                    lo = decode_vector((unsigned char*)&compr_monotig_color[pos][0], compr_monotig_color_size[pos], color_number, record_counts);
                }
                if (record_counts) {
                    qcounts = count_string_to_count_vector(lo, color_number * 2);
                } else {
                    qcounts = {};
                    vector<uint8_t> qcounts8 = count_string_to_count_vector8(lo, color_number);
                    copy(qcounts8.begin(), qcounts8.end(), back_inserter(qcounts));
                }
                delete[] lo;
                lastV = qcounts;
                if (not qcounts.empty()) {
                    query_counts.push_back(qcounts);
                } else {
                    query_counts.push_back(vec);
                }
            }
        } else {
            query_counts.push_back(vec);
        }
    }
    delete[] color;
    in.close();
}

// for all queried k-mers, get the colors/counts in vector<vector<uint16_t>>& query_counts
void get_colors_counts(vector<int64_t>& kmer_ids, bool record_counts, uint64_t color_number, vector<int64_t>& kmers_colors, vector<vector<uint16_t>>& query_counts, vector<unsigned char*>& compr_monotig_color, vector<unsigned>& compr_monotig_color_size, vector<long>& position_in_file)
{
    vector<uint16_t> counts;
    int64_t lastId(-1);
    vector<uint16_t> qcounts, lastV;
    for (int64_t i(0); i < kmer_ids.size(); ++i) {
        if (kmer_ids[i] >= 0) {
            // no need to compute for a kmer if it comes from the same monotig than the previous, just copy the result
            if (kmer_ids[i] == lastId) {
                lastId = kmer_ids[i];
                query_counts.push_back(lastV);
            } else {
                qcounts = get_count_monotig(compr_monotig_color[kmer_ids[i]], compr_monotig_color_size[kmer_ids[i]], color_number, record_counts);
                lastV = qcounts;
                if (not qcounts.empty())
                    query_counts.push_back(qcounts);
            }
        }
    }
}

// compute a string that sums up the count(s) for each dataset
vector<uint> write_count_output(bool average, bool record_counts, vector<vector<uint16_t>>& query_counts, uint64_t color_number, vector<string>& toW, vector<string>& color_counts, uint k_size)
{
    vector<uint> covered_positions;
    string nc("*");
    vector<string> last(color_number, "*");
    for (auto&& c : query_counts) {
        for (uint color(0); color < c.size(); ++color) {
            nc = to_string(c[color]);
            if (nc == "0")
                nc = "*";
            if (last[color] != nc) {
                toW[color] += nc + ":";
                last[color] = nc;
            }
        }
    }
    uint covered_pos = 0;

    for (uint color(0); color < color_number; ++color) {
        string out_str("");
        uint total_for_average = 0, out_uint = 0;
        vector<pair<pair<uint, uint>, string>> coords;

        for (uint c(0); c < query_counts.size(); ++c) {
            string val_str = to_string(query_counts[c][color]);
            if (val_str == "0")
                val_str = "*";
            if (c == 0) {
                coords.push_back({ { 0, 0 }, val_str });
            } else {

                if (val_str == coords.back().second) {
                    coords.back().first.second++;
                } else {
                    uint new_coord = coords.back().first.second + 1;
                    coords.push_back({ { new_coord, new_coord }, val_str });
                }
            }
            nc = to_string(query_counts[c][color]);
        }
        uint cov_positions(query_counts.size() + k_size - 1);
        for (auto&& coo : coords) {
            if (average) {
                if (coo.second != "*") {
                    out_uint += (coo.first.second-coo.first.first+1)*stoi(coo.second);
                }
                total_for_average += (coo.first.second-coo.first.first+1);
            } else {
                out_str += to_string(coo.first.first) + "-" + to_string(coo.first.second) + ":" + coo.second + ",";
            }
            if (coo.second == "*") {
                cov_positions -= (coo.first.second - coo.first.first + 1);
            }
        }
        if (average) {
            out_uint = out_uint/total_for_average;
            out_str = to_string(out_uint);
        } else {
            out_str.pop_back(); //remove last comma
        }
        color_counts.push_back(out_str);
        covered_positions.push_back(cov_positions * 100 / (query_counts.size() + k_size - 1));
    }
    return covered_positions;
}

void write_results_above_threshold(string& toWrite, vector<vector<uint16_t>>& query_counts, uint64_t color_number, vector<string>& toW, vector<string>& color_counts, string& header, bool record_counts, uint threshold, string& line, uint k, vector<uint>& covered_positions)
{
    vector<double_t> percent(color_number, 0);

    toWrite += header.substr(0, 50);

    for (uint cp(0); cp < covered_positions.size(); ++cp) {
        if (covered_positions[cp] >= threshold) {
            if (record_counts) {
                toWrite += "\t" + color_counts[cp];
            } else {
                toWrite += "\t" + to_string(covered_positions[cp]);
            }
        } else {
            toWrite += "\t*";
        }
    }
    toWrite += "\n";
}

void write_output(vector<int64_t>& kmers_colors, string& toWrite, bool average, bool record_counts, vector<vector<uint32_t>>& query_unitigID, vector<vector<uint32_t>>& query_unitigID_tmp, uint64_t& color_number, string& header, string& line, uint k, uint threshold, vector<vector<uint16_t>>& query_counts, vector<vector<uint8_t>>& query_colors)
{
    vector<string> color_counts;
    vector<string> toW(color_number, "");
    vector<uint> covered_positions = write_count_output(average, record_counts, query_counts, color_number, toW, color_counts, k);
    write_results_above_threshold(toWrite, query_counts, color_number, toW, color_counts, header, record_counts, threshold, line, k, covered_positions);
}

void doQuery(string& input, string& name, kmer_Set_Light& ksl, uint64_t& color_number, uint k, bool average, bool record_counts, uint threshold, vector<vector<uint32_t>>& query_unitigID, uint nb_threads, vector<unsigned char*>& compr_monotig_color, vector<unsigned>& compr_monotig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb, uint64_t nb_monotig, vector<long>& position_in_file, string& fof)
{

    ifstream query_file(input);
    ofstream out(name);
    string initW;
    init_outputfile(initW, fof); //todo give fof file (todo retain fof file name in attributes) and output
    out << initW << endl;
    uint64_t num_seq(0);
    string qline;
    mutex mm;
    vector<string> lines;
    vector<vector<uint32_t>> query_unitigID_tmp;
    // FOR EACH LINE OF THE QUERY FILE
    bool first(true);
    while (not query_file.eof()) {
#pragma omp parallel num_threads(nb_threads)
        {
#pragma omp critical(i_file)
            {
                //~ uint i(0);
                lines = getLineFasta_buffer2(&query_file, 4000, k);
            }
            uint i;
#pragma omp for ordered
            for (i = (0); i < lines.size(); i += 2) {
                uint j(i);

                if (i % 1000 == 0)
                    cout << "-";
                string toWrite;
                string header;
                while (j < i + 2) {
                    string line = lines[j];
                    if (j % 2 == 1) {
                        vector<int64_t> kmers_colors;
                        vector<string> color_counts;
                        vector<int64_t> kmer_ids;
                        kmer_ids = ksl.get_rank_query(line);
                        vector<vector<uint16_t>> query_counts;
                        vector<vector<uint8_t>> query_colors;
                        get_colors_counts_query_eq_classes(kmer_ids, color_number, query_counts, compr_monotig_color, compr_monotig_color_size, position_in_file, record_counts, query_colors, do_query_on_disk, rd_file);
                        mm.lock();
                        write_output(kmers_colors, toWrite, average, record_counts, query_unitigID, query_unitigID_tmp, color_number, header, line, k, threshold, query_counts, query_colors);
                        mm.unlock();
                    } else {
                        header = line;
                    }
                    j++;
                }
#pragma omp ordered
                mm.lock();
                if (toWrite != header + "\n") {
                    out << toWrite;
                }
                mm.unlock();
            }
        }
        lines = {};
    }
    cout << endl;
    query_file.close();
    out.close();
}

void query_by_file(uint& counter, string& entry, kmer_Set_Light& ksl, uint64_t& color_number, uint k, bool average, bool record_counts, uint threshold, string& output, uint nb_threads, vector<unsigned char*>& compr_monotig_color, vector<unsigned>& compr_monotig_color_size, bool do_query_on_disk, string& rd_file, long eq_class_nb, uint64_t nb_monotig, vector<long>& position_in_file, string& fof)
{
    string outName(output + "/out_query_Reindeer_P" + to_string(threshold) + "_" + get_file_name(entry).substr(0, 50) + "_" + to_string(counter) + ".out");
    cout << "Result will be written in " << outName << endl;
    vector<vector<uint32_t>> query_unitigID(color_number, { 0 });
    doQuery(entry, outName, ksl, color_number, k, average, record_counts, threshold, query_unitigID, nb_threads, compr_monotig_color, compr_monotig_color_size, do_query_on_disk, rd_file, eq_class_nb, nb_monotig, position_in_file, fof);
    counter++;
}

template <class T>
void Reindeer_Index<T>::perform_query(kmer_Set_Light& ksl, uint threshold, string& query, vector<long>& position_in_file, bool average)
{
    uint counter(0), patience(0);
    string entry;
    if (not query.empty()) {
        if (exists_test(query)) {
            high_resolution_clock::time_point t121 = high_resolution_clock::now();
            query_by_file(counter, query, ksl, nb_colors, k, average, record_counts, threshold, output, threads, compressed_monotig_color, compressed_monotig_color_sizes, do_query_on_disk, matrix_name, nb_eq_class, nb_monotig, position_in_file, fof_file);
            high_resolution_clock::time_point t13 = high_resolution_clock::now();
            duration<double> time_span13 = duration_cast<duration<double>>(t13 - t121);
            cout << "Query done: " << time_span13.count() << " seconds." << endl;
        } else {
            cout << "Empty query file !" << endl;
        }
    } else {

        while (true) {
            char str[256];
            cout << "Enter the name of the query file (and optionally -P option separated by a space): ";

            cin.get(str, 256);
            cin.clear();
            cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            //~ entry = str;
            vector<string> entries = split_utils(str, ' ');
            entry = entries[0];
            if (entry == "") {
                if (patience > 0) {
                    cout << "See you soon !" << endl;
                    exit(0);
                } else {
                    cout << "Empty query file !  Type return again if you  want to quit" << endl;
                    patience++;
                }
            } else {
                uint thresh(threshold);
                patience = 0;
                if (exists_test(entry)) {
                    if (entries.size() > 1) {
                        if (entries[1] != "") {
                            thresh = stoi(entries[1]);
                        }
                    }
                    cout << "Queried file: " << entry << " option -P: " << thresh << endl;
                    high_resolution_clock::time_point t121 = high_resolution_clock::now();
                    query_by_file(counter, entry, ksl, nb_colors, k, average, record_counts, thresh, output, threads, compressed_monotig_color, compressed_monotig_color_sizes, do_query_on_disk, matrix_name, nb_eq_class, nb_monotig, position_in_file, fof_file);
                    memset(str, 0, 255);
                    high_resolution_clock::time_point t13 = high_resolution_clock::now();
                    duration<double> time_span13 = duration_cast<duration<double>>(t13 - t121);
                    cout << "Query done: " << time_span13.count() << " seconds." << endl;
                } else {
                    cout << "The entry is not a file" << endl;
                }
            }
        }
    }
}
