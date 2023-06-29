#include "monotig.hpp"

using namespace std;
using namespace chrono;

uint64_t MAX_ABUNDANCE_DISCRETE;
vector<uint8_t> abundance_discretization;

// formula to retrieve the estimated abundance from a 8 bits quantized count
uint16_t kmer_Set_Light::abundance_at(uint8_t index)
{
    if (index < abundance_discretization.size()) {
        return floorf((abundance_discretization[index] + abundance_discretization[index + 1]) / 2.0);
    } else {
        return 0;
    }
}

//initialisation of the table : count to discrete count
void kmer_Set_Light::init_discretization_scheme()
{
    MAX_ABUNDANCE_DISCRETE = 65536;
    abundance_discretization.resize(257);
    int total = 0;
    abundance_discretization[0] = 0;
    int idx = 1;
    for (int ii = 1; ii <= 70; ii++, idx++) {
        total += 1;
        abundance_discretization[idx] = total;
    }

    for (int ii = 1; ii <= 15; ii++, idx++) {
        total += 2;
        abundance_discretization[idx] = total;
    }

    for (int ii = 1; ii <= 40; ii++, idx++) {
        total += 10;
        abundance_discretization[idx] = total;
    }
    for (int ii = 1; ii <= 25; ii++, idx++) {
        total += 20;
        abundance_discretization[idx] = total;
    }
    for (int ii = 1; ii <= 40; ii++, idx++) {
        total += 100;
        abundance_discretization[idx] = total;
    }
    for (int ii = 1; ii <= 25; ii++, idx++) {
        total += 200;
        abundance_discretization[idx] = total;
    }
    for (int ii = 1; ii <= 40; ii++, idx++) {
        total += 1000;
        abundance_discretization[idx] = total;
    }
    abundance_discretization[256] = total;
}

// discretization of one count
uint8_t kmer_Set_Light::return_count_bin(uint16_t abundance)
{
    int idx;
    if (abundance >= MAX_ABUNDANCE_DISCRETE) {
        //~ _nb_abundances_above_precision++;
        //std::cout << "found abundance larger than discrete: " << abundance << std::endl;
        idx = abundance_discretization.size() - 2;
    } else {
        //get first cell strictly greater than abundance
        auto up = upper_bound(abundance_discretization.begin(), abundance_discretization.end(), abundance);
        up--; // get previous cell
        idx = up - abundance_discretization.begin();
    }
    return idx;
}

// parse bcalm
uint16_t kmer_Set_Light::parseCoverage_bin(const string& str)
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
    return return_count_bin((uint16_t)stof(str.substr(pos + 5, i)));
}

// tests if 2 vectors have the same null positions
bool equal_nonull(const vector<uint16_t>& V1, const vector<uint16_t>& V2)
{
    if (V1.empty()) {
        return false;
    }
    if (V2.size() != V1.size()) {
        return false;
    }
    for (uint i(0); i < V1.size(); ++i) {
        if ((V1[i] == 0) ^ (V2[i] == 0)) {
            return false;
        }
    }
    return true;
}

// deprecated
bool kmer_Set_Light::similar_count(const vector<uint16_t>& V1, const vector<uint16_t>& V2)
{
    if (V1.empty()) {
        return false;
    }
    if (V2.size() != V1.size()) {
        return false;
    }
    for (uint i(0); i < V1.size(); ++i) {
        if (((double)max(V1[i], V2[i]) / ((double)min(V1[i], V2[i])) > max_divergence_count)) {
            return false;
        }
    }
    return true;
}

// RLE on vector
string compress_vector(const vector<uint16_t>& V)
{
    string res;
    unsigned char comp[V.size() * 2 + 1024];
    uint32_t rle_length = trlec((unsigned char*)&V[0], V.size() * 2, comp);
    res.assign((const char*)comp, rle_length);
    return res;
}

uint64_t max_size_bucket(0);
uint64_t min_size_bucket(10000000000);
uint64_t max_size_superbucket(0);
uint64_t min_size_superbucket(10000000000);

// builds blight index from fof
void kmer_Set_Light::construct_index_fof(const string& input_file, vector<uint64_t>& kmers_by_file, const string& tmp_dir, int colormode)
{
    omp_set_nested(2);
    if (not tmp_dir.empty()) {
        working_dir = tmp_dir + "/";
    }
    color_mode = colormode;
    if (color_mode == 2) {
        init_discretization_scheme();
    }
    if (m1 < m2) {
        cout << "n should be inferior to m" << endl;
        exit(0);
    }
    if (m2 < m3) {
        cout << "s should be inferior to n" << endl;
        exit(0);
    }
    nuc_minimizer = new uint32_t[minimizer_number.value()];
    current_pos = new uint64_t[minimizer_number.value()];
    start_bucket = new uint64_t[minimizer_number.value()];
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    ifstream infof(input_file);
    vector<string> fnames;
    //STACK SUPERBUCKETS
    while (not infof.eof()) {
        string file;
        getline(infof, file);
        if (not file.empty() and exists_test(file)) {
            fnames.push_back(file);
        }
    }
    // files for the _blout that contain monotigs
    vector<ofstream*> out_blout;
    create_super_buckets_list(fnames, kmers_by_file);
    reset();

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    cout << "Partition created	" << intToString(read_kmer) << " kmers read " << endl;
    duration<double> time_span12 = duration_cast<duration<double>>(t2 - t1);
    cout << time_span12.count() << " seconds." << endl;

    string monotig_folder(working_dir + "monotig_files/");
    int systRet;
    if (not dirExists(monotig_folder)) {
        systRet = system(("mkdir " + monotig_folder).c_str());
    }
    {
#pragma omp parallel for num_threads(coreNumber)
        for (uint i_superbuckets = 0; i_superbuckets < number_superbuckets.value(); ++i_superbuckets) {
            string in_name(working_dir + "_super_k_mers" + to_string(i_superbuckets) + ".gz");
            string out_name(monotig_folder + "_blout" + to_string(i_superbuckets) + ".gz");
            merge_super_buckets_mem(in_name, fnames.size(), out_name, 1, colormode);
            remove((in_name).c_str());
            cout << "-" << flush;
        }
    }
    initialize_buckets();
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    cout << "\nMonotigs computed, now regular indexing start" << endl;
    duration<double> time_span32 = duration_cast<duration<double>>(t3 - t2);
    cout << time_span32.count() << " seconds." << endl;
    read_super_buckets_reindeer(monotig_folder + "_blout");
    delete[] nuc_minimizer;
    delete[] start_bucket;
    delete[] current_pos;

    high_resolution_clock::time_point t5 = high_resolution_clock::now();
    duration<double> time_span53 = duration_cast<duration<double>>(t5 - t3);
    cout << "Index created: " << time_span53.count() << " seconds." << endl;
    duration<double> time_spant = duration_cast<duration<double>>(t5 - t1);
    cout << "The whole indexing took me " << time_spant.count() << " seconds." << endl;
}

// from the vector of pairs that stores colors in the kmers in the super buckets (during monotig construction), builds the real color vector (a vector of size c with c datasets)
vector<uint8_t> get_color_vector(const vector<pair<uint16_t, uint16_t>>& V, uint64_t number_color)
{
    vector<uint8_t> res(number_color, 0);
    for (uint64_t i(0); i < V.size(); ++i) {
        res[V[i].first] = (uint8_t)V[i].second;
    }
    return res;
}
// from the vector of pairs that stores colors in the kmers in the super buckets (during monotig construction), builds the real color vector (a vector of size c with c datasets)
vector<uint16_t> get_count_vector(const vector<pair<uint16_t, uint16_t>>& V, uint64_t number_color)
{
    vector<uint16_t> res(number_color, 0);
    for (uint64_t i(0); i < V.size(); ++i) {
        res[V[i].first] = V[i].second;
    }
    return res;
}

// build monotigs within a super bucket
void kmer_Set_Light::merge_super_buckets_mem(const string& input_file, uint64_t number_color, string& out_name, uint64_t number_pass, int colormode)
{
    zstr::ofstream* out = new zstr::ofstream(out_name); // todo lz4
    bool toobig(false);
    for (uint pass(0); pass < number_pass; ++pass) {
        vector<robin_hood::unordered_node_map<kmer, kmer_context>> min2kmer2context(minimizer_number.value() / number_superbuckets.value());
        uint64_t mms = min2kmer2context.size();
        vector<int32_t> minimizers(mms, -1);
        zstr::ifstream in(input_file);
        uint64_t inserted_elements(0);

        // #pragma omp parallel num_threads(coreNumber)
        {
            string line, sequence;
            uint32_t color;
            uint16_t coverage;
            uint32_t min_int;
            uint8_t sizel(0);
            bool stop(false);
            while (not in.eof() and in.good() and not toobig) {
                // #pragma omp critical
                {
                    in.read(reinterpret_cast<char*>(&min_int), 4);
                    in.read(reinterpret_cast<char*>(&color), 4);
                    in.read(reinterpret_cast<char*>(&coverage), 2);
                    in.read(reinterpret_cast<char*>(&sizel), 1);
                    sequence.resize(sizel);
                    in.read(reinterpret_cast<char*>(&sequence[0]), sizel);
                    if (in.eof()) {
                        stop = true;
                    } else {
                        minimizers[min_int % mms] = (min_int);
                    }
                }
                if (not stop) {
                    uint64_t indice(min_int % mms);
                    if ((number_pass != 1) and (indice % number_pass) != pass) {
                        continue;
                    }
                    kmer seq(str2num(sequence.substr(0, k))), rcSeq(rcb(seq)), canon(min_k(seq, rcSeq));
                    canon = (min_k(seq, rcSeq));
                    positions_mutex[indice % 4096].lock();
                    if (min2kmer2context[indice].count(canon) == 0) {
                        min2kmer2context[indice][canon] = { false, {} };
                        // #pragma omp atomic
                        inserted_elements++;
                    }
                    min2kmer2context[indice][canon].count.push_back({ color, coverage });
                    uint64_t sks = sequence.size();
                    for (uint i(0); i + k < sks; ++i) {
                        updateK(seq, sequence[i + k]);
                        updateRCK(rcSeq, sequence[i + k]);
                        canon = (min_k(seq, rcSeq));
                        if (min2kmer2context[indice].count(canon) == 0) {
                            min2kmer2context[indice][canon] = { false, {} };
                            // #pragma omp atomic
                            inserted_elements++;
                        }
                        min2kmer2context[indice][canon].count.push_back({ color, coverage });
                    }
                    if (inserted_elements > 100000000) {
                        // #pragma omp critical
                        toobig = true;
                    }
                    positions_mutex[indice % 4096].unlock();
                }
            }
        }
        if (not toobig) {

            get_monocolor_minitigs_mem(min2kmer2context, out, minimizers, number_color, colormode);
        }
    }
    number_pass *= 2;
    if (toobig) {
        merge_super_buckets_mem(input_file, number_color, out_name, number_pass, colormode);
    }
    //~ out->close();
    delete out;
}

string nucleotides("ACGT");

// look for a next kmer to be merged in the current monotig
//NB updateK : build the next kmer by rolling to the next nucleotide
kmer kmer_Set_Light::select_good_successor(const robin_hood::unordered_node_map<kmer, kmer_context>& kmer2context, const kmer& start)
{
    kmer canon = canonize(start, k);
    if (kmer2context.count(canon) == 0) {
        return -1;
    }
    kmer_context kc(kmer2context.at(canon));
    for (uint64_t i(0); i < 4; ++i) {
        kmer target = start;
        updateK(target, nucleotides[i]);
        kmer targetc = canonize(target, k);
        if (kmer2context.count(targetc) != 0) {
            if (kmer2context.at(targetc).isdump) {
                continue;
            }
            if (kmer2context.at(targetc).count == kc.count) {
                return target;
            }
        }
    }
    return -1;
}

void kmer_Set_Light::write_buffer_count(vector<string>& buffers, zstr::ofstream* out, vector<uint16_t>& headerV, string& seq2dump, int32_t minimi)
{
    string tmp_buffer("");
    uint n = headerV.size() * 2;
    uint nn(n + 4096);
    unsigned char comp[nn];
    unsigned char* in = (unsigned char*)&headerV[0];
    unsigned compr_header_size = trlec(in, n, comp);
    tmp_buffer += ">";
    tmp_buffer.append(reinterpret_cast<char*>(&minimi), sizeof(int32_t)); //minimizer
    tmp_buffer.append(reinterpret_cast<char*>(&compr_header_size), sizeof(unsigned)); //size of compressed rle colors
    tmp_buffer.append((char*)&comp[0], compr_header_size); //rle colors
    tmp_buffer += "\n" + seq2dump + "\n"; // sequence
    buffers[minimi % number_superbuckets.value()] += tmp_buffer;
#pragma omp critical(monocolorFile)
    if (buffers[minimi % number_superbuckets.value()].size() > 8000) {
        {
            *out << buffers[minimi % number_superbuckets.value()];
        }
        buffers[minimi % number_superbuckets.value()].clear();
    }
}

void kmer_Set_Light::write_buffer_color(vector<string>& buffers, zstr::ofstream* out, vector<uint8_t>& headerV, string& seq2dump, int32_t minimi)
{
    string tmp_buffer("");
    uint n = headerV.size();
    uint nn(n + 4096);
    unsigned char comp[nn];
    unsigned char* in = (unsigned char*)&headerV[0];
    unsigned compr_header_size = trlec(in, n, comp);
    //~ unsigned char *decoded;
    //~ decoded = new unsigned char [nn + 4096];
    //~ unsigned char *comp_written0;
    //~ comp_written0 = new unsigned char[compr_header_size];
    //~ unsigned sz = trled(comp_written0, compr_header_size, decoded, headerV.size());
    tmp_buffer += ">";
    tmp_buffer.append(reinterpret_cast<char*>(&minimi), sizeof(int32_t)); //minimizer
    tmp_buffer.append(reinterpret_cast<char*>(&compr_header_size), sizeof(unsigned)); //size of compressed rle colors
    tmp_buffer.append((char*)&comp[0], compr_header_size); //rle colors
    tmp_buffer += "\n" + seq2dump + "\n"; // sequence
    buffers[minimi % number_superbuckets.value()] += tmp_buffer;
#pragma omp critical(monocolorFile)
    if (buffers[minimi % number_superbuckets.value()].size() > 8000) {
        {
            *out << buffers[minimi % number_superbuckets.value()];
        }
        buffers[minimi % number_superbuckets.value()].clear();
    }
}

//writes the final monotigs
void kmer_Set_Light::get_monocolor_minitigs_mem(vector<robin_hood::unordered_node_map<kmer, kmer_context>>& min2kmer2context, zstr::ofstream* out, const vector<int32_t>& mini, uint64_t number_color, int colormode)
{
    // #pragma omp parallel num_threads(coreNumber)
    {
        vector<string> buffers(mini.size(), "");
        string sequence, seq2dump, compact;
        uint64_t ms = min2kmer2context.size();
        // #pragma omp for schedule(static,ms/coreNumber)
        for (uint i_set = (0); i_set < ms; i_set++) {
            for (auto& it : min2kmer2context[i_set]) {
                sort(it.second.count.begin(), it.second.count.end());
            }
            for (auto& it : min2kmer2context[i_set]) {
                if (not it.second.isdump) {
                    it.second.isdump = true;
                    auto colorV2dump(it.second.count);
                    kmer start = it.first;
                    seq2dump = kmer2str(start);
                    for (uint64_t step(0); step < 2; step++) {
                        while (true) {
                            kmer next = select_good_successor(min2kmer2context[i_set], start);
                            if (next == (kmer)-1) {
                                break;
                            }
                            compact = compaction(seq2dump, kmer2str(next), false);
                            if (compact.empty()) {
                                break;
                            }
                            min2kmer2context[i_set].at(canonize(next, k)).isdump = true;
                            seq2dump = compact;
                            start = next;
                        }
                        start = rcb(it.first);
                        seq2dump = revComp(seq2dump);
                    }
// mandatory updates for the index
#pragma omp atomic
                    nuc_minimizer[mini[i_set]] += seq2dump.size();
#pragma omp atomic
                    all_buckets[mini[i_set]].skmer_number++;
#pragma omp atomic
                    all_mphf[mini[i_set] / number_bucket_per_mphf].mphf_size += seq2dump.size() - k + 1;
                    all_mphf[mini[i_set] / number_bucket_per_mphf].empty = false;
#pragma omp atomic
                    total_nuc_number += seq2dump.size();
                    // compress colors in monotigs headers
                    vector<uint16_t> headerV_count;
                    vector<uint8_t> headerV_color;
                    if (mini[i_set] >= 0) {
                        if (color_mode == 0) // colors
                        {
                            headerV_color = get_color_vector(colorV2dump, number_color);
                            write_buffer_color(buffers, out, headerV_color, seq2dump, mini[i_set]);
                        } else // counts or others
                        {
                            headerV_count = get_count_vector(colorV2dump, number_color);
                            write_buffer_count(buffers, out, headerV_count, seq2dump, mini[i_set]);
                        }
                    }
                }
            }
        }
#pragma omp critical(monocolorFile)
        for (uint i_set = (0); i_set < ms; i_set++) {
            if (mini[i_set] >= 0) {
                *out << buffers[mini[i_set] % number_superbuckets.value()] << flush;
            }
        }
        //~ buffer.clear();
    }
}

uint16_t kmer_Set_Light::parseCoverage(const string& str)
{
    if (color_mode == 0) {
        return parseCoverage_bool(str);
    }
    if (color_mode == 1) {
        return parseCoverage_exact(str);
    }
    if (color_mode == 2) {
        return parseCoverage_bin(str);
    }
    return parseCoverage_log2(str);
}

// hashes the kmers from all the input files in super bucket files
void kmer_Set_Light::create_super_buckets_list(const vector<string>& input_files, vector<uint64_t>& kmers_by_file)
{
    kmers_by_file.resize(input_files.size());
    fill(kmers_by_file.begin(), kmers_by_file.end(), 0);
    struct rlimit rl;
    getrlimit(RLIMIT_NOFILE, &rl);
    rl.rlim_cur = number_superbuckets.value() + 10 + coreNumber;

    int res = setrlimit(RLIMIT_NOFILE, &rl);
    if (res == -1) {
        std::cerr << "error setting rlimit to " << number_superbuckets.value() + 10 + coreNumber << std::endl;
        exit(1);
    }
    // atomic<uint64_t> total_nuc_number(0);

    vector<ostream*> out_files;
    for (uint64_t i(0); i < number_superbuckets; ++i) {
        auto out = new zstr::ofstream(working_dir + "_super_k_mers" + to_string(i) + ".gz", zstr::ofstream::app); //todo lz4
        out_files.push_back(out);
    }
    omp_lock_t lock[number_superbuckets.value()];
    for (uint64_t i = 0; i < number_superbuckets; i++) {
        omp_init_lock(&(lock[i]));
    }

#pragma omp parallel num_threads(coreNumber)
    {
#pragma omp for
        for (uint32_t i_file = 0; i_file < input_files.size(); ++i_file) {
            auto inUnitigsread = new zstr::ifstream(input_files[i_file]); //todo lz4
            if (not inUnitigsread->good()) {
                delete inUnitigsread;
                continue;
            }

            string ref, useless;
            vector<string> buffer(number_superbuckets.value());
            minimizer_type old_minimizer, minimizer;
            while (not inUnitigsread->eof()) {
                uint64_t occurrence_kmer = 0;
                ref = useless = "";
                getline(*inUnitigsread, useless);
                getline(*inUnitigsread, ref);
                if (ref.size() < k) {
                    ref = "";
                } else {
                    occurrence_kmer = ref.size() - k + 1;
#pragma omp atomic
                    read_kmer += ref.size() - k + 1;
                }
                //FOREACH UNITIG
                if (not ref.empty() and not useless.empty()) {
                    old_minimizer = minimizer = minimizer_number.value();
                    uint64_t last_position(0);
                    //FOREACH KMER
                    kmer seq(str2num(ref.substr(0, k)));
                    uint64_t position_min;
                    uint64_t min_seq = (str2num(ref.substr(k - minimizer_size_graph, minimizer_size_graph))), min_rcseq(rcbc(min_seq, minimizer_size_graph)), min_canon(min(min_seq, min_rcseq));
                    minimizer = regular_minimizer_pos(seq, position_min);
                    old_minimizer = minimizer;
                    uint64_t hash_min = unrevhash(minimizer);
                    uint64_t i(0);
                    for (; i + k < ref.size(); ++i) {
                        updateK(seq, ref[i + k]);
                        updateM(min_seq, ref[i + k]);
                        updateRCM(min_rcseq, ref[i + k]);
                        min_canon = (min(min_seq, min_rcseq));
                        uint64_t new_h = unrevhash(min_canon);
                        //THE NEW mmer is a MINIMIZER
                        if (new_h < hash_min) {
                            minimizer = (min_canon);
                            hash_min = new_h;
                            position_min = i + k - minimizer_size_graph + 1;
                        } else {
                            //the previous minimizer is outdated
                            if (i >= position_min) {
                                minimizer = regular_minimizer_pos(seq, position_min);
                                hash_min = unrevhash(minimizer);
                                position_min += (i + 1);
                            }
                        }
                        if (old_minimizer != minimizer) {
                            old_minimizer = (revhash(old_minimizer) % minimizer_number);
                            uint16_t cov(parseCoverage(useless));
                            uint8_t size_sk(i - last_position + k);
                            buffer[old_minimizer / bucket_per_superBuckets.value()].append(reinterpret_cast<char*>(&old_minimizer), 4);
                            buffer[old_minimizer / bucket_per_superBuckets.value()].append(reinterpret_cast<char*>(&i_file), 4);
                            buffer[old_minimizer / bucket_per_superBuckets.value()].append(reinterpret_cast<char*>(&cov), 2);
                            buffer[old_minimizer / bucket_per_superBuckets.value()].append(reinterpret_cast<char*>(&size_sk), 1);
                            buffer[old_minimizer / bucket_per_superBuckets.value()] += ref.substr(last_position, i - last_position + k);
                            if (buffer[old_minimizer / bucket_per_superBuckets.value()].size() > 80000) {
                                omp_set_lock(&(lock[((old_minimizer)) / bucket_per_superBuckets.value()]));
                                *(out_files[((old_minimizer)) / bucket_per_superBuckets.value()]) << buffer[old_minimizer / bucket_per_superBuckets.value()];
                                omp_unset_lock(&(lock[((old_minimizer)) / bucket_per_superBuckets.value()]));
                                buffer[old_minimizer / bucket_per_superBuckets.value()].clear();
                            }
                            last_position = i + 1;
                            old_minimizer = minimizer;
                        }
                    }
                    if (ref.size() - last_position > k - 1) {
                        old_minimizer = (revhash(old_minimizer) % minimizer_number);
                        uint16_t cov(parseCoverage(useless));
                        uint8_t size_sk(ref.size() - last_position);
                        buffer[old_minimizer / bucket_per_superBuckets.value()].append((char*)&old_minimizer, 4);
                        buffer[old_minimizer / bucket_per_superBuckets.value()].append((char*)&i_file, 4);
                        buffer[old_minimizer / bucket_per_superBuckets.value()].append((char*)&cov, 2);
                        buffer[old_minimizer / bucket_per_superBuckets.value()].append((char*)&size_sk, 1);
                        buffer[old_minimizer / bucket_per_superBuckets.value()] += ref.substr(last_position);
                        if (buffer[old_minimizer / bucket_per_superBuckets.value()].size() > 80000) {
                            omp_set_lock(&(lock[((old_minimizer)) / bucket_per_superBuckets.value()]));
                            *(out_files[((old_minimizer)) / bucket_per_superBuckets.value()]) << buffer[old_minimizer / bucket_per_superBuckets.value()];
                            omp_unset_lock(&(lock[((old_minimizer)) / bucket_per_superBuckets.value()]));
                            buffer[old_minimizer / bucket_per_superBuckets.value()].clear();
                        }
                        // count total occurence of this monotig
                        occurrence_kmer *= cov;
                    }
                }
                kmers_by_file[i_file] += occurrence_kmer;
            }
            for (uint64_t i(0); i < number_superbuckets.value(); ++i) {
                if (not buffer[i].empty()) {
                    omp_set_lock(&(lock[i]));
                    *(out_files[i]) << buffer[i];
                    omp_unset_lock(&(lock[i]));
                }
            }
            delete inUnitigsread;
        }
    }
    for (uint64_t i(0); i < number_superbuckets; ++i) {
        *out_files[i] << flush;
        delete (out_files[i]);
    }
}
