
#include "build_index.hpp"
#include "reindeer.hpp"

using namespace chrono;
using namespace std;

// read colors from bcalm headers
vector<uint8_t> get_colors_monotigs(string& line)
{
    vector<uint8_t> colors;
    vector<string> colors_monotig = split_utils(line, ':');
    uint value;
    for (uint c(1); c < colors_monotig.size(); ++c) // convert the bit string to a bit vector
    {
        value = stoi(colors_monotig[c]);
        if (value > 0)
            colors.push_back(1);
        else
            colors.push_back(0);
    }
    return colors;
}

// read counts from bcalm headers
vector<uint16_t> get_counts_monotigs(string& line)
{
    vector<uint16_t> counts;
    uint pred(0);
    for (uint c(1); c < line.size(); ++c) // convert the bit string to a bit vector
    {
        if (line[c] == ':') {
            if (pred != 0) {
                counts.push_back(stoi(line.substr(pred, c - pred)));
            }
            pred = c + 1;
        }
    }
    counts.push_back(stoi(line.substr(pred)));
    return counts;
}

void get_header_monotig_file(zstr::ifstream& in, string& header)
{
    unsigned compressed_header_size;
    header = "";
    unsigned char* trash; //todo optim
    trash = new unsigned char[5];
    in.read((char*)trash, 5); // ">" and minimizer
    delete[] trash;
    in.read(reinterpret_cast<char*>(&compressed_header_size), sizeof(unsigned)); // size of colors/counts with rle
    header.append(reinterpret_cast<char*>(&compressed_header_size), sizeof(unsigned));
    if (header.capacity() < sizeof(unsigned) + compressed_header_size)
        header.reserve(header.size() + sizeof(unsigned) + compressed_header_size + 1024); //todo ok ?
    header.resize(sizeof(unsigned) + compressed_header_size);
    in.read((char*)(&header[sizeof(unsigned)]), compressed_header_size); //rle
    char c = in.get(); // \n
}

// dispatch count vectors in files. Similar counts go in similar files
template <class T>
void Reindeer_Index<T>::write_matrix_in_bucket_files(kmer_Set_Light* ksl)
{
    //create bucket files for partitionning the compressed counts -> finding eq classes
    vector<ofstream*> all_files;
    for (uint i(0); i < 1000; ++i) {
        ofstream* out = new ofstream(reindeer_index_files + "/matrix_bucket_" + to_string(i) + ".gz");
        all_files.push_back(out);
    }

    string output_file_name;
    vector<string> monotig_files_names;
    string monotig_folder(monotig_files + "/");
    get_all_blout(monotig_folder, monotig_files_names); //todo .lz4
    ofstream out_info(matrix_eqc_info_file);
    mutex mm;
    uint nb_treated_monotigs(0);
    uint i;
#pragma omp parallel for num_threads(threads)
    for (i = 0; i < monotig_files_names.size(); ++i) //loop on _blout files
    {
        string fname;
        fname = monotig_files_names[i];
        if (fname != "." and fname != "..") {
            fname = monotig_files + "/" + fname;
            zstr::ifstream monotigs_file(fname);
            if (not exists_test(fname)) {
                cerr << "File problem\n";
                continue;
            }
            vector<int64_t> monotig_id;
            vector<uint16_t> counts;
            vector<uint8_t> colors;
            string header(5, ' '), monotig, buffer;
            header.reserve(4096);
            // record count vector for each monotig at index given by the mphf
            monotigs_file.peek();
            while (not monotigs_file.eof()) {
                get_header_monotig_file(monotigs_file, header);
                getline(monotigs_file, monotig);
                if (monotig.empty() or header.empty()) {
                    continue;
                }
                monotig_id.clear();
                if (monotig[0] == 'A' or monotig[0] == 'C' or monotig[0] == 'G' or monotig[0] == 'T') {
#pragma omp atomic
                    ++nb_monotig;
                    // get index from MPHF
                    monotig_id = ksl->get_rank_query(monotig.substr(0, k)); // all kmers have the same id so we only query one
                    if ((not monotig_id.empty()) and monotig_id.back() >= 0) {
                        mm.lock();
                        dump_compressed_vector_bucket(monotig_id.back(), all_files, header);
                        nb_treated_monotigs++;
                        mm.unlock();
                    }
                }
                monotigs_file.peek();
            }
            if (dele_monotig_file)
                remove(fname.c_str());
        }
    }
    for (uint i(0); i < all_files.size(); ++i) {
        all_files[i]->flush();
        all_files[i]->close();
    }
    out_info.write(reinterpret_cast<char*>(&nb_monotig), sizeof(uint64_t));
    out_info.write(reinterpret_cast<char*>(&k), sizeof(uint)); // in info: 1/nb monotigs, 2/k 3/record option 4/nb eq classes, 5/nb colors
    uint val(1);
    if (!record_counts) {
        val = 0;
    } else {
        if (quantize) {
            val = 2;
        } else if (do_log) {
            val = 3;
        } else {
            val = 1;
        }
    }
    out_info.write(reinterpret_cast<char*>(&val), sizeof(uint));
    // compute final equivalence class and write them
    write_eq_class_matrix(all_files, &out_info);
    // remove bucket files
    for (uint i(0); i < all_files.size(); ++i) {
        string name(reindeer_index_files + "/matrix_bucket_" + to_string(i) + ".gz");
        remove(&name[0]);
        delete all_files[i];
    }
    out_info.write(reinterpret_cast<char*>(&do_query_on_disk), sizeof(bool));
    out_info.close();
}

// color using monotig file: either build and dump the color matrix during the index construction, or load it during the query
template <class T>
void Reindeer_Index<T>::do_coloring(kmer_Set_Light* ksl)
{
    vector<string> file_names;
    if (not color_load_file.empty()) //query
    {
        if (not do_query_on_disk) {
            uint64_t color_number;
            uint64_t monotig_number;
            compressed_monotig_color = load_compressed_vectors(color_load_file, compressed_monotig_color_sizes, color_number, monotig_number, nb_eq_class);
        }
    } else //indexing
    {
        if (exists_test(fof)) {
            ifstream fofin(fof);
            string file_name;
            while (not fofin.eof()) {
                getline(fofin, file_name);
                if (not file_name.empty()) {
                    if (exists_test(file_name)) {
                        file_names.push_back(file_name);
                    } else {
                        cerr << "[ERROR] " << file_name << " is not here" << endl;
                    }
                }
            }
        } else {
            cerr << "[ERROR] File of file problem" << endl;
        }

        write_matrix_in_bucket_files(ksl);
    }
}

// load dumped index(+colors)
template <class T>
kmer_Set_Light* Reindeer_Index<T>::load_rle_index()
{
    kmer_Set_Light* ksl = new kmer_Set_Light(reindeer_index_files + "/reindeer_index.gz");
    do_coloring(ksl);
    if (dele_monotig_file) {
        string cmd("rm -rf " + reindeer_index_files + "/monotig_files");
        int sysRet(system(cmd.c_str()));
    }
    return ksl;
}

// build index from new file
template <class T>
void Reindeer_Index<T>::build_index(kmer_Set_Light* ksl)
{
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    bool dont_dump(false);
    int color_mode;
    if (record_counts)
        if (quantize)
            color_mode = 2;
        else if (do_log)
            color_mode = 3;
        else
            color_mode = 1;
    else
        color_mode = 0;
    if (not dirExists(monotig_files)) {
        cout << "#Monotigs and index constuction..." << endl;
        // apply monotig merge (-> MMM) with rule regarding colors or counts
        // color 0, count 1, quantize 2, log 3
        ksl->construct_index_fof(fof, reindeer_index_files, color_mode);
    } else {
        cerr << "[Warning] monotig files (monotig_files) were found in output dir, I will use them and I won't delete them" << endl;
        if (not exists_test(index_file)) {
            // color 0, count 1, quantize 2, log 3
            string m_folder(monotig_files + "/");
            string fof_blout(do_fof(m_folder, reindeer_index_files));
            ksl->construct_index_fof(fof_blout, reindeer_index_files, color_mode); //todo + .lz4
            //todo could be optimized: blight needs to update some initialized variables and could take these directly as input super kmers
            remove((reindeer_index_files + "/fof_blout").c_str());
        } else {
            dont_dump = true;
            cerr << "[Warning] index file (reindeer_index.gz) was found in output dir, I will use it and I won't delete it" << endl;
            ksl = new kmer_Set_Light(index_file);
        }
    }

    if (!dont_dump) {
        cout << "#Dumping index..." << endl;
        ksl->dump_disk(index_file);
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        duration<double> time_span12 = duration_cast<duration<double>>(t2 - t1);
        cout << "Index written on disk: " << time_span12.count() << " seconds." << endl;
    }
    cout << "#Building colors and equivalence classes matrix to be written on disk..." << endl;
    high_resolution_clock::time_point t3 = high_resolution_clock::now();
    do_coloring(ksl);
    high_resolution_clock::time_point t4 = high_resolution_clock::now();
    duration<double> time_span34 = duration_cast<duration<double>>(t4 - t3);
    cout << "Matrix done: " << time_span34.count() << " seconds." << endl;
    if (dele_monotig_file) {
        string cmd("rm -r " + monotig_files);
        int sysRet(system(cmd.c_str()));
    }
}
