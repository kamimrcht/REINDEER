#include "reindeer.hpp"

using namespace std;
using namespace chrono;

// constructor for index construction
template <class T>
Reindeer_Index<T>::Reindeer_Index(uint pk, string& pfof, bool precord_counts, string& preindeer_index_files, uint pthreads, bool pdo_query_on_disk, bool pquantize, bool pdo_log, uint pm1, uint pm3, bool dele_tmp, bool poutput_monotigs)
{
    // MPHF options
    m2 = 10;
    bit = 0;
    ex = 0;
    m1 = pm1;
    m3 = pm3;

    k = pk; //size of k-mers

    // options
    threads = pthreads; // number of threads
    record_counts = precord_counts; // if true, record counts else presence/absence //todo replace by color mode
    do_query_on_disk = pdo_query_on_disk; // if true the full index is dumped on the disk, else it is rebuilt and in ram
    output_monotigs = poutput_monotigs;

    // other counts modes
    quantize = pquantize; // record quantization
    do_log = pdo_log; // record log

    // color parameters
    nb_eq_class = 0;

    //files
    reindeer_index_files = preindeer_index_files;
    color_dump_file = reindeer_index_files + "/reindeer_matrix";
    monotig_files = reindeer_index_files + "/monotig_files";
    index_file = reindeer_index_files + "/reindeer_index.gz";
    matrix_eqc_info_file = reindeer_index_files + "/reindeer_matrix_eqc_info.txt";
    matrix_eqc_file = reindeer_index_files + "/reindeer_matrix_eqc";
    matrix_eqc_position_file = reindeer_index_files + "/reindeer_matrix_eqc_position";
    dele_monotig_file = dele_tmp;

    fof = pfof;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    ksl = new kmer_Set_Light(k, m1, m2, m3, threads, bit);
    int systemRet;
    // BUILD THE INDEX
    nb_colors = get_color_number(fof);
    build_index();
    high_resolution_clock::time_point t12 = high_resolution_clock::now();
    duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
    cout << "Index building and Coloration done total: " << time_span12.count() << " seconds." << endl;
    uint64_t mem(getMemorySelfMaxUsed());
    cout << "Max Memory used: " << mem << endl;
}

template <class T>
void Reindeer_Index<T>::read_info()
{
    ifstream info_file(matrix_name + "_info"); //todo change
    if (!info_file.is_open()) {
        cout << "Can't open an index file" << endl;
        exit(1);
    }
    uint record_option;
    //get nb of monotigs, nb of eq_classes, nb_of colors

    string line;
    getline(info_file, line);
    size_t colon_pos = line.find(':');
    nb_monotig = stoull(line.substr(colon_pos + 1));

    getline(info_file, line);
    colon_pos = line.find(':');
    k = stoul(line.substr(colon_pos + 1));

    getline(info_file, line);
    colon_pos = line.find(':');
    record_option = stoul(line.substr(colon_pos + 1));

    getline(info_file, line);
    colon_pos = line.find(':');
    nb_eq_class = stol(line.substr(colon_pos + 1));

    getline(info_file, line);
    colon_pos = line.find(':');
    nb_colors = stoull(line.substr(colon_pos + 1));

    getline(info_file, line);
    colon_pos = line.find(':');
    string do_query_on_disk_string = line.substr(colon_pos + 1);
    do_query_on_disk = (do_query_on_disk_string != "0");

    kmers_by_file.resize(nb_colors);
    fill(kmers_by_file.begin(), kmers_by_file.end(), make_pair("",0));
    for (int i = 0; i < nb_colors; i++){
        getline(info_file, line);
        colon_pos = line.find(":");
        kmers_by_file[i]=(make_pair(line.substr(0, colon_pos),stoull(line.substr(colon_pos + 1))));
    }

    if (record_option == 1) {
        record_counts = true;
    }
}

// constructor for query
template <class T>
Reindeer_Index<T>::Reindeer_Index(string& poutput, string& poutput_query, uint pthreads, bool dele_tmp)
{
    color_load_file = poutput;
    reindeer_index_files = poutput;
    output = poutput_query;
    threads = pthreads;
    // QUERY //
    //check if loading directory exists and all reindeer files are present
    color_load_file = getRealPath("reindeer_matrix_eqc", poutput);
    matrix_name = color_load_file;
    matrix_eqc_info_file = matrix_name + "_info";
    matrix_eqc_position_file = matrix_name + "_position";
    dele_monotig_file = dele_tmp;
}

template <class T>
void Reindeer_Index<T>::load_index(){
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    read_info();
    cout << "\n#Loading index..." << endl;
    std::ofstream index_loading_semaphore(reindeer_index_files + "/index_loading"); // begin semaphore
    //~ long eq_class_nb(0);
    bool quantize = false, log = false;
    this->ksl = load_rle_index();
    this->position_in_file = get_position_vector_query_disk();
    high_resolution_clock::time_point t12 = high_resolution_clock::now();
    duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
    cout << "Index loading total: " << time_span12.count() << " seconds." << endl;
    std::remove(string(output + "/index_loading").c_str()); // end semaphore
}

template <class T>
void Reindeer_Index<T>::querying(string query, uint threshold, string& output_format){
    this->output_format = output_format;
    this->threshold = threshold;
    cout << "\n#Computing query..." << endl;
    high_resolution_clock::time_point tnew = high_resolution_clock::now();
    perform_query(query);
    high_resolution_clock::time_point tnew2 = high_resolution_clock::now();
    duration<double> time_spannew2 = duration_cast<duration<double>>(tnew2 - tnew);
    cout << "Querying sequences took " << time_spannew2.count() << " seconds in total." << endl;
    uint64_t mem(getMemorySelfMaxUsed());
    cout << "Max Memory used: " << mem << endl;
    //~ return 0;
}

template <class T>
void Reindeer_Index<T>::print_Reindeer()
{
    cout << "m2 " << m2 << endl;
    cout << "bit " << bit << endl;
    cout << "ex " << ex << endl;
    cout << "m1 " << m1 << endl;
    cout << "m3 " << m3 << endl;
    cout << "k " << 31 << endl;
    cout << "threads " << threads << endl; // number of threads
    cout << "record_counts " << record_counts << endl; // if true, record counts else presence/absence //todo replace by color mode
    cout << "color_mode " << color_mode << endl;
    cout << "do_query_on_disk " << do_query_on_disk << endl;
    cout << "quantize " << quantize << endl; // record quantization
    cout << "do_log " << do_log << endl; // record log
    cout << "output_format" << output_format << endl; // output_format

    // color variables
    cout << "nb_colors " << nb_colors << endl; // number of samples
    cout << "nb_monotig " << nb_monotig << endl;
    cout << "nb_eq_class " << nb_eq_class << endl;

    //files
    cout << "matrix_eqc_info_file" << matrix_eqc_info_file << endl;
    cout << "matrix_eqc_position_file" << matrix_eqc_position_file << endl;
    cout << "matrix_eqc_file" << matrix_eqc_file << endl;
    cout << "color_load_file " << color_load_file << endl;
    cout << "color_dump_file " << color_dump_file << endl;
    cout << "fof " << fof << endl;
    cout << "reindeer_index_files " << reindeer_index_files << endl;
    cout << "index_file " << index_file << endl;
    cout << "output " << output << endl;
    cout << "monotig_files " << monotig_files << endl;
    cout << "dumped_index " << dumped_index << endl;
    cout << "matrix_name " << matrix_name << endl;
}
