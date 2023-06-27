#include "reindeer.hpp"


using namespace std;
using namespace chrono;


// MPHF options

uint m2(10);
uint c(1);
uint bit(0);
uint ex(0);

void reindeer_index(uint k, string& fof,  string& color_dump_file, bool record_counts, string& output, bool output_monotigs, string& color_load_file, uint threads,  bool do_query_on_disk, bool quantize, bool do_log, uint m1, uint m3)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	kmer_Set_Light ksl(k,m1,m2,m3,threads,bit);
	int systemRet;
	// BUILD THE INDEX
	uint64_t nb_colors(get_color_number(fof));
	build_index(k, m1, m2, m3, bit, color_load_file, color_dump_file, fof, record_counts, &ksl, threads,  output, output_monotigs, do_query_on_disk, quantize, do_log, nb_colors);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index building and Coloration done total: "<< time_span12.count() << " seconds."<<endl;
	uint64_t mem(getMemorySelfMaxUsed());
	cout<<"Max Memory used: "<< mem  << endl;
}




uint reindeer_query(string& output,string& output_query, uint threshold,  string& query, uint threads,  bool do_query_on_disk)
{
	// QUERY //
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	string fof;
	string color_dump_file("");
	string color_load_file;
	string matrix_name;

	//check if loading directory exists and all reindeer files are present
	if (dirExists(output))
	{
		if ( not (exists_test(output + "/reindeer_matrix_eqc.gz") and exists_test(output + "/reindeer_matrix_eqc_info") and exists_test(output + "/reindeer_index.gz") and exists_test(output + "/reindeer_matrix_eqc_position.gz")))
		{
			cerr << "[ERROR] REINDEER index files are missing. Stopped."<< endl; 
			return 0;
		}
	} else 
	{
		cerr << "[ERROR] REINDEER index directory is missing (or path in -l is wrong). Stopped." << endl;
		return 0;
	}

	if (do_query_on_disk)
	{
		color_load_file = getRealPath("reindeer_matrix_eqc", output);
		matrix_name = color_load_file;
	}
	else
	{
		color_load_file = getRealPath("reindeer_matrix_eqc.gz", output);
		size_t wo_ext = color_load_file.find_last_of("."); 
		matrix_name = color_load_file.substr(0, wo_ext); 
	}
	uint64_t nb_colors, nb_monotig ;
	uint k, record_option;
	long eq_class_nb;
	read_info(k, nb_monotig, eq_class_nb, nb_colors, record_option, matrix_name);
	bool record_counts(0);
	if (record_option == 1)
	{
		record_counts = true;
	}
	vector<unsigned char*> compr_monotig_color;
	vector<unsigned> compr_monotig_color_sizes;
	cout << "\n#Loading index..." << endl;
    std::ofstream index_loading_semaphore(output_query + "/index_loading"); // begin semaphore
	//~ long eq_class_nb(0);
	bool quantize = false, log = false;
	kmer_Set_Light* ksl = load_rle_index(k, color_load_file, color_dump_file, fof, record_counts,  threads,  output, compr_monotig_color, compr_monotig_color_sizes, do_query_on_disk,  eq_class_nb, nb_colors,   quantize,  log);
	vector<long> position_in_file;
	string position_file_name(matrix_name + "_position.gz");
	get_position_vector_query_disk(position_in_file,  position_file_name,nb_monotig);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index loading total: "<< time_span12.count() << " seconds."<<endl;
    std::remove(string(output_query + "/index_loading").c_str()); // end semaphore
	cout << "\n#Computing query..." << endl;
	high_resolution_clock::time_point tnew = high_resolution_clock::now();
	perform_query(*ksl, nb_colors, k, record_counts,  threshold,  query, output_query, threads,  compr_monotig_color, compr_monotig_color_sizes, do_query_on_disk, matrix_name, eq_class_nb, nb_monotig, position_in_file);
	high_resolution_clock::time_point tnew2 = high_resolution_clock::now();
	duration<double> time_spannew2 = duration_cast<duration<double>>(tnew2 - tnew);
	cout<<"Querying sequences took "<< time_spannew2.count() << " seconds in total."<<endl;
	uint64_t mem(getMemorySelfMaxUsed());
	cout<<"Max Memory used: "<< mem  << endl;
	return 0;
}


