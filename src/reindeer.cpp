#include "reindeer.hpp"


using namespace std;
using namespace chrono;

void build_reindeer_index(reindeer_index& index_values)
{
	// clock
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	
	// build mphf (BLight functions)
	kmer_Set_Light ksl(index_values.k, index_values.m1, index_values.m2, index_values.m3, index_values.threads, index_values.bit);
	// build index
	build_index(index_values, &ksl);
	
	// log
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index building and Coloration done total: "<< time_span12.count() << " seconds."<<endl;
	uint64_t mem(getMemorySelfMaxUsed());
	cout<<"Max Memory used: "<< mem  << endl;
}




uint reindeer_query(reindeer_index& index_values,string& output_query, uint threshold,  string& query)
{
	// QUERY //
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	string matrix_name;

	//check if loading directory exists and all reindeer files are present
	if (dirExists(index_values.output))
	{
		if ( not (exists_test(index_values.output + "/reindeer_matrix_eqc_info") and exists_test(index_values.output + "/reindeer_index.gz") and exists_test(index_values.output + "/reindeer_matrix_eqc_position.gz"))) //todo changer
		{
			cerr << "[ERROR] REINDEER index files are missing. Stopped."<< endl; 
			return 0;
		}
	} else 
	{
		cerr << "[ERROR] REINDEER index directory is missing (or path in -l is wrong). Stopped." << endl;
		return 0;
	}
	cout << "wwwwww" << index_values.do_query_on_disk << endl;
	if (index_values.do_query_on_disk)
	{
		if ( not (exists_test(index_values.output + "/reindeer_matrix_eqc") ))//todo changer + lignes d'apres
		{
			cerr << "[ERROR] REINDEER index files are missing. Stopped."<< endl; 
			return 0;
		}
		else
		{
		index_values.color_load_file = getRealPath("reindeer_matrix_eqc", index_values.output);
		matrix_name = index_values.color_load_file;
		}
	}
	else
	{
		if ( not (exists_test(index_values.output + "/reindeer_matrix_eqc.gz") ))
		{
			cerr << "[ERROR] REINDEER index files are missing. Stopped."<< endl; 
			return 0;
		}
		else
		{
		index_values.color_load_file = getRealPath("reindeer_matrix_eqc.gz", index_values.output);
		size_t wo_ext = index_values.color_load_file.find_last_of("."); 
		matrix_name = index_values.color_load_file.substr(0, wo_ext); 
		}
	}
	uint  record_option;
	read_info(index_values, matrix_name, record_option);
	
	vector<unsigned char*> compr_monotig_color;
	vector<unsigned> compr_monotig_color_sizes;
	cout << "\n#Loading index..." << endl;
    std::ofstream index_loading_semaphore(output_query + "/index_loading"); // begin semaphore
	//~ long eq_class_nb(0);
	bool quantize = false, log = false;
	kmer_Set_Light* ksl = load_rle_index(index_values, compr_monotig_color, compr_monotig_color_sizes);
	vector<long> position_in_file;
	string position_file_name(matrix_name + "_position.gz");
	get_position_vector_query_disk(position_in_file,  position_file_name, index_values.nb_monotig);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index loading total: "<< time_span12.count() << " seconds."<<endl;
    std::remove(string(output_query + "/index_loading").c_str()); // end semaphore
	cout << "\n#Computing query..." << endl;
	high_resolution_clock::time_point tnew = high_resolution_clock::now();
	perform_query(*ksl, index_values.nb_colors, index_values.k, index_values.record_counts,  threshold,  query, output_query, index_values.threads,  compr_monotig_color, compr_monotig_color_sizes, index_values.do_query_on_disk, matrix_name, index_values.nb_eq_class, index_values.nb_monotig, position_in_file);//todo change
	high_resolution_clock::time_point tnew2 = high_resolution_clock::now();
	duration<double> time_spannew2 = duration_cast<duration<double>>(tnew2 - tnew);
	cout<<"Querying sequences took "<< time_spannew2.count() << " seconds in total."<<endl;
	uint64_t mem(getMemorySelfMaxUsed());
	cout<<"Max Memory used: "<< mem  << endl;
	return 0;
}


