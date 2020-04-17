#include "reindeer.hpp"


using namespace std;
using namespace chrono;



// MPHF options
uint m1(10);
uint m2(10);
uint m3(5);
uint c(1);
uint bit(0);
uint ex(0);

void reindeer_index(uint k, string& fof,  string& color_dump_file, bool record_counts, string& output, string& color_load_file, uint threads,  bool do_query_on_disk, bool quantize, bool do_log)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint16_t c;
	kmer_Set_Light ksl(k,m1,m2,m3,threads,bit);
	int systemRet;
	// BUILD THE INDEX
	uint64_t nb_colors(get_color_number(fof));
	build_index(k, m1, m2, m3, c, bit, color_load_file, color_dump_file, fof, record_counts,  ksl, threads,  output, do_query_on_disk, quantize, do_log, nb_colors);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index building and Coloration done total: "<< time_span12.count() << " seconds."<<endl;
}




void reindeer_query(string& output,string& output_query, bool record_counts, uint threshold,  string& query, uint threads,  bool do_query_on_disk)
{
	// QUERY //
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	string fof;
	string color_dump_file("");
	string color_load_file;
	string matrix_name;
	uint64_t nb_colors, nb_monotig ;
	uint k, record_option;
	long eq_class_nb;
	read_info(k, nb_monotig, eq_class_nb, nb_colors, record_option, matrix_name);
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
	
	vector<unsigned char*> compr_monotig_color;
	vector<unsigned> compr_monotig_color_sizes;
	cout << "\n#Loading index..." << endl;
	//~ long eq_class_nb(0);
	bool quantize, log;
	kmer_Set_Light* ksl = load_rle_index(k, color_load_file, color_dump_file, fof, record_counts,  threads,  output, compr_monotig_color, compr_monotig_color_sizes, do_query_on_disk,  eq_class_nb, nb_colors,   quantize,  log);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"#Index loading total: "<< time_span12.count() << " seconds."<<endl;
	cout << "\n#Computing query..." << endl;
	high_resolution_clock::time_point tnew = high_resolution_clock::now();
	perform_query(*ksl, nb_colors, k, record_counts,  threshold,  query, output_query, threads,  compr_monotig_color, compr_monotig_color_sizes, do_query_on_disk, matrix_name, eq_class_nb, nb_monotig);
	high_resolution_clock::time_point tnew2 = high_resolution_clock::now();
	duration<double> time_spannew2 = duration_cast<duration<double>>(tnew2 - tnew);
	cout<<"#Querying sequences took "<< time_spannew2.count() << " seconds in total."<<endl;
}

