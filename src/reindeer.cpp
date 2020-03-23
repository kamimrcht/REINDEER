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

void reindeer_index(uint k, string& fof,  string& color_dump_file, bool record_counts, bool record_reads, string& output, string& color_load_file, uint threads, bool exact, bool do_query_on_disk, bool quantize, bool do_log)
{
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	uint64_t color_number;
	uint16_t c;
	kmer_Set_Light ksl(k,m1,m2,m3,threads,bit);
	int systemRet;
	// BUILD THE INDEX
	uint16_t nb_colors(get_color_number(fof));
	//~ cout << "ttttttttttttttttttt" << nb_colors << endl;
	build_index(k, m1, m2, m3, c, bit, color_load_file, color_dump_file, fof, record_counts, record_reads, color_number, ksl, threads, exact, output, do_query_on_disk, quantize, do_log, nb_colors);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index building and Coloration done total: "<< time_span12.count() << " seconds."<<endl;
}




void reindeer_query(uint k, string& output,string& output_query, bool record_counts, bool record_reads, uint threshold, string& bgreat_paths_fof, string& query, uint threads, bool exact, bool do_query_on_disk)
{
	// QUERY //
	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	string fof;
	//~ string fof( getRealPath("graphs.lst", output));
	string color_dump_file("");
	string color_load_file;
	string matrix_name;
	string color_nb_file(getRealPath("reindeer_matrix_eqc_nb_colors", output));
	uint16_t nb_colors ;
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

	string nb_eq_class_file(getRealPath("reindeer_matrix_eqc_nb_class", output));

	uint64_t color_number; //virer
	//~ uint64_t color_number(get_color_number(fof));
	
	vector<unsigned char*> compr_minitig_color;
	vector<unsigned> compr_minitig_color_sizes;

	cout << "\nLoading index.." << endl;
	long eq_class_nb(0);
	kmer_Set_Light* ksl = load_rle_index(k, color_load_file, color_dump_file, fof, record_counts, record_reads, color_number, threads, exact, output, compr_minitig_color, compr_minitig_color_sizes, do_query_on_disk, nb_eq_class_file, color_nb_file, eq_class_nb, nb_colors);
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Index loading total: "<< time_span12.count() << " seconds."<<endl;
	//~ cout << nb_colors << endl;
	cout << "\nComputing query..." << endl;
	high_resolution_clock::time_point tnew = high_resolution_clock::now();
	perform_query(*ksl, nb_colors, k, record_counts,  record_reads,  threshold, bgreat_paths_fof, query, output_query, threads, exact, compr_minitig_color, compr_minitig_color_sizes, do_query_on_disk, matrix_name, eq_class_nb);
	high_resolution_clock::time_point tnew2 = high_resolution_clock::now();
	duration<double> time_spannew2 = duration_cast<duration<double>>(tnew2 - tnew);
	cout<<"Querying sequence took "<< time_spannew2.count() << " seconds."<<endl;
}

