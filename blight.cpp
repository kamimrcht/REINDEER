#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>

#include "ksl.h"



using namespace std;
using namespace chrono;




int main(int argc, char ** argv){
	cout<<"START"<<endl;
	kmer_Set_Light ksl(31,8,4,8);

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	//~ ksl.create_super_buckets("lambda_virus.unitigs.fa");
	//~ ksl.create_super_buckets("swag");
	//~ ksl.create_super_buckets("ecoli.31.unitigs.fa");
	ksl.create_super_buckets("cel.31.unitigs.fa");
	high_resolution_clock::time_point t12 = high_resolution_clock::now();
	duration<double> time_span12 = duration_cast<duration<double>>(t12 - t1);
	cout<<"Super bucket created: "<< time_span12.count() << " seconds."<<endl;
	ksl.read_super_buckets("_out");
	high_resolution_clock::time_point t13 = high_resolution_clock::now();
	duration<double> time_span13 = duration_cast<duration<double>>(t13 - t12);
	cout<<"Bucket loaded: "<< time_span13.count() << " seconds."<<endl;
	ksl.create_mphf();
	high_resolution_clock::time_point t14 = high_resolution_clock::now();
	duration<double> time_span14 = duration_cast<duration<double>>(t14 - t13);
	cout<<"Mphf constructed: "<< time_span14.count() << " seconds."<<endl;
	ksl.fill_positions();
	high_resolution_clock::time_point t15 = high_resolution_clock::now();
	duration<double> time_span15 = duration_cast<duration<double>>(t15- t14);
	cout<<"Position  filled: "<< time_span15.count() << " seconds."<<endl;

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "The whole indexing took me " << time_span.count() << " seconds."<< endl;


	//~ ksl.multiple_query("lambda_virus.unitigs.fa");
	//~ ksl.multiple_query("ecoli.31.unitigs.fa");
	ksl.multiple_query("cel.31.unitigs.fa");
	//~ ksl.multiple_query("swage");

	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
	cout << "The whole query took me " << time_span2.count() << " seconds.";
	cout << endl;



	return 0;
}

