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
	kmer_Set_Light ksl(31,4, 2,8);

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	//~ ksl.create_super_buckets("lambda_virus.unitigs.fa");
	ksl.create_super_buckets("ecoli.31.unitigs.fa");
	//~ ksl.create_super_buckets("swag");
	cout<<"Super bucket created"<<endl;
	ksl.read_super_buckets("_out");
	cout<<"Bucket loaded"<<endl;
	ksl.create_mphf();
	cout<<"Mphf constructed"<<endl;
	ksl.fill_positions();
	cout<<"Position  filled"<<endl;

	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "It took me " << time_span.count() << " seconds.";
	cout << endl;

	//~ cin.get();

	ksl.multiple_query("lambda_virus.unitigs.fa");
	//~ ksl.multiple_query("ecoli.31.unitigs.fa");

	high_resolution_clock::time_point t3 = high_resolution_clock::now();
	duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
	cout << "It took me " << time_span2.count() << " seconds.";
	cout << endl;



	return 0;
}

