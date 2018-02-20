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



int main(int argc, char ** argv){
	cout<<"go"<<endl;
	kmer_Set_Light ksl(6,4, 2,2);
	ksl.create_super_buckets("lambda_virus.unitigs.fa");
	//~ ksl.create_super_buckets("swag");
	cout<<"super bucket created"<<endl;
	ksl.read_super_buckets("_out");
	cout<<"bucket loaded"<<endl;
	ksl.create_mphf();
	cout<<"mphf constructed"<<endl;
	ksl.fill_positions();
	cout<<"position  filled"<<endl;
	ksl.multiple_query("GTCAGCGAGGACGGGTATCCGGTTTCCGTCTTCACGGACTTCGTTGCTTTCCAGTTTAGCAATACGCTTACTCCCATCCGAGATAACACCTTCGTAATACTCACGCTGCTCGTTGAGTTTTGATTTTGCTGTTTCAAGCTCAACACGCAGTTTCCCTACTGTTAGCGCAATATCCTCGTTCTCCTGGTCGCGGCGTTTGATGT");
	//~ ksl.multiple_query("CATGCATGCTAGCTGCATGCTAGCTGCATGCATGCTAGCTGCA");
	//~ ksl.multiple_query("ATTCGATGC");
	return 0;
}

