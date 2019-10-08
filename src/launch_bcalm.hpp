#ifndef BCA
#define BCA
#include <stdio.h>
#include <stdlib.h>
#include "utils.hpp"

using namespace std;

void get_list_graphs_fof(string& fof, vector<string>& list_graphs)
{
	ifstream file(fof);
	string sample;
	while(not file.eof())
	{
		getline(file, sample);
		if (sample.empty()){break;}
		list_graphs.push_back(sample);
	}
}

void bcalm_launcher_single(string& input_fof, uint k, uint t, string& output_bcalm)
{
	int systemRet;
	systemRet=system(("mkdir -p "+ output_bcalm).c_str());
	if(systemRet == -1){}
	string input, output;
	vector<string> list_graphs;
	get_list_graphs_fof(input_fof, list_graphs);
	ofstream file_list_graphs( output_bcalm + "/graphs.lst");
	for (uint i(0); i < list_graphs.size(); ++i)
	{
		input = list_graphs[i];
		output = split(list_graphs[i], '.')[0];
		systemRet = system(("./bin/bcalm -in " + input + " -kmer-size " + to_string(k) + " -abundance-min 2  -nb-cores " + to_string(t) + " -out " + output).c_str());
		if(systemRet == -1){
			// failed
		}
		systemRet = system(("mv " + output + ".unitigs.fa " + output_bcalm).c_str());
		if(systemRet == -1){
			// failed
		}
		file_list_graphs << output_bcalm + "/" + output + ".unitigs.fa" << endl;
	}
}


void bcalm_launcher_union(string& graph_list, uint k, uint t, string& output_bcalm)
{

	int systemRet;
	systemRet = system(("./bin/bcalm -in " + graph_list + " -kmer-size " + to_string(k) + " -abundance-min  1 -nb-cores " + to_string(t) + " -out union_graph").c_str());
	if(systemRet == -1){}
	systemRet = system(("mv union_graph.unitigs.fa " + output_bcalm).c_str());
	if(systemRet == -1){}
}

#endif
