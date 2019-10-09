#ifndef BCA
#define BCA
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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

string bcalm_launcher_single(string& input_fof, uint k, uint t, string& main_output, string& output_bcalm)
{
	int systemRet;
	string outputDir(main_output + "/" +output_bcalm);
	systemRet=system(("mkdir -p "+ outputDir).c_str());
	//~ if(systemRet == -1){}
	string input, output;
	vector<string> list_graphs;
	get_list_graphs_fof(input_fof, list_graphs);
	ofstream file_list_graphs( main_output + "/graphs.lst");
	string cmd(""), filename("");
	char* f;
	for (uint i(0); i < list_graphs.size(); ++i)
	{
		input = list_graphs[i];
		output = split(split(list_graphs[i], '.')[0], '/').back();
		cmd = "./bin/bcalm -in " + input + " -kmer-size " + to_string(k) + " -abundance-min 2  -nb-cores " + to_string(t) + " -out " + output;
		cout << cmd << endl;
		systemRet = system(cmd.c_str());
		if(systemRet == -1){
			// failed
		}
		systemRet = system(("mv " + output + ".unitigs.fa " +outputDir).c_str());
		if(systemRet == -1){
			// failed
		}
		filename = outputDir + "/" +output + ".unitigs.fa";
		f = &(filename[0]);
		char* full_path = getcwd(f, filename.size());
		file_list_graphs << string(full_path) << "/" << outputDir << "/" << output + ".unitigs.fa" << endl;
		//~ cout << string(full_path) << "/" << outputDir << "/" << output + ".unitigs.fa" << endl;
	}
	return main_output + "/graphs.lst";
}


string bcalm_launcher_union(string& graph_list, uint k, uint t, string& main_output, string& output_bcalm)
{

	int systemRet;
	string outputDir(main_output + "/" +output_bcalm);
	systemRet=system(("mkdir -p "+ outputDir).c_str());
	string cmd("./bin/bcalm -in " + graph_list + " -kmer-size " + to_string(k) + " -abundance-min  1 -nb-cores " + to_string(t) + " -out union_graph");
	cout << cmd << endl;
	systemRet = system((cmd).c_str());
	if(systemRet == -1){}
	systemRet = system(("mv union_graph.unitigs.fa " + outputDir).c_str());
	if(systemRet == -1){}
	return outputDir + "/union_graph.unitigs.fa";
}

void bcalm_cleanup()
{
	int systemRet;
	systemRet=system("rm -f *glue* *.h5");
}


string getRealPaths(string& fof, string& main_output)
{
	string input,names;
	vector<string> list_graphs;
	get_list_graphs_fof(fof, list_graphs);
	ofstream file_list_graphs( main_output + "/graphs.lst");
	string cmd(""), filename("");
	char* f;
	for (uint i(0); i < list_graphs.size(); ++i)
	{
		input = list_graphs[i];
		names = list_graphs[i];
		f = &(input[0]);
		char* full_path = getcwd(f, names.size());
		string fp(full_path);
		file_list_graphs << fp << "/" << names << endl;

	}
	return main_output + "/graphs.lst";
}

#endif
