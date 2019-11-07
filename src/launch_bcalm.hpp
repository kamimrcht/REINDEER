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

string bcalm_launcher_single(string& input_fof, uint k, uint threads, string& main_output, string& output_bcalm)
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
		//~ output = split(split(list_graphs[i], '.')[0], '/').back();
		output = split_utils(split_utils(list_graphs[i], '/').back(), '.')[0];
		cmd = "./bin/bcalm -in " + input + " -kmer-size " + to_string(k) + " -abundance-min 2  -nb-cores " + to_string(threads) + " -out " + output;
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
		char *symlinkpath = &filename[0];
		char actualpath [PATH_MAX+1];
		char *ptr;
		ptr = realpath(symlinkpath, actualpath);
		string full_path(ptr);
		file_list_graphs << full_path << endl;
		//~ cout << string(full_path) << "/" << outputDir << "/" << output + ".unitigs.fa" << endl;
	}
	return main_output + "/graphs.lst";
}


string bcalm_launcher_union(string& graph_list, uint k,  uint threads, string& main_output, string& output_bcalm)
{

	int systemRet;
	string outputDir(main_output + "/" +output_bcalm);
	systemRet=system(("mkdir -p "+ outputDir).c_str());
	string cmd("./bin/bcalm -in " + graph_list + " -kmer-size " + to_string(k) + " -abundance-min  1 -nb-cores " + to_string(threads) + " -out union_graph");
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
		//~ names = list_graphs[i];
		//~ f = &(input[0]);
		//~ char* full_path = getcwd(f, names.size());
		//~ string fp(full_path);
		//~ cout << fp << endl;
		//~ file_list_graphs << fp << "/" << names << endl;
		char *symlinkpath = &input[0];
		char actualpath [PATH_MAX+1];
		char *ptr;
		ptr = realpath(symlinkpath, actualpath);
		string rp(ptr);
		file_list_graphs << rp  << endl;

	}
	return main_output + "/graphs.lst";
}

#endif
