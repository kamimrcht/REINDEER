
#include "utils.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifndef BCA
#define BCA
using namespace std;

void get_list_graphs_fof(string& fof, vector<string>& list_graphs);

string bcalm_launcher_single(string& input_fof, uint k, uint threads, string& main_output, string& output_bcalm, bool PE);

string bcalm_launcher_union(string& graph_list, uint k, uint threads, string& main_output, string& output_bcalm);
void bcalm_cleanup();

string getRealPaths(string& fof, string& main_output);
#endif
