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
#include "graphLight.h"



using namespace std;



int main(int argc, char ** argv){
	graphLight GL("lol.fa", 12, 5);
	cout<<GL.getUnitig(0)<<endl;
	cout<<GL.getUnitig(1)<<endl;
	cout<<GL.getUnitig(2)<<endl;
	cout<<GL.getUnitig(3)<<endl;
	return 0;
}

