#include "graph.h"
#include "wtime.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <pthread.h>
#define N 256*256
using namespace std;



int main(int args, char *argv[]) {
//	pthread_t thd1;
//	pthread_t *thd = new pthread_t[GPU_NUM];
	std::cout<<"Input format: ./exe graph-file-name"
						<<" (json formated file)\n";

	if(args != 2) return -1;
	string json_file 	= argv[1];
	graph *graph_d 		= new graph	(json_file); 
double t0=wtime();
	graph_d->validation();
double t1=wtime();
	cout<<"total time = "<<t1-t0<<" secondes"<<endl;

	return 0;
}
