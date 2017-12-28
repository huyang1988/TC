//Graph format: Json based format
//Storage format: 
//struct{
//		int: src_ver
//		Arr: [ver_0|ver_1|ver_2|...]
//		Int: num_conn_ver
//	}
/* main.cu */
#include "graph.h"
#include "graph.c"
#include "wtime.h"
#include "sort.cu"
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
	graph *graph_d 
		= new graph	(json_file); 
//	mygraph=graph_d;
// test sort part!
//	graph_d->sort();
//	graph_d->reduce();
//	graph_d->reverse_rank_by_degree();
double t0=wtime();
//	graph_d->rank_by_degree();
	graph_d->preproc();
double t1=wtime();
cout<<"pre-processing time = "<<t1-t0<<" secondes"<<endl;

t0=wtime();
	graph_d->rank_by_degree();
//	graph_d->preproc();
t1=wtime();
cout<<"orientation time = "<<t1-t0<<" secondes"<<endl;
//	graph_d->partition();
	double time=0;
//for(int i=0;i<5;i++){

	return 0;
}
