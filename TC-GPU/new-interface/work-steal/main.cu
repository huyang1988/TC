//Graph format: Json based format
//Storage format: 
//struct{
//		int: src_ver
//		Arr: [ver_0|ver_1|ver_2|...]
//		Int: num_conn_ver
//	}
/* main.cu */
#include "graph.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <pthread.h>
#define N 256*256
using namespace std;



int main(int args, char *argv[]) {
//	pthread_t thd1;
	std::cout<<"Input format: ./exe graph-file-name"
						<<" (json formated file)\n";

	if(args != 2) return -1;
	string json_file 	= argv[1];
	graph *graph_d 
		= new graph	(json_file); 
	
	cout<<"GPU  NUMBER = "<<GPU_NUM<<endl;
	cout<<"PART NUMBER = "<<PART_NUM<<endl;

	mygraph=graph_d;
	cout<<"rank by degree\n";

	double total_t=0;
int r=2;
for(int n=0; n<r; n++){
	//initiation, clean the footprint of last executionfor(int i=0; i<GPU_NUM+1; i++){
	for(int i=0; i<DEV_NUM; i++){
		graph_d->ds_complete[i]=0;
		graph_d->ds_help[i]=0;
	}

#pragma omp parallel for
	for(int i=0; i<PART_NUM * graph_d->ChunkNum; i++){
		graph_d->ds_status[i]=0;
	}

	double t0=wtime();
//	index_t total = 0;

//int tid;
omp_set_nested(1);
{
#pragma omp parallel for num_threads(DEV_NUM) schedule(static)
	for(int i=0; i<DEV_NUM; i++){
		if(i<GPU_NUM){
			graph_d->gpuProc(i);
		}
	}
}

	
//	for(int i=0; i<GPU_NUM+1; i++){
//		total+= graph_d->count[i];
//	}
	graph_d->reduceResult();
	double t1=wtime();
	cout<<"total count "<<graph_d->count[0]<<"\n";
	cout<<"total time  "<<t1-t0<<" seconds\n";
	total_t += t1-t0;
}
cout<<"merge average time of 5 round = "<<total_t/r<<endl;
	return 0;
}
