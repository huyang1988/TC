//Graph format: Json based format
//Storage format: 
//struct{
//		int: src_ver
//		Arr: [ver_0|ver_1|ver_2|...]
//		Int: num_conn_ver
//	}
/* main.cu */
#include "graph.h"
#include "wtime.h"
#include "scan.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#include <stdint.h>
#include <sys/stat.h>

#define N 256*256
using namespace std;



int main(int argc, char **argv) {

	int myid, numprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	cout<<"total machine # is "<<numprocs<<endl;
	cout<<"my rank is "<<myid<<endl;
	
	

//	pthread_t thd1;
	std::cout<<"Input format: ./exe graph-file-name"
						<<" (json formated file)\n";

	if(argc != 2) return -1;
	string json_file 	= argv[1];
	graph *graph_d 
		= new graph	(json_file); 
	
	graph_d->preproc();
	
	index_t proc_count;
	double total_t=0;

int r=3;
for(int n=0; n<r; n++){
	MPI_Barrier(MPI_COMM_WORLD);
	double t0=wtime();
	gpuProc(graph_d,0,myid,numprocs);
/*
omp_set_nested(1);
{
#pragma omp parallel for num_threads(DEV_NUM)
	for(int i=0; i<DEV_NUM; i++){
//		tid = omp_get_thread_num();
		if(i<GPU_NUM){
			gpuProc(graph_d,i,myid,numprocs);
		}
//		else if(i == GPU_NUM){
//			graph_d->cpuProc(myid);
//		}
		else cout<<"tid = "<<i<<endl;
	}
}
*/
	
	graph_d->reduceResult();
	proc_count = graph_d->count[0];
//cout<<"rank "<<myid<<" count "<<proc_count<<endl; 
	
	double t2=wtime();
cout<<"rank "<<myid<<" time "<<t2-t0<<endl; 


	index_t global_count=0;

	MPI_Reduce(
			&proc_count,
			&global_count,
			1,
			MPI_LONG,
			MPI_SUM,
			0,
			MPI_COMM_WORLD
			);
	MPI_Barrier(MPI_COMM_WORLD);
	double t1=wtime();
	if(myid==0){
		cout<<"total count from all machine = "<<global_count<<"\n";
		cout<<"total time of all machine = "<<t1-t0<<"\n";
	}
	total_t += t1-t0;
}
if(myid==0){
	cout<<"average time of 5 round = "<<total_t/r<<endl;
}



	MPI_Finalize();
	return 0;
}
