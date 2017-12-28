#include "graph.h"
//#include "graph.c"
#include "worker.h"
//#include "worker.c"
#include "wtime.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include <omp.h>

#include <stdint.h>
#include <sys/stat.h>

#define N 256*256
using namespace std;



int main(int argc, char **argv) {

//	int myid, numprocs;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

//	cout<<"total machine # is "<<numprocs<<endl;
//	cout<<"my rank is "<<myid<<endl;
	
	

	std::cout<<"Input format: ./exe graph-file-name (json formated file)\n";
	if(argc != 2) return -1;
	string json_file 	= argv[1];

	
	
	
	
//	graph *graph_d 
//		= new graph	(json_file); 
	/*
	CSR *csr = new CSR(json_file);
	csr->CSR_load(0,0);
	csr->CSR_free();
	csr->CSR_load(0,1);
	csr->CSR_free();
	csr->CSR_load(1,0);
	csr->CSR_free();
	csr->CSR_load(1,1);
	csr->CSR_free();

	*/
	/*
	Edge_chunk * ec = new Edge_chunk(json_file);
	ec-> Edge_open(0,0);
	for(int i = 0; i<ec->chunk_num; i++){
//	for(int i = 0; i<5; i++){
		cout<<"chunk "<<i<<endl;
		ec->Edge_reload();
		for(int j =0; j<ec->chunk_size; j++){
			cout<<ec->edge[j].A<<" - "<<ec->edge[j].B<<"   ";
		}
		cout<<endl;
	}
	ec-> Edge_free();
*/


double t0=wtime();
	worker *local = new worker(json_file);
	local -> work();
double t1=wtime();

cout<<"time = "<<t1-t0<<endl;	
cout<<"time of io = "<<local->time_io<<endl;	
cout<<"time of tc = "<<local->time_compute<<endl;	



//	MPI_Finalize();
	return 0;
}

