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

//#define N 256*256
#define WORKTAG 1
#define DIETAG 2

using namespace std;

void master(graph* g);
void slave(graph* g);



int main(int argc, char **argv) {

	int myrank, nprocs;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	cout<<"total machine # is "<<nprocs<<endl;
	cout<<"my rank is "<<myrank<<endl;
	
	

//	pthread_t thd1;
	std::cout<<"Input format: ./exe graph-file-name"
						<<" (json formated file)\n";

	if(argc != 2) return -1;
	string json_file 	= argv[1];
	graph *graph_d 
		= new graph	(json_file); 
	
	
	index_t proc_count;
	double total_t=0;

	MPI_Barrier(MPI_COMM_WORLD);

	if (myrank == 0) {
		master(graph_d);
	} else {
		slave(graph_d);
	}



	MPI_Finalize();
	return 0;
}

void master(graph* g)
{
	index_t proc_count;
	index_t global_count=0;
	cout<<"master"<<endl;
	int ntasks;
//	double       result;
	index_t	result;
	MPI_Status     status;
	MPI_Comm_size(MPI_COMM_WORLD,   &ntasks);// get the total process number
// Seed the slaves.
	index_t max_work = g->ChunkNum * PART_NUM;
cout<<"total work number = "<<max_work<<endl;	
	index_t rem_work = max_work;
	index_t rem_result = max_work;	
	index_t work = 0;
	
	double t0=wtime();
	for (int rank = 1; rank < ntasks; rank++) {
//		work = ;		// get_next_work_request ;
		MPI_Send(&work, //         message buffer 
			1,              // one data item 
			MPI_LONG,        // data item is an integer 
			rank,           // destination process rank 
			WORKTAG,        // user chosen message tag 
			MPI_COMM_WORLD
		);// always use this
//	cout<<"send work to proc "<<rank<<endl;	
		work++;
		rem_work--;
		if(rem_work ==0){break;}
	}

// * * Receive a result from any slave and dispatch a new work
// * * request work requests have been exhausted.
// 	work = // get_next_work_request ;
	while (rem_work>0){//( valid new work request ) {
		MPI_Recv(&result,	//        message buffer 
			1,              	// one data item 
			MPI_LONG,     	// of type double real 
			MPI_ANY_SOURCE, 	// receive from any sender 
			MPI_ANY_TAG,    	// any type of message 
			MPI_COMM_WORLD, 	// always use this 
			&status
		);       	// received message in

		global_count+=result;

//cout<<"recv work from proc "<<status.MPI_SOURCE<<endl;
		MPI_Send(&work, 1, MPI_LONG, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);
//cout<<"send work to proc "<<status.MPI_SOURCE<<endl;
		work ++; 			// get_next_work_request ;
		rem_work--;
		rem_result--;
	}

// * * Receive results for outstanding work requests.
///*
	for (int i = 1; i < ntasks; i++) {
		MPI_Recv(
				&result,    //result value
				1, 
				MPI_LONG, //result type
				MPI_ANY_SOURCE,
				MPI_ANY_TAG, 
				MPI_COMM_WORLD, 
				&status
		);
//	cout<<"rcv work to proc "<<status.MPI_SOURCE<<endl;
	}
	
	double t1=wtime();
//*/
// * * Tell all the slaves to exit.
	for (int rank = 1; rank < ntasks; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}

	cout<<"total count from all machine = "<<global_count<<"\n";
	cout<<"total time of all machine = "<<t1-t0<<"\n";
//collect result
	

}



//--------------------------------------


void slave(graph* g)
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	int GPU_id = 0;//set GPU for this process
	
	index_t proc_count = 0;
	index_t global_count=0;
	cout<<"slave"<<endl;
	
	int	P_curr = -1;
	index_t	Chunk_id;
	index_t ChunkNum = g->ChunkNum;
//	double              result;
	index_t	result = 0;
	int	work;

	MPI_Status	status;
double t0=wtime();
	for (;;) {
		MPI_Recv(&work, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
// * * Check the tag of the received message.
		if (status.MPI_TAG == DIETAG) {
			if(P_curr!=-1){
				gpuReduce(g,GPU_id);
			}
double t1=wtime();
cout<<"proc"<<myrank<<" time = "<<t1-t0<<"\n";
			return;
		}

//work region	
//cout<<"slave recv "<<work<<endl;
		int P = (int)work/ChunkNum;
		Chunk_id = work%ChunkNum;

		if(P_curr != P){//load part P into GPU memory
			if(P_curr!=-1){
				gpuReduce(g,GPU_id);
			}
				
			initDevice(g,GPU_id,P);
			P_curr = P;
		}

		DeviceCompute(g,GPU_id,Chunk_id);
		result = g->count[0];
	
//work region
//		result = // do the work ;
		MPI_Send(&result, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
	}
}
