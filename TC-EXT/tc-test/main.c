#include "graph.h"
#include "worker.h"
#include "wtime.h"

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

#define WORKTAG 1
#define DIETAG 2
using namespace std;

void master(worker* g);
void slave(worker* g);

int main(int argc, char **argv) {
	

	
	

	std::cout<<"Input format: ./exe graph-file-name (json formated file)\n";
	if(argc != 2) return -1;
	string json_file 	= argv[1];

	cout<<json_file<<endl;


double t0=wtime();
	worker *local = new worker(json_file);
	local -> work();
double t1=wtime();

cout<<"time = "<<t1-t0<<endl;	
cout<<"time of io = "<<local->time_io<<endl;	
cout<<"time of tc = "<<local->time_compute<<endl;	
	
	return 0;
}

void master(worker* g)
{
	index_t global_count=0;
	cout<<"master"<<endl;
	int ntasks;
//	double       result;
	index_t	result;
	MPI_Status     status;
	MPI_Comm_size(MPI_COMM_WORLD,   &ntasks);// get the total process number
// Seed the slaves.
	index_t max_work = g->subproblem;
cout<<"total work number = "<<max_work<<endl;	
	index_t rem_work = max_work;
	index_t rem_result = max_work;	
	index_t work = 0;
	
	double t0=wtime();
	for (int rank = 1; rank < ntasks; rank++) {
		MPI_Send(&work, //         message buffer 
			1,              // one data item 
			MPI_LONG,        // data item is an integer 
			rank,           // destination process rank 
			WORKTAG,        // user chosen message tag 
			MPI_COMM_WORLD
		);// always use this
//cout<<"send work to proc "<<rank<<endl;	
		work++;
		rem_work--;
		if(rem_work ==0){break;}
	}

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
		global_count+=result;
//	cout<<"rcv work to proc "<<status.MPI_SOURCE<<endl;
	}
	
	double t1=wtime();
	for (int rank = 1; rank < ntasks; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
	}

	cout<<"total count from all machine = "<<global_count<<"\n";
	cout<<"total time of all machine = "<<t1-t0<<"\n";
//collect result
	

}




void slave(worker* g)
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
//	int GPU_id = 0;//set GPU for this process
	
	cout<<"slave"<<endl;
	
	index_t	result = 0;
	int	work;

	MPI_Status	status;
double t0=wtime();
	for (;;) {
		MPI_Recv(&work, 1, MPI_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (status.MPI_TAG == DIETAG) {
double t1=wtime();
cout<<"proc"<<myrank<<" time = "<<t1-t0<<"\n";
			return;
		}

//work region	
cout<<"slave recv "<<work<<endl;

		g->proc(work);	
		result = g->count;
	
//work region
//		result = // do the work ;
		MPI_Send(&result, 1, MPI_LONG, 0, 0, MPI_COMM_WORLD);
	}
}

