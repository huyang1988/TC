#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <stdint.h>
#include <sys/stat.h>
inline off_t fsize(const char *filename) {
	struct stat st; 
	if (stat(filename, &st) == 0)
		return st.st_size;
	return -1; 
}


int main(int argc, char **argv)
{
	int myid, numprocs;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//find my filename
	char filename[256];
	sprintf(filename,"hang_graph_%d_of_%d.bin",myid, numprocs);

	int nmemb=fsize(filename)/sizeof(int);
	int *data=new int[nmemb];
	
	////////////////////////////
	//for(int i=0;i<numprocs;i++)
	//	data[i]=rand()+myid;
	//FILE *file=fopen(filename,"wb");
	//fwrite(data,sizeof(int), numprocs,file);
	/////////////////////////////////////////////////////
	

	//open and read data from my file
	FILE *file=fopen(filename, "rb");
	if(file==NULL)
	{
		std::cout<<"wrong open @ "<<myid<<"\n"; 
		exit(-1);
	}
	fread(data,sizeof(int),nmemb,file);
	fclose(file);

	//check the correctness
	for(int i=0;i<numprocs;i++)
	{
		if(i==myid)
		{
			std::cout<<"Process "<<i<<": ";
			for(int j=0;j<numprocs;j++)
				std::cout<<data[j]<<" ";
			std::cout<<"\n";
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	//finish
	MPI_Finalize();
	return 0;
}
