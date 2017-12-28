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

	PART_NUM = 1 + (int)sqrt((graph_d->edge_count-1)/PartitionSize);
//	PART_NUM = (int)sqrt((graph_d->edge_count-1)/PartitionSize);
	cout<<"PART NUMBER = "<<PART_NUM<<endl;

	cout<<"rank by degree\n";
	
	graph_d->vertical_partition();
	
//	graph_d->part_validation();

	graph_d->further_partition();
	graph_d->edge_2d();

//	graph_d->write_back();

	return 0;
}
