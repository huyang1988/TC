#include "graph.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <pthread.h>
#include <cmath>
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
	graph_d->further_partition();
//cout<<"further partition"<<endl;
	graph_d->edge_2d();
cout<<"edge partition"<<endl;

//	graph_d->write_back();

	return 0;
}
