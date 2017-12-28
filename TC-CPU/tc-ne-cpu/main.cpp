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
	//pthread_t *thd = new pthread_t[GPU_NUM];
	std::cout<<"Input format: ./exe graph-file-name"
						<<" (json formated file)\n";

	if(args != 2) return -1;
	string json_file 	= argv[1];
	graph *graph_d 
		= new graph	(json_file); 
	cout<<"no edge reduction\n";
//	mygraph=graph_d;
// test sort part!
//	graph_d->sort();
//	cout<<"naive reduce\n";
//	graph_d->reduce();
//	graph_d->reverse_rank_by_degree();
//	cout<<"rank by degree\n";
//	graph_d->rank_by_degree();
//	cout<<"heap\n";
//	graph_d->reheap();
	cout<<"calculate\n";

	double time = 0;
for(int i=0;i<1;i++){
double t0=wtime();
//	graph_d->validation();
	graph_d->bs();
//	graph_d->tc();
double t1=wtime();
	time+=t1-t0;
	cout<<"total time = "<<t1-t0<<" secondes"<<endl;
}
//	time = time/3;
//	cout<<"average time of 3 round = "<<time<<" secondes"<<endl;
/*
	graph_d->triangle_count();
	graph_d->validation();

	int err=0;
	int count1=graph_d->count[0];
	int count2=0;
	for(int i=0; i<graph_d->vert_count; i++){
		count2+= graph_d->valid[i];
	}
	err = count1-count2;
	printf("count1 = %d, count2 = %d\n",count1,count2);
	printf("err number = %d\n",err);
*/
// test scan part!'

	return 0;
}
