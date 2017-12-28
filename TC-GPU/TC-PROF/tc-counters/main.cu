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
	pthread_t *thd = new pthread_t[GPU_NUM];
	std::cout<<"Input format: ./exe graph-file-name"
						<<" (json formated file)\n";

	if(args != 2) return -1;
	string json_file 	= argv[1];
	graph *graph_d 
		= new graph	(json_file); 
	mygraph=graph_d;
// test sort part!
//	graph_d->sort();
//	graph_d->reduce();
	graph_d->rank_by_degree();
//	graph_d->reverse_rank_by_degree();
	graph_d->partition();
//	graph_d->validation();
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

	total_count1=0;
	total_count2=0;

	double t0=wtime();
	int id[GPU_NUM];
	for(int i=0; i<GPU_NUM;i++){
//		cout<<"loop "<<i<<"\n";
		id[i] = i;
		pthread_create(&thd[i],NULL,part_scan,&id[i]);
//		part_scan(&id[i]);
	}

	for(int i=0; i<GPU_NUM;i++){
		pthread_join(thd[i],NULL);
	}

	long int total=0;
	for(int i=0;i<GPU_NUM;i++){
		total+= graph_d->count[i];
	}
	double t1=wtime();
	cout<<"total count "<<total<<"\n";
	cout<<"total time  "<<t1-t0<<" seconds\n";
	
	cout<<"total mem reads  "<<total_count1<<"\n";
	cout<<"total divergence "<<total_count2<<"\n";

	return 0;
}
