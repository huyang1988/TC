//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...


//#include "graph.h"
#include "comm.h"
#include <fstream>
#include <omp.h>

#define FILE_NOT_EXIST	1
#define FILE_EXIST	0
using namespace std;

graph::graph(
	string jsonfile)//,
{
	cout<<"read from folder "<<jsonfile<<endl;
	
	string s_begin = jsonfile+"/begin.bin";
	string s_adj = jsonfile+"/adjacent.bin";
	string s_head = jsonfile+"/head.bin";
	string s_degree = jsonfile+"/degree.bin";

	char* begin_file = const_cast<char*>(s_begin.c_str());
	char* adj_file = const_cast<char*>(s_adj.c_str());
	char* head_file = const_cast<char*>(s_head.c_str());
	char* degree_file = const_cast<char*>(s_degree.c_str());

	vert_count = fsize(begin_file)/sizeof(index_t) - 1;
	edge_count = fsize(head_file)/sizeof(vertex_t);

	cout<<"vert:"<< vert_count<<"  edge: "<<edge_count<<endl;


	FILE *pFile= fopen(adj_file,"rb");
	adj_list = (vertex_t *)malloc(fsize(adj_file));
	fread(adj_list,sizeof(vertex_t),edge_count,pFile);
	fclose(pFile);

	FILE *pFile1= fopen(head_file,"rb");
	head_list = (vertex_t *)malloc(fsize(head_file));
	fread(head_list,sizeof(vertex_t),edge_count,pFile1);
	fclose(pFile1);

//	FILE *pFile2= fopen(degree_file,"rb");
//	adj_card = (index_t *)malloc(fsize(degree_file));
//	fread(adj_card,sizeof(index_t),vert_count,pFile2);
//	fclose(pFile2);

	FILE *pFile3 = fopen(begin_file,"rb");
	beg_pos = (index_t *)malloc(fsize(begin_file));
	fread(beg_pos,sizeof(index_t),vert_count+1,pFile3);
	fclose(pFile3);
	count = (index_t *)malloc(256*256*sizeof(index_t));
//	valid = (int *)malloc(vert_count*sizeof(int));
	ds_count  =   (long*)malloc(GPU_NUM * GPU_NUM *sizeof(long));	
	ds_status =   (int*)malloc(GPU_NUM *GPU_NUM * sizeof(int));
	ds_complete = (int*)malloc(GPU_NUM * sizeof(int));
	ds_help   =   (int*)malloc(GPU_NUM * sizeof(int));

	for(int i=0; i<GPU_NUM*GPU_NUM; i++){
		ds_status[i]=0;
		ds_count[i]=0;
	}
	for(int i=0; i<GPU_NUM; i++){
		ds_complete[i]=0;
		ds_help[i]=0;
	}
}

/*
void graph::validation(){
	int mycount=0;
	for(int i=0; i<upperEdgeCount; i++){
		int U=upperHead[i];
		int V=upperAdj[i];
		int m=upperDegree[U];
		int n=upperDegree[V];

		int *u = &upperAdj[upperBegin[U]];
		int *v = &upperAdj[upperBegin[V]];

		int u1=0;
		int v1=0;
		while(u1<m && v1<n){
			int x=u[u1];
			int y=v[v1];
			if(x<y){
				u1++;
			}
			else if(x>y){
				v1++;
			}
			else if(x==y){
				u1++;
				v1++;
				mycount++;
			}
		}
	}
	printf("validation version tc = %d\n",mycount);
}
*/

void graph::triangulation(){
	index_t mycount=0;
#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(dynamic,512)
	for(index_t i=0; i<edge_count; i++){
		int U=head_list[i];
		int V=adj_list[i];
		int m=beg_pos[U+1]-beg_pos[U];
		int n=beg_pos[V+1]-beg_pos[V];

		int *u = &adj_list[beg_pos[U]];
		int *v = &adj_list[beg_pos[V]];

		int u1=0;
		int v1=0;
		while(u1<m && v1<n){
			int x=u[u1];
			int y=v[v1];
			if(x<y){
				u1++;
			}
			else if(x>y){
				v1++;
			}
			else if(x==y){
				u1++;
				v1++;
				mycount++;
				break;
			}
		}
	}
	printf("count of tc_edge = %d\n",mycount);

}
