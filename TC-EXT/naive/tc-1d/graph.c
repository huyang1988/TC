#include "graph.h"
#include "comm.h"
#include <fstream>
#include <omp.h>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include "graph.h"
using namespace std;

graph::graph(
	string jsonfile)//,
{
	filename = jsonfile;
	cout<<"read from folder "<<jsonfile<<endl;
//read metadata
	string s_meta = filename + "/metadata";
	char* meta_file = const_cast<char*>(s_meta.c_str());
	size_t meta_size = fsize(meta_file)/sizeof(vertex_t);
	FILE *mFile= fopen(meta_file,"rb");
	vertex_t *meta = (vertex_t *)malloc(fsize(meta_file));
	fread(meta,sizeof(vertex_t),meta_size,mFile);
	fclose(mFile);
	PART_NUM = meta[0];





	partEdgeCount 	= new index_t [PART_NUM];

	for(int i = 0; i<PART_NUM; i++){

		stringstream ss;
		ss<<i;
		string str = ss.str();

		string s_begin = jsonfile+"/begin-"+str;
		string s_adj = jsonfile+"/adjacent-"+str;

		char* begin_file = const_cast<char*>(s_begin.c_str());
		char* adj_file = const_cast<char*>(s_adj.c_str());

		vert_count 	 = fsize(begin_file)/sizeof(index_t)-1;	
cout<<"vert count = "<<vert_count<<endl;	
		partEdgeCount[i] = fsize(adj_file)/sizeof(vertex_t);
cout<<"part "<<i<<" edge count = "<<partEdgeCount[i]<<endl;
	}
	
	
	
	string s_edge = jsonfile+"/edge";
	char* edge_file = const_cast<char*>(s_edge.c_str());
	
		
	upperEdgeCount = fsize(edge_file)/sizeof(Edge);
cout<<"upperEdgeCount = "<<upperEdgeCount<<endl;
	OrientedEdge = new Edge[upperEdgeCount];
	FILE *pFile1 = fopen(edge_file,"rb");
	fread(OrientedEdge, sizeof(Edge), upperEdgeCount, pFile1);

	
	buffer = new Edge_chunk(filename);
	currentBegin = new index_t [vert_count+1];

//	string s_meta = filename + "/metadata";
//	size_t meta_size = fsize(meta_file)/sizeof(vertex_t);
//	FILE *mFile = fopen(meta_file,"rb");
//	vertex_t * meta = (vertex_t *)malloc (fsize(meta_file));
//	fread(meta,sizeof(vertex_t),meta_size,mFile);
//	fclose(mFile);
//
	

}




void graph::part_validation(){
//cout<<"validation "<<upperEdgeCount<<endl;
//
	long int global_count=0;
	for(int p=0; p<PART_NUM; p++){
		
		stringstream ss;
		ss<<p;
		string str = ss.str();

		string s_begin = filename+"/begin-"+str;
		string s_adj = filename+"/adjacent-"+str;

		char* begin_file = const_cast<char*>(s_begin.c_str());
		char* adj_file = const_cast<char*>(s_adj.c_str());

		
		FILE *pFile = fopen(adj_file,"rb");
		currentAdj = (vertex_t*)malloc(partEdgeCount[p]*sizeof(vertex_t));
		fread(currentAdj, sizeof(vertex_t), partEdgeCount[p], pFile);
		fclose(pFile);

		FILE *pFile3 = fopen(begin_file,"rb");
		fread(currentBegin, sizeof(index_t), vert_count+1 , pFile3);
		fclose(pFile3);
	

		buffer->Edge_open();

		for(int c=0; c<buffer->chunk_num;c++){
			index_t mycount=0;

			buffer->Edge_reload();

#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(dynamic,1024)
			for(index_t i=0; i<buffer->chunk_size; i++){
				//get the edge A-B from global graph
//				vertex_t A=OrientedEdge[i].A;
//				vertex_t B=OrientedEdge[i].B;
				vertex_t A = buffer->edge[i].A;
				vertex_t B = buffer->edge[i].B;


				//degree of current partition
				index_t m=currentBegin[A+1]-currentBegin[A];
				index_t n=currentBegin[B+1]-currentBegin[B];

//cout<<"edge: "<<A<<"-"<<B<<" "<<endl;
		//cout<<"degree: "<<m<<" "<<n<<endl;

				vertex_t *a = &currentAdj[currentBegin[A]];
				vertex_t *b = &currentAdj[currentBegin[B]];

				vertex_t u1=0;
				vertex_t v1=0;
				while(u1<m && v1<n){
					vertex_t x=a[u1];
					vertex_t y=b[v1];
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
//cout<<"find triangle "<<A<<"-"<<B<<"-"<<x<<endl;
					}
				}

			}
			global_count += mycount;

		
		}
		free(currentAdj);
		buffer->Edge_free();

	}
	cout<<"merge version tc = "<<global_count<<endl;
}
