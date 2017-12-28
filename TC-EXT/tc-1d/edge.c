#include "edge.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include "comm.h"
//#include "wtime.h"
#include <fstream>
#include <omp.h>

#define FILE_NOT_EXIST	1
#define FILE_EXIST	0
using namespace std;

Edge_chunk::Edge_chunk(string Dir_name){
	dir_name = Dir_name;
	chunk_id = -1;
	edge_size = 0;
	chunk_num = 0;
	edge = (Edge *)malloc(BufferSize*sizeof(Edge));

}

//void Edge_chunk::Edge_open(int row, int col){
void Edge_chunk::Edge_open(){

	string s_edge = dir_name + "/edge";
	char* edge_file = const_cast<char*>(s_edge.c_str());


	edge_size = fsize(edge_file)/sizeof(Edge);
	chunk_num = (edge_size-1)/BufferSize + 1;
	
	pFile = fopen(edge_file,"rb");

}

void Edge_chunk::Edge_reload(){

	chunk_id++;
	if(chunk_id == (chunk_num-1) ){
		chunk_size = edge_size - (chunk_num-1)*BufferSize;
	}
	else{
		chunk_size = BufferSize;
	}
	fread(edge,sizeof(Edge),chunk_size,pFile);

//	for(int j =0; j<chunk_size; j++){
//		cout<<edge[j].A<<"  -  "<<edge[j].B<<endl;
//	}
}



void Edge_chunk::Edge_free(){
//cout<<"start free!"<<endl;
	chunk_id = -1;
//	chunk_size = 0;
	fclose(pFile);	
	pFile = NULL;
}








