//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...

#include "graph.h"
//#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include<fcntl.h>
#include <stdio.h>
#include "comm.h"
#include <omp.h>

using namespace std;


CSR::CSR(string Dir_name){
	dir_name = Dir_name;
	row_id = -1;
	col_id = -1;

	beg_size = 0;
	adj_size = 0;

}

void CSR::CSR_load(int row, int col){
	stringstream ss;
	ss << row;
	ss << "-";
	ss << col;
	string postfix = ss.str();
	

	string s_begin = dir_name + "/begin" + postfix;
	string s_adj = dir_name + "/adjacent" + postfix;
	
	char* begin_file = const_cast<char*>(s_begin.c_str());
	char* adj_file = const_cast<char*>(s_adj.c_str());

//cout<<begin_file<<endl;
//cout<<adj_file<<endl;

	beg_size = fsize(begin_file)/sizeof(index_t);
	adj_size = fsize(adj_file)/sizeof(vertex_t);
	
//cout<<"vert:"<< beg_size-1<<"  edge: "<< adj_size<<endl;

//allocation

	int fd;	
//	FILE *pFile= fopen(adj_file,"rb");
	adj = (vertex_t *)malloc(fsize(adj_file));
//	fread(adj,sizeof(vertex_t),adj_size,pFile);
//	fclose(pFile);
	fd = open(adj_file , O_RDONLY);
	int rt = read(fd, adj, fsize(adj_file));
	//cout<<"fd = "<<fd<<" : "<<rt ;
       	perror("Error:");
	cout<<endl;
	close(fd);
	
	FILE *pFile3 = fopen(begin_file,"rb");
	begin = (index_t *)malloc(fsize(begin_file));
	fread(begin,sizeof(index_t),beg_size,pFile3);
	fclose(pFile3);
		

	index_t verify_adj_size = begin[beg_size-1];
//	cout<<"adj list size = "<<adj_size<<" while last of begin array = "<<verify_adj_size<<endl;
	if(verify_adj_size!=adj_size){
		cout<<"CSR partition "<<row<<"-"<<col<<" edge count does not match" <<endl;
		cout<<"adj list size = "<<adj_size<<" while last of begin array = "<<verify_adj_size<<endl;
		cout<<"beg list size = "<<beg_size<<endl;
	}

//test
	/*
cout<<"CSR "<<row<<"-"<<col<<endl;
	for(vertex_t i = 0; i<beg_size; i++){
		cout<<begin[i]<<" ";
	}
	cout<<endl;
	for(index_t i = 0; i<adj_size; i++){
		cout<<adj[i]<<" ";
	}
	cout<<endl;
	*/
}

void CSR::CSR_free(){
	row_id = -1;
	col_id = -1;
	beg_size = 0;
	adj_size = 0;
	free(begin);
	free(adj);
}





Edge_chunk::Edge_chunk(string Dir_name){
	dir_name = Dir_name;
	row_id = -1;
	col_id = -1;
	chunk_id = -1;
	edge_size = 0;
	chunk_num = 0;
	edge = (Edge *)malloc(BufferSize*sizeof(Edge));

}

void Edge_chunk::Edge_open(int row, int col, int k){
	stringstream ss;
	ss << row;
	ss << "-";
	ss << col;
	ss<<"-";
	ss<<k;
	string postfix = ss.str();
	string s_edge = dir_name + "/edge" + postfix;
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
	row_id = -1;
	col_id = -1;
	chunk_id = -1;
//	chunk_size = 0;
	fclose(pFile);	
	pFile = NULL;
}








