//graph.h
//Graph format: Json based format: [src_id, src_weigh,[[connected_ver_0, edge_weight],[connected_ver_1, edge_weight],[connected_ver_2, edge_weight]]]
//Storage format: 
//struct{
//		int: src_ver
//		Arr: [ver_0|ver_1|ver_2|...]
//		Int: num_conn_ver
//	}
#ifndef	GRAPH_H
#define	GRAPH_H

#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include "comm.h"

class graph{
	
	//variable
public:
	vertex_t 	vert_count;
	vertex_t	*adj_list;
	vertex_t	*head_list;
//	index_t	*adj_card;
	index_t	*beg_pos;
	//after sort
	vertex_t	*upperAdj;
	vertex_t	*upperHead;
	index_t	*upperBegin;
	index_t	*upperDegree;
	index_t	upperEdgeCount;
	
	index_t		edge_count;

	//after partition
	vertex_t**	partAdj;
	vertex_t**	partHead;
	index_t**	partBegin;
	index_t**	partDegree;
	index_t*	partEdgeCount;
	
	index_t		*count;	
	int 		*valid;


	//gpu data
//	GPU_data *gdata;

	//constructor
public:
	graph() {};
	graph(	std::string filename);//,

	void info();
	void triangle_count();
	void sort();
	void partition();
	void* scan(void*data);

	void process_graph();//produce data to run CMU code

	void validation();
	void bsvalidation();

};

#include "graph.c"
//#include "kernel.cu"
#include "sort.c"
#include "scan.cu"
#endif
