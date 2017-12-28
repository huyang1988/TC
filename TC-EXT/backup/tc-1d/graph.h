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
#include "wtime.h"
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
	
	Edge		*OrientedEdge;
	index_t	upperEdgeCount;
	
	index_t		edge_count;

	vertex_t**	partAdj;
	index_t**	partBegin;
	index_t*	partEdgeCount;
	
	index_t		*count;	
	int 		*valid;


public:
	graph() {};
	graph(	std::string filename);//,

	void part_validation();
	

};

//#include "kernel.cu"
//#include "scan.cu"
#endif
