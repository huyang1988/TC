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
//	vertex_t	*upperAdj;
//	vertex_t	*upperHead;
	
	Edge		*OrientedEdge;

	index_t	*upperBegin;
	index_t	*upperDegree;
	index_t	upperEdgeCount;
	
	index_t		edge_count;

	//1-d vertical partition
	vertex_t	**partAdj;
	vertex_t	**partHead;
	index_t		**partBegin;
	index_t		*partEdgeCount;
	
	index_t		*count;	

	//2-d partitioned data structures
	//edge
	Edge		**p2Edge;//Edge[k] use k=i*N+j represnet E[i][j]
	index_t		*p2EdgeCount;

	//csr
	index_t		**p2CSR;//number of edges of CSR partition
	
	vertex_t	***p2Adj;//p2Adj[i][j] is a pointer to partAdj
	index_t		***p2Begin; //p2Begin[i][j] is an array of begin position 

	vertex_t	*metadata;
	



	//constructor
public:
	graph() {};
	graph(	std::string filename);//,
	void part_validation();
	void preproc();		//CSR 2d partitioner

	void vertical_partition();	//output: 1-D vertical partition partAdj, partBegin
	void further_partition();	//input: vertical partitions, memory size, output: partition number P2, metadata, 2D CSR
	void edge_2d();		//edge 2d partitioner according to metadata
	void write_back();	//write back the partitioned bit files to disk


};

#include "graph.c"
#include "sort.c"
#endif
