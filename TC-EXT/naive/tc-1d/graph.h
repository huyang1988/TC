#ifndef	GRAPH_H
#define	GRAPH_H
#include "wtime.h"
#include <string>
#include <iostream>
#include "comm.h"
#include "edge.h"

class graph{
	
	//variable
public:

	std::string	filename;
	index_t	upperEdgeCount;
	vertex_t 	vert_count;
	
	Edge		*OrientedEdge;

	Edge_chunk	*buffer;

	
	index_t		edge_count;

	vertex_t *	currentAdj;
	index_t	*	currentBegin;

//	vertex_t**	partAdj;
//	index_t**	partBegin;
	index_t*	partEdgeCount;
	
	index_t		*count;	
	int 		*valid;


public:
	graph() {};
	graph(	std::string filename);//,

	void part_validation();
	

};
#endif
