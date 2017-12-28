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
	int PART_NUM;

	index_t		upperEdgeCount;
	vertex_t 	vert_count;
	index_t*	partEdgeCount;


	Edge_chunk	*buffer;
	vertex_t *	currentAdj;
	index_t	*	currentBegin;

	
	double 		time_compute;
	double		time_io;


public:
	graph() {};
	graph(	std::string filename);//,

	void part_validation();
	

};
#endif
