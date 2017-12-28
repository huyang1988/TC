#ifndef	EDGE_H
#define	EDGE_H

#include "comm.h"



class Edge_chunk{
public:
	std::string 	dir_name;
//	int		col_id;
//	int		row_id;
	int		chunk_id;
	index_t	chunk_size;//the size of next chunk to reload
	int		chunk_num;

	index_t		edge_size;//the size of this edge list partition

	Edge		*edge;
	FILE		*pFile;

public:
	Edge_chunk(std::string Dir_name);
	void Edge_open();
	void Edge_reload();
	void Edge_free();

};




#endif
