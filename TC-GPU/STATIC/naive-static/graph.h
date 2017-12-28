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

	//after partition
	vertex_t**	partAdj;
	vertex_t**	partHead;
	index_t**	partBegin;
	index_t*	partEdgeCount;
	
	index_t		*count;	
	int 		*valid;
//dynaic scheduling
	index_t		ChunkNum;	// = roof(upperEdgeCount/BufferSize), the number of chunks for workload edge list

	index_t		**ds_count;	//ds_count[P][i] is the count of partition P chunk i	
	int		**ds_status;	//ds_status[P][i] = 0: no worker;  1: working or complete
	int		*ds_progress;	//ds_progress[p]: the ID of last chunk need worker

	double		copy_time;

	//gpu data
	GPU_data *gdata;

	//constructor
public:
	graph() {};
	graph(	std::string filename);//,

//	void validation();
//	void sort();
//	void reduce();
//	void rank_by_degree();
//	void reverse_rank_by_degree();
//	void partition();
//	void* scan(void*data);
	
	void preproc();
	void reduceResult();

	void initDevice(int GPU_id,int Part_id);
	void DeviceCompute(int GPU_id,index_t Chunk_id);
	void gpuReduce(int GPU_id);
	void gpuProc(int GPU_id);

	void cpuProc();
	void cpuCompute(int Part_id, index_t Chunk_id);

};

#include "graph.c"
//#include "kernel.cu"
#include "sort.cu"
#include "scan.cu"
#endif
