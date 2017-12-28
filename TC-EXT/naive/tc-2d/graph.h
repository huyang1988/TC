#ifndef	GRAPH_H
#define	GRAPH_H

#include "comm.h"



class CSR{
public:
	std::string dir_name;
	int col_id;
	int row_id;
	vertex_t	beg_size;
	index_t		adj_size;
	index_t		*begin;
	vertex_t	*adj;

public:
	CSR(std::string Dir_name);
	void CSR_load(int row,int col);
	void CSR_free();


};



class Edge_chunk{
public:
	std::string 	dir_name;
	int		col_id;
	int		row_id;
	int		chunk_id;
	index_t	chunk_size;//the size of next chunk to reload
	int		chunk_num;

	index_t		edge_size;//the size of this edge list partition

	Edge		*edge;
	FILE		*pFile;

public:
	Edge_chunk(std::string Dir_name);
	void Edge_open(int row,int col);
	void Edge_reload();
	void Edge_free();

};



/*
class graph{
	
	//variable
public:
	vertex_t 	vert_count;
	vertex_t	*adj_list;
	vertex_t	*head_list;
	index_t	*beg_pos;
	
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

	//gpu data
	GPU_data *gdata;

	//constructor
public:
	graph() {};
	graph(	std::string filename);//,
	
	void preproc();
	void reduceResult();

//	void initDevice(int GPU_id,int Part_id);
//	void DeviceCompute(int GPU_id,index_t Chunk_id);
//	void gpuReduce(int GPU_id);
//	void gpuProc(int GPU_id);

	void cpuProc(int P);
	void cpuCompute(int Part_id, index_t Chunk_id);

};
*/

#endif
