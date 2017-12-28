#ifndef	WORKER_H
#define	WORKER_H

#include "comm.h"
#include "wtime.h"
#include "graph.h"


typedef struct GPUData{


	vertex_t * adj1;
	index_t  * beg1;

	vertex_t * adj2;
	index_t  * beg2;
	
	Edge     * buffer;

	index_t  * count;
	index_t  host_count;

	long	 currentChunkSize;
	

} GPU_data;


class worker{

public:
	CSR		*csr1;
	CSR		*csr2;
	Edge_chunk	*buffer;
	
	double time_io;
	double time_compute;

	int		N1;
	int		N2;
	vertex_t	*offset;
	
	index_t		count;
	long		edge_error;
	//local proportion ??
	GPU_data gdata;
	
public:
	worker(std::string s);
	void work();
	void proc(int i, int j, int k);// edge list(i,j), CSR(i,k) and CSR(i,k)

	void init_gpu();
	void call_gpu();
	void free_gpu();
};



#endif

