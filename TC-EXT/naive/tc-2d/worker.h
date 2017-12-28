#ifndef	WORKER_H
#define	WORKER_H

#include "comm.h"
#include "wtime.h"
#include "graph.h"

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
	//local proportion ??
	
public:
	worker(std::string s);
	void work();
	void proc(int i, int j, int k);// edge list(i,j), CSR(i,k) and CSR(i,k)
};



#endif

