#include "cuUtil.cu"
#include "scan.h"
#include "comm.h"
#include "wtime.h"
#include "graph.h"
#include "worker.h"
#include <stdio.h>
#include "iostream"
#define max_thd 256 
#define max_block 256
using namespace std;

__global__ void warp_binary_kernel
(
	//input data
	vertex_t	*adj1,
	index_t		*beg1,
	vertex_t	*adj2,
	index_t		*beg2,
	Edge*		buffer,
	//parameters
	index_t		bufferSize,
	vertex_t	offsetA,
	vertex_t	offsetB,
	//output
	index_t*	count	
)
{
	//phase 1, partition
	index_t tid = (threadIdx.x + blockIdx.x * blockDim.x)/32;
	index_t mycount=0;
	__shared__ index_t local[max_thd];

	int i = threadIdx.x%32;
	int p = threadIdx.x/32;

	while(tid<bufferSize){
		vertex_t A = buffer[tid].A - offsetA;
		vertex_t B = buffer[tid].B - offsetB;
		index_t m = beg1[A+1]-beg1[A];//degree[A];
		index_t n = beg2[B+1]-beg2[B];//degree[B];
//if(i==0) printf("A %d B %d\n");
		vertex_t* a = &(adj1[beg1[A]]);
		vertex_t* b = &(adj2[beg2[B]]);
		
		index_t tempd;
		vertex_t *tempa;	
		if(m<n){
			tempa = a;
			a = b;
			b = tempa;
			tempd = m;
			m = n;
			n = tempd;
		}


		
	//initial cache
		local[p*32+i]=a[i*m/32];	
		__syncthreads();
			
	//search
		int j=i;
		while(j<n){
			vertex_t X = b[j];
			vertex_t Y;
			//phase 1: cache
			int bot = 0;
			int top = 32;
			int r;
			while(top>bot+1){
				r = (top+bot)/2;
				Y = local[p*32+r];
				if(X==Y){
					mycount++;
					bot = top + 32;
//printf("find A %d B %d C %d\n",A,B,X);
				}
				if(X<Y){
					top = r;
				}
				if(X>Y){
					bot = r;
				}
			}
			//phase 2
			bot = bot*m/32;
			top = top*m/32 -1;
			while(top>=bot){
				r = (top+bot)/2;
				Y = a[r];
				if(X==Y){
					mycount++;
//printf("find A %d B %d C %d\n",A,B,X);
				}
				if(X<=Y){
					top = r-1;
				}
				if(X>=Y){
					bot = r+1;
				}
			}
			j += 32;
		
		}
//		tid += GPU_NUM* blockDim.x*gridDim.x/32;
		tid += blockDim.x*gridDim.x/32;
		__syncthreads();
	}

	__syncthreads();
	//reduce
	local[threadIdx.x] = mycount;
	__syncthreads();
	if(threadIdx.x==0){
		index_t val=0;
		for(int i=0; i<blockDim.x; i++){
			val+= local[i];
		}
		count[blockIdx.x]=val;
//		count[blockIdx.x]+=val;
	}
	__syncthreads();

}

/*
__global__ void warp_binary_kernel
(	//vertex_t*	head,
	//vertex_t*	adj,
	Edge*		workload,
	vertex_t*	adj_list,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,
	index_t*	count
)
{
	//phase 1, partition
	index_t tid = (threadIdx.x + blockIdx.x * blockDim.x)/32 + Ns;
	index_t mycount=0;
	__shared__ index_t local[max_thd];

	int i = threadIdx.x%32;
	int p = threadIdx.x/32;

	while(tid<Ne){
		vertex_t A = workload[tid].A;
		vertex_t B = workload[tid].B;
		index_t m = begin[A+1]-begin[A];//degree[A];
		index_t n = begin[B+1]-begin[B];//degree[B];
//if(i==0) printf("A %d B %d\n");
		index_t temp;	
		if(m<n){
			temp = A;
			A = B;
			B = temp;
			temp = m;
			m = n;
			n = temp;
		}

		vertex_t* a = &(adj_list[begin[A]]);
		vertex_t* b = &(adj_list[begin[B]]);
		
	//initial cache
		local[p*32+i]=a[i*m/32];	
		__syncthreads();
			
	//search
		int j=i;
		while(j<n){
			vertex_t X = b[j];
			vertex_t Y;
			//phase 1: cache
			int bot = 0;
			int top = 32;
			int r;
			while(top>bot+1){
				r = (top+bot)/2;
				Y = local[p*32+r];
				if(X==Y){
					mycount++;
					bot = top + 32;
//printf("find A %d B %d C %d\n",A,B,X);
				}
				if(X<Y){
					top = r;
				}
				if(X>Y){
					bot = r;
				}
			}
			//phase 2
			bot = bot*m/32;
			top = top*m/32 -1;
			while(top>=bot){
				r = (top+bot)/2;
				Y = a[r];
				if(X==Y){
					mycount++;
//printf("find A %d B %d C %d\n",A,B,X);
				}
				if(X<=Y){
					top = r-1;
				}
				if(X>=Y){
					bot = r+1;
				}
			}
			j += 32;
		
		}
//		tid += GPU_NUM* blockDim.x*gridDim.x/32;
		tid += blockDim.x*gridDim.x/32;
		__syncthreads();
	}

	__syncthreads();
	//reduce
	local[threadIdx.x] = mycount;
	__syncthreads();
	if(threadIdx.x==0){
		index_t val=0;
		for(int i=0; i<blockDim.x; i++){
			val+= local[i];
		}
//		count[blockIdx.x]=val;
		count[blockIdx.x]+=val;
	}
	__syncthreads();

}
*/

__global__ void init_count(index_t* count)
{
	int tid = threadIdx.x;
	count[tid] = 0;
}

__global__ void reduce_kernel(index_t* count)
{
	index_t val = 0;
	for(int i=0; i<max_block; i++){
		val += count[i];
	}
	count[0] = val;
}


//---------------------------------------- cpu function--------------------
//------------------------------------------------------------------

void worker::init_gpu(){

	H_ERR(cudaSetDevice(1));

	H_ERR(cudaMalloc(&gdata.adj1, csr1->adj_size*sizeof(vertex_t)) );
	H_ERR(cudaMalloc(&gdata.beg1, csr1->beg_size*sizeof(index_t)) );

	H_ERR(cudaMalloc(&gdata.adj2, csr2->adj_size*sizeof(vertex_t)) );
	H_ERR(cudaMalloc(&gdata.beg2, csr2->beg_size*sizeof(index_t)) );
	
	H_ERR(cudaMalloc(&gdata.buffer, BufferSize*sizeof(Edge)) );
	H_ERR(cudaMalloc(&gdata.count,    max_block*sizeof(index_t)) );
	
	H_ERR(cudaDeviceSynchronize() );
	init_count <<<1,max_thd>>>(gdata.count);

	H_ERR(cudaMemcpy(gdata.adj1, csr1->adj, csr1->adj_size*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(gdata.beg1, csr1->begin, csr1->beg_size*sizeof(index_t), cudaMemcpyHostToDevice) );

	H_ERR(cudaMemcpy(gdata.adj2, csr2->adj, csr2->adj_size*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(gdata.beg2, csr2->begin, csr2->beg_size*sizeof(index_t), cudaMemcpyHostToDevice) );


	

}



void worker::call_gpu(){
	H_ERR(cudaMemcpy(gdata.buffer, buffer->edge, buffer->chunk_size*sizeof(Edge), cudaMemcpyHostToDevice) );
	H_ERR(cudaDeviceSynchronize() );

//	cout<<"csr 1 id = "<<csr1->row_id<<", offsets 1: "<<offset[csr1->row_id]<<endl;
//	cout<<"csr 2 id = "<<csr2->row_id<<", offsets 2: "<<offset[csr2->row_id]<<endl;
	//compute
	warp_binary_kernel<<<max_block,max_thd>>>
	(	
		gdata.adj1,	
		gdata.beg1,	
		gdata.adj2,	
		gdata.beg2,	
		gdata.buffer,
		buffer->chunk_size,
		offset[csr1->row_id],
		offset[csr2->row_id],
		gdata.count	
	);

	//write the result of this chunk back
	H_ERR(cudaDeviceSynchronize() );
	index_t tempcount[max_block];
	index_t mycount=0;
	H_ERR(cudaMemcpy(tempcount, gdata.count, max_block*sizeof(index_t), cudaMemcpyDeviceToHost));
	for(int i=0; i<max_block; i++){ mycount += tempcount[i];}
	gdata.host_count = mycount;

}

void worker::free_gpu(){
	H_ERR(cudaFree(gdata.adj1) );
	H_ERR(cudaFree(gdata.beg1) );
	H_ERR(cudaFree(gdata.adj2) );
	H_ERR(cudaFree(gdata.beg2) );
	H_ERR(cudaFree(gdata.buffer) );
	H_ERR(cudaFree(gdata.count) );
}




