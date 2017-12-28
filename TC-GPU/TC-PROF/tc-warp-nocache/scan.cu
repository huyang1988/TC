//scan.cu
//#include "kernel.cu"
#include "comm.h"
#include "wtime.h"
#include "iostream"
#define max_thd 256 
#define max_block 256 
#define thread_limit 256 
#define block_limit 1024 

#define GPU_COWORKER 1 
graph * mygraph;
long	total_count;

__global__ void block_binary_kernel
(	vertex_t*	head,
	vertex_t*	adj,
	vertex_t*	adj_list,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,

	long int*	counter_1,
	long int*	counter_2,

	index_t*	count
)
{
	int p = threadIdx.x/32;
	long counter1=0;
	long counter2=0;
	//phase 1, partition
	index_t tid = Ns + (threadIdx.x + blockIdx.x * blockDim.x)/ max_thd;
	int i = threadIdx.x% max_thd;
	index_t mycount=0;
//	__shared__ vertex_t cache[256];
	__shared__ index_t local[max_thd];

	while(tid<Ne){
		vertex_t A = head[tid];
		vertex_t B = adj[tid];
		index_t m = begin[A+1]-begin[A];//degree[A];
		index_t n = begin[B+1]-begin[B];//degree[B];



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
	
//		local[i]=a[i*m/max_thd];	
		__syncthreads();

	counter1 += 8;
	//search
		int j=i;
		while(j<n){
			vertex_t X = b[j];
	counter1++;
			vertex_t Y;
			//phase 1: cache
			int bot = 0;
			int top = max_thd;
			int r;
/*
			while(top>bot+1){
				
		__syncthreads();
				warp_path[3*p]=0;
				warp_path[3*p+1]=0;
				warp_path[3*p+2]=0;
		__syncthreads();

				r = (top+bot)/2;
				Y = local[r];
				if(X==Y){
					mycount++;
					bot = top + max_thd;
					warp_path[3*p]=1;
				}
				if(X<Y){
					top = r;
					warp_path[3*p+1]=1;
				}
				if(X>Y){
					bot = r;
					warp_path[3*p+2]=1;
				}

				int k=0;
				if(warp_path[3*p]!=0){
					k++;
				}
				if(warp_path[3*p+1]!=0){
					k++;
				}
				if(warp_path[3*p+2]!=0){
					k++;
				}
		counter2 +=k;

			}
*/
			//phase 2
//			bot = bot*m/max_thd;
//			top = top*m/max_thd -1;
			bot = 0;
			top = m-1;
			while(top>=bot){
		

				r = (top+bot)/2;
				Y = a[r];
	counter1++;
				if(X==Y){
					mycount++;
				}
				if(X<=Y){
					top = r-1;
				}
				if(X>=Y){
					bot = r+1;
				}
				
			
			}
			j += max_thd;
		
		}
		tid += GPU_COWORKER * gridDim.x*blockDim.x/ max_thd;
		__syncthreads();
	}

	//reduce
	__syncthreads();
	local[threadIdx.x] = mycount;
	__syncthreads();
	if(threadIdx.x==0){
		index_t val=0;
		for(int i=0; i<blockDim.x; i++){
			val+= local[i];
		}
		count[blockIdx.x]+=val;
//		count[blockIdx.x]=val;
	}
	counter_1[blockDim.x*blockIdx.x+threadIdx.x]+=counter1;
	counter_2[blockDim.x*blockIdx.x+threadIdx.x]+=counter2;
}

__global__ void warp_binary_kernel
(	vertex_t*	head,
	vertex_t*	adj,
	vertex_t*	adj_list,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,

	long int*	counter_1,
	long int*	counter_2,

	index_t*	count
)
{
	long counter1=0;
	long counter2=0;
	//phase 1, partition
	index_t tid = (threadIdx.x + blockIdx.x * blockDim.x)/32 + Ns;
	index_t mycount=0;
	__shared__ index_t local[max_thd];

	int i = threadIdx.x%32;
	int p = threadIdx.x/32;

	while(tid<Ne){
		vertex_t A = head[tid];
		vertex_t B = adj[tid];
		index_t m = begin[A+1]-begin[A];//degree[A];
		index_t n = begin[B+1]-begin[B];//degree[B];

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
//		local[p*32+i]=a[i*m/32];	
	counter1+=8;
		__syncthreads();
			
	//search
		int j=i;
		while(j<n){
			vertex_t X = b[j];
	counter1++;
			vertex_t Y;

			//phase 1: cache
			int bot = 0;
			int top = 32;
			int r;
/*
			while(top>bot+1){
		__syncthreads();
				warp_path[3*p]=0;
				warp_path[3*p+1]=0;
				warp_path[3*p+2]=0;
		__syncthreads();
				r = (top+bot)/2;
				Y = local[p*32+r];
				if(X==Y){
					mycount++;
					bot = top + 32;
					warp_path[3*p]=1;
				}
				if(X<Y){
					top = r;
					warp_path[3*p+1]=1;
				}
				if(X>Y){
					bot = r;
					warp_path[3*p+2]=1;
				}
				int k=0;
				if(warp_path[3*p]!=0){
					k++;
				}
				if(warp_path[3*p+1]!=0){
					k++;
				}
				if(warp_path[3*p+2]!=0){
					k++;
				}
		counter2 +=k;
			}
*/
			//phase 2
			bot = 0;
			top = m -1;
			while(top>=bot){
				r = (top+bot)/2;
				Y = a[r];
	counter1++;
				if(X==Y){
					mycount++;
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
		tid += GPU_COWORKER* blockDim.x*gridDim.x/32;
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
	}
	__syncthreads();
	counter_1[blockDim.x*blockIdx.x+threadIdx.x]=counter1;
	counter_2[blockDim.x*blockIdx.x+threadIdx.x]=counter2;

}

//----------------------------------------------------------------------------------------

__global__ void reduce_kernel_count(index_t* count)
{
	index_t val = 0;
	for(int i=0; i<max_block*max_block; i++){
		val += count[i];
	}
	count[0] = val;
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

void* part_scan(void * data){

	index_t thd_count=0;	
	int GPU_id = *(int*)data;
	int i = GPU_id;
//	cout<<"GPU id = "<<GPU_id<<"\n";
	cudaSetDevice(GPU_id);
	H_ERR(cudaDeviceSynchronize() );

	vertex_t*	dev_adj;
	index_t*	dev_begin;
	index_t*	dev_count;

	index_t partEdgeCount = mygraph->partEdgeCount[i];
	vertex_t vert_count = mygraph->vert_count;
	vertex_t* partAdj = mygraph->partAdj[i];
	vertex_t* partHead= mygraph->partHead[i];
	index_t* partBegin  = mygraph->partBegin[i];
	index_t* count    = mygraph->count;

	H_ERR(cudaMalloc(&dev_adj, partEdgeCount*sizeof(vertex_t)) );
	H_ERR(cudaMalloc(&dev_begin,  (vert_count+1)*sizeof(index_t)) );
	H_ERR(cudaMalloc(&dev_count,    max_block*sizeof(index_t)) );

		
	
	index_t* block_offset;
	H_ERR(cudaMalloc(&block_offset, max_block*sizeof(index_t)) );
	
	H_ERR(cudaMemcpy(dev_adj,    partAdj, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(dev_begin,  partBegin,  (vert_count+1)*sizeof(index_t),  cudaMemcpyHostToDevice) );
	
long 	counter_1_cpu=0;
//long 	counter_2_cpu=0;
long int tmp_counter1,tmp_counter2;
long int*       counter_1;//counter for memory read
long int*       counter_2;//counter for divergence
H_ERR(cudaMalloc(&counter_1,    max_thd*max_block*sizeof(long int)) );
H_ERR(cudaMalloc(&counter_2,    max_thd*max_block*sizeof(long int)) );

	double time2=wtime();
	for(int j=0; j<PART_NUM; j++){	
		index_t totalEdgeCount = mygraph->partEdgeCount[j];
		vertex_t* 	head = mygraph->partHead[j];
		vertex_t* 	adj  = mygraph->partAdj[j];

		vertex_t*	src_head;
		vertex_t*	src_adj;
		
		H_ERR(cudaMalloc(&src_head, totalEdgeCount*sizeof(vertex_t)) );
		H_ERR(cudaMalloc(&src_adj,  totalEdgeCount*sizeof(vertex_t)) );
		
		H_ERR(cudaMemcpy(src_adj,    adj, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
		H_ERR(cudaMemcpy(src_head,   head, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
		
		//

	double time1=wtime();
		H_ERR(cudaDeviceSynchronize() );

		
		warp_binary_kernel<<<max_block,max_thd>>>
		(	src_head,
			src_adj,
			dev_adj,
			dev_begin,
			0,
			totalEdgeCount,
			
			counter_1,
			counter_2,

			dev_count
		);
		H_ERR(cudaDeviceSynchronize() );
		
		reduce_kernel <<<1,1>>>(dev_count);
		H_ERR(cudaDeviceSynchronize() );
		
		H_ERR(cudaMemcpy(&count[i], dev_count, sizeof(index_t), cudaMemcpyDeviceToHost));
		thd_count += count[i];
		
		
		reduce_kernel_count <<<1,1>>>(counter_1);
		H_ERR(cudaDeviceSynchronize() );
//		reduce_kernel_count <<<1,1>>>(counter_2);
//		H_ERR(cudaDeviceSynchronize() );
		//long int tmp_counter1,tmp_counter2;
		H_ERR(cudaMemcpy(&tmp_counter1, counter_1, sizeof(long), cudaMemcpyDeviceToHost));
//		H_ERR(cudaMemcpy(&tmp_counter2, counter_2, sizeof(long), cudaMemcpyDeviceToHost));
		counter_1_cpu += tmp_counter1;
//		counter_2_cpu += tmp_counter2;
	
	
		H_ERR(cudaFree(src_head) );
		H_ERR(cudaFree(src_adj) );
//		H_ERR(cudaFree(src_begin) );
//		cout<<"GPU "<<i<<" part "<<j<<"\n";
	}
	
	double time4 = wtime();
	count[i] = thd_count;
//	cout<<"gpu "<<i<<" binary count="<<count[i]<<"\n";
//	cout<<"time = "<<time4-time2<<" seconds"<<endl;

//	cout<<"counter for mem_read   = "<<counter_1_cpu<<endl;
//	cout<<"counter for divergence = "<<counter_2_cpu<<endl;
	total_count += counter_1_cpu;

	H_ERR(cudaFree(dev_adj) );
	H_ERR(cudaFree(dev_begin) );
	
	H_ERR(cudaFree(block_offset) );
	H_ERR(cudaFree(dev_count) );
	return NULL;	
}


