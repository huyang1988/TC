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
#define warp_thd 256 
#define block_thd 256
graph * mygraph;
//#define GPU_NUM 4 
//int     NUM_ONE_TASK;

__global__ void block_binary_kernel
(	vertex_t*	head,
	vertex_t*	adj,
	vertex_t*	adj_list,
//	index_t*	degree,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,
	index_t*	count
)
{
	//phase 1, partition
	index_t tid = Ns + (threadIdx.x + blockIdx.x * blockDim.x)/ block_thd;
	int i = threadIdx.x% block_thd;
	index_t mycount=0;
//	__shared__ int cache[256];
//	__shared__ int offset[256];
	__shared__ index_t local[block_thd];

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
	
		local[i]=a[i*m/block_thd];	
		__syncthreads();

	//search
		int j=i;
		while(j<n){
			vertex_t X = b[j];
			vertex_t Y;
			//phase 1: cache
			int bot = 0;
			int top = block_thd;
			int r;
			while(top>bot+1){
				r = (top+bot)/2;
				Y = local[r];
				if(X==Y){
					mycount++;
					bot = top + block_thd;
				}
				if(X<Y){
					top = r;
				}
				if(X>Y){
					bot = r;
				}
			}
			//phase 2
//			bot = bot*k;
//			top = top*k;
			bot = bot*m/block_thd;
			top = top*m/block_thd -1;
			while(top>=bot){
				r = (top+bot)/2;
				Y = a[r];
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
			j += block_thd;
		
		}
		tid += GPU_COWORKER * gridDim.x*blockDim.x/ block_thd;
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
}

__global__ void warp_binary_kernel
(	vertex_t*	head,
	vertex_t*	adj,
	vertex_t*	adj_list,
//	index_t*	degree,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,
	index_t*	count
)
{
	//phase 1, partition
	index_t tid = (threadIdx.x + blockIdx.x * blockDim.x)/32 + Ns;
	index_t mycount=0;
	__shared__ index_t local[warp_thd];

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
				}
				if(X<Y){
					top = r;
				}
				if(X>Y){
					bot = r;
				}
			}
			//phase 2
//			bot = bot*k;
//			top = top*k;
			bot = bot*m/32;
			top = top*m/32 -1;
			while(top>=bot){
				r = (top+bot)/2;
				Y = a[r];
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

}

//----------------------------------------------------------------------------------------

__global__ void classify_kernel	//step 1: classify the edge list into different arrays
(	vertex_t* adj_list,
	vertex_t* head_list,
	index_t* begin,
	index_t  N,		//inputs
	index_t* small_num,
	index_t* mid_num,
	index_t* large_num
	//outputs: small/large head, adjacent, and number by thread
)
{
	int tid = threadIdx.x +blockIdx.x*blockDim.x;
	index_t bin_size = (N-1)/(blockDim.x*gridDim.x)+1;
	index_t thd_base = tid*bin_size;		//start point of threads space
	index_t small_offset=0;
	index_t mid_offset=0;
	index_t large_offset=0;
	
	//temp variables
	vertex_t head;
	vertex_t adj;
	index_t m;
	index_t n;
	for(index_t i=0; i<bin_size; i++){
		index_t id = thd_base + i;
		if(id<N){
			head = head_list[id];
			adj  = adj_list[id];
			m = begin[head+1]-begin[head];//degree[head];
			n = begin[adj+1]-begin[adj];//degree[adj];
			if(m<n){
				n=m;
			}
			if(n<thread_limit){
				small_offset++;
			}
			else if(n<block_limit){	//could be more then 2 catigories
//			else{
				mid_offset++;
			}
			else {	//could be more then 2 catigories
				large_offset++;
			}
		}
	}
	small_num[tid] = small_offset;
	mid_num[tid]   = mid_offset;
	large_num[tid] = large_offset;

}

__global__ void prefix_kernel_1	//this prefix scan function could be easier for data size is always 256*256
(	
 	index_t*	data,
	index_t*	block_offset
)
{
		
	//step 1: each block do prefix sum inside
	int tid = threadIdx.x +blockIdx.x*blockDim.x;

	__shared__ index_t temp_in[256];
	temp_in[threadIdx.x] = data[tid];
	__syncthreads();

	index_t val=0;
	for(int i=0; i<=threadIdx.x; i++){
		val += temp_in[i];
	}


	__syncthreads();
	
	if(threadIdx.x==255){
		block_offset[blockIdx.x] = val;
		
	}
	data[tid] = val;
	__syncthreads();
	
}

__global__ void prefix_kernel_2	
(	
	index_t*	block_offset
)
{
	//step 2: collect each block's offset and do prefix for this set
	__shared__ index_t temp_in[256];
	temp_in[threadIdx.x] = block_offset[threadIdx.x];
	__syncthreads();
	index_t val=0;
	for(int i=0; i<threadIdx.x; i++){
		val += temp_in[i];
	}
//		val = temp_in[threadIdx.x];
	block_offset[threadIdx.x] = val;
	__syncthreads();
	
}

__global__ void prefix_kernel_3	
(	
	index_t*	data,
	index_t*	block_offset
)
{
	//step 3: update by adding block offset
	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	index_t val = data[tid];
	index_t offset = block_offset[blockIdx.x];
	val += offset;

	data[tid] = val;
	__syncthreads();
}

__global__ void collect_kernel
(	vertex_t* 	adj_list,
	vertex_t* 	head_list,
	index_t* 	begin,
	index_t	N,
	index_t* 	small_num,
	index_t* 	mid_num,
	index_t* 	large_num,
	index_t 	N1,
	index_t	N2,
	vertex_t*	dest_head,
	vertex_t*	dest_adj
)
{
	int tid = threadIdx.x +blockIdx.x*blockDim.x;
	index_t bin_size = (N-1)/(blockDim.x*gridDim.x)+1;
	index_t thd_base = tid*bin_size;		//start point of threads space


	index_t thd_base_small = 0;
	index_t thd_base_mid   = N1;
	index_t thd_base_large = N1+N2;
	if(tid!=0){
		thd_base_small = small_num[tid-1];
		thd_base_mid   = N1 + mid_num[tid-1];
		thd_base_large = N1 + N2 + large_num[tid-1];
	}
	
	//temp variables
	vertex_t head;
	vertex_t adj;
	index_t m;
	index_t n;
	index_t small_offset = thd_base_small;
	index_t mid_offset   = thd_base_mid;
	index_t large_offset = thd_base_large;
	for(index_t i=0; i<bin_size; i++){
		index_t id = thd_base + i;
		if(id<N){
			head = head_list[id];
			adj  = adj_list[id];
			m = begin[head+1]-begin[head];//degree[head];
			n = begin[adj+1]-begin[adj];//degree[adj];
			if(m<n){
				n=m;
			}
			if(n<thread_limit){
				dest_head[small_offset] = head;
				dest_adj [small_offset] = adj;
				small_offset++;
			}
			else if(n<block_limit){	//could be more then 2 catigories
//			else{
				dest_head[mid_offset] = head;
				dest_adj [mid_offset] = adj;
				mid_offset++;
			}
			else {	//could be more then 2 catigories
				dest_head[large_offset] = head;
				dest_adj [large_offset] = adj;
				large_offset++;
			}
		}
	}
}


__global__ void reduce_kernel2(index_t* count)
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
	vertex_t*	dev_head;
//	index_t*	dev_degree;
	index_t*	dev_begin;
	index_t*	dev_count;

	index_t partEdgeCount = mygraph->partEdgeCount[i];
	vertex_t vert_count = mygraph->vert_count;
	vertex_t* partAdj = mygraph->partAdj[i];
	vertex_t* partHead= mygraph->partHead[i];
//	index_t* partDegree = mygraph->partDegree[i];
	index_t* partBegin  = mygraph->partBegin[i];
	index_t* count    = mygraph->count;

	H_ERR(cudaMalloc(&dev_adj, partEdgeCount*sizeof(vertex_t)) );
	H_ERR(cudaMalloc(&dev_head, partEdgeCount*sizeof(vertex_t)) );
//	H_ERR(cudaMalloc(&dev_degree, vert_count*sizeof(index_t)) );
	H_ERR(cudaMalloc(&dev_begin,  vert_count*sizeof(index_t)) );
	H_ERR(cudaMalloc(&dev_count,    max_block*sizeof(index_t)) );

		
	
	index_t* block_offset;
	H_ERR(cudaMalloc(&block_offset, max_block*sizeof(index_t)) );
	
	H_ERR(cudaMemcpy(dev_adj,    partAdj, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(dev_head,   partHead, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
//	H_ERR(cudaMemcpy(dev_degree, partDegree, vert_count*sizeof(index_t),  cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(dev_begin,  partBegin,  vert_count*sizeof(index_t),  cudaMemcpyHostToDevice) );

	for(int j=0; j<PART_NUM; j++){	
		index_t totalEdgeCount = mygraph->partEdgeCount[j];
		vertex_t* 	head = mygraph->partHead[j];
		vertex_t* 	adj  = mygraph->partAdj[j];
//		index_t* 	degree = mygraph->partDegree[j];
		vertex_t*	classified_head;
		vertex_t*	classified_adj;
		
		index_t*	small_num;
		index_t*	mid_num;
		index_t*	large_num;

		vertex_t*	src_head;
		vertex_t*	src_adj;
//		index_t*	src_degree;
		
		H_ERR(cudaMalloc(&small_num, max_thd*max_block*sizeof(index_t)) );
		H_ERR(cudaMalloc(&mid_num,   max_thd*max_block*sizeof(index_t)) );
		H_ERR(cudaMalloc(&large_num, max_thd*max_block*sizeof(index_t)) );
		H_ERR(cudaMalloc(&src_head, totalEdgeCount*sizeof(vertex_t)) );
		H_ERR(cudaMalloc(&src_adj,  totalEdgeCount*sizeof(vertex_t)) );
//		H_ERR(cudaMalloc(&src_degree,  vert_count*sizeof(index_t)) );
		
		H_ERR(cudaMemcpy(src_adj,    adj, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
		H_ERR(cudaMemcpy(src_head,   head, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
//		H_ERR(cudaMemcpy(src_degree, degree, vert_count*sizeof(index_t), cudaMemcpyHostToDevice) );
		
		H_ERR(cudaMalloc(&classified_head, totalEdgeCount*sizeof(vertex_t)) );
		H_ERR(cudaMalloc(&classified_adj,  totalEdgeCount*sizeof(vertex_t)) );
		//

	double time1=wtime();
		H_ERR(cudaDeviceSynchronize() );

		
		classify_kernel <<<max_block,max_thd>>>(
					src_adj,
					src_head,
					dev_begin,
					totalEdgeCount,
					small_num,
					mid_num,
					large_num
					);
		H_ERR(cudaDeviceSynchronize() );

		//test for prefix sum

		prefix_kernel_1 <<<max_block,max_thd>>>(small_num, block_offset);
		H_ERR(cudaDeviceSynchronize() );
		prefix_kernel_2 <<<1,max_thd>>>(block_offset);
		H_ERR(cudaDeviceSynchronize() );
		prefix_kernel_3 <<<max_block,max_thd>>>(small_num, block_offset);
		H_ERR(cudaDeviceSynchronize() );

		prefix_kernel_1 <<<max_block,max_thd>>>(mid_num, block_offset);
		H_ERR(cudaDeviceSynchronize() );
		prefix_kernel_2 <<<1,max_thd>>>(block_offset);
		H_ERR(cudaDeviceSynchronize() );
		prefix_kernel_3 <<<max_block,max_thd>>>(mid_num, block_offset);
		H_ERR(cudaDeviceSynchronize() );
		
		prefix_kernel_1 <<<max_block,max_thd>>>(large_num, block_offset);
		H_ERR(cudaDeviceSynchronize() );
		prefix_kernel_2 <<<1,max_thd>>>(block_offset);
		H_ERR(cudaDeviceSynchronize() );
		prefix_kernel_3 <<<max_block,max_thd>>>(large_num, block_offset);
		H_ERR(cudaDeviceSynchronize() );

	index_t N1,N2,N3;	
		H_ERR(cudaMemcpy(&N1 ,  &small_num[65535] , sizeof(index_t), cudaMemcpyDeviceToHost) );
		H_ERR(cudaMemcpy(&N2 , &mid_num[65535] , sizeof(index_t), cudaMemcpyDeviceToHost) );
		H_ERR(cudaMemcpy(&N3 ,  &large_num[65535]   , sizeof(index_t), cudaMemcpyDeviceToHost) );

		H_ERR(cudaDeviceSynchronize() );
	//	cout<<"N1 = "<<N1<<"\n";
	//	cout<<"N2 = "<<N2<<"\n";
	//	cout<<"N3 = "<<N3<<"\n";
		
		collect_kernel <<<max_block,max_thd>>>(
					src_adj,
					src_head,
					dev_begin,
					totalEdgeCount,
					small_num,
					mid_num,
					large_num,
					N1,
					N2,
					classified_head,
					classified_adj
					);
		H_ERR(cudaDeviceSynchronize() );


	//double time2=wtime();


		warp_binary_kernel<<<max_block,warp_thd>>>
		(	classified_head,
			classified_adj,
			dev_adj,
//			dev_degree,
			dev_begin,
			0,
			N1,
			dev_count
		);
		H_ERR(cudaDeviceSynchronize() );

	//double time3 = wtime();
	/*
	//	int* 	large_head = &(classified_head[N1]);
	//	int*	large_adj  = &(classified_adj[N1]);
		block_binary_kernel<<<max_block,block_thd>>>
		(	classified_head,
			classified_adj,
			dev_adj,
//			dev_degree,
			dev_begin,
			N1+GPU_id*256,
			N1+N2,
			dev_count
		);
		H_ERR(cudaDeviceSynchronize() );
	*/
	//double time4 = wtime();
		block_binary_kernel<<<max_block,max_thd>>>
		(	classified_head,
			classified_adj,
			//dev_head,
			//dev_adj,

			dev_adj,
//			dev_degree,
			dev_begin,
			N1,
	//		0 + GPU_id*256,
			totalEdgeCount,
			dev_count
		);
		H_ERR(cudaDeviceSynchronize() );
		
		reduce_kernel2 <<<1,1>>>(dev_count);
		H_ERR(cudaDeviceSynchronize() );
		
		H_ERR(cudaMemcpy(&count[i], dev_count, sizeof(index_t), cudaMemcpyDeviceToHost));
		thd_count += count[i];
		
		H_ERR(cudaFree(small_num) );
		H_ERR(cudaFree(large_num) );
		H_ERR(cudaFree(classified_head) );
		H_ERR(cudaFree(classified_adj) );
		H_ERR(cudaFree(src_head) );
		H_ERR(cudaFree(src_adj) );
//		H_ERR(cudaFree(src_degree) );
		cout<<"GPU "<<i<<" part "<<j<<"\n";
	}

	count[i] = thd_count;
	cout<<"gpu binary count="<<count[i]<<"\n";
	H_ERR(cudaFree(dev_adj) );
	H_ERR(cudaFree(dev_head) );
//	H_ERR(cudaFree(dev_degree) );
	H_ERR(cudaFree(dev_begin) );
	
	H_ERR(cudaFree(block_offset) );
	H_ERR(cudaFree(dev_count) );
	return NULL;	
}


