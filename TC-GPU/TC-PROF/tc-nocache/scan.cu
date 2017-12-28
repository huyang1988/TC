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
			if(n<thread_limit && n>0){
				small_offset++;
			}
			else if(n>0){	//could be more then 2 catigories
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
			if(n<thread_limit && n>0){
				dest_head[small_offset] = head;
				dest_adj [small_offset] = adj;
				small_offset++;
			}
			else if(n>0){	//could be more then 2 catigories
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
	vertex_t*	dev_head;
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
	H_ERR(cudaMalloc(&dev_begin,  (vert_count+1)*sizeof(index_t)) );
	H_ERR(cudaMalloc(&dev_count,    max_block*sizeof(index_t)) );

		
	
	index_t* block_offset;
	H_ERR(cudaMalloc(&block_offset, max_block*sizeof(index_t)) );
	
	H_ERR(cudaMemcpy(dev_adj,    partAdj, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(dev_head,   partHead, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
//	H_ERR(cudaMemcpy(dev_degree, partDegree, vert_count*sizeof(index_t),  cudaMemcpyHostToDevice) );
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
		vertex_t*	classified_head;
		vertex_t*	classified_adj;
		
		index_t*	small_num;
		index_t*	mid_num;
		index_t*	large_num;

		vertex_t*	src_head;
		vertex_t*	src_adj;
		
		H_ERR(cudaMalloc(&small_num, max_thd*max_block*sizeof(index_t)) );
		H_ERR(cudaMalloc(&mid_num,   max_thd*max_block*sizeof(index_t)) );
		H_ERR(cudaMalloc(&large_num, max_thd*max_block*sizeof(index_t)) );
		H_ERR(cudaMalloc(&src_head, totalEdgeCount*sizeof(vertex_t)) );
		H_ERR(cudaMalloc(&src_adj,  totalEdgeCount*sizeof(vertex_t)) );
		
		H_ERR(cudaMemcpy(src_adj,    adj, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
		H_ERR(cudaMemcpy(src_head,   head, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
		
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




		warp_binary_kernel<<<max_block,max_thd>>>
		(	classified_head,
			classified_adj,
			dev_adj,
			dev_begin,
			0,
			N1,
			
			counter_1,
			counter_2,

			dev_count
		);
		H_ERR(cudaDeviceSynchronize() );

		block_binary_kernel<<<max_block,max_thd>>>
		(	classified_head,
			classified_adj,
			dev_adj,
			dev_begin,
			N1,
			N1+N2,
	//		0 + GPU_id*256,
	//		totalEdgeCount,
			
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
	
	
		H_ERR(cudaFree(small_num) );
		H_ERR(cudaFree(large_num) );
		H_ERR(cudaFree(classified_head) );
		H_ERR(cudaFree(classified_adj) );
		H_ERR(cudaFree(src_head) );
		H_ERR(cudaFree(src_adj) );
//		H_ERR(cudaFree(src_begin) );
		cout<<"GPU "<<i<<" part "<<j<<"\n";
	}
	
	double time4 = wtime();
	count[i] = thd_count;
	cout<<"gpu "<<i<<" binary count="<<count[i]<<"\n";
	cout<<"time = "<<time4-time2<<" seconds"<<endl;

	cout<<"counter for mem_read   = "<<counter_1_cpu<<endl;
//	cout<<"counter for divergence = "<<counter_2_cpu<<endl;
	total_count += counter_1_cpu;

	H_ERR(cudaFree(dev_adj) );
	H_ERR(cudaFree(dev_head) );
//	H_ERR(cudaFree(dev_degree) );
	H_ERR(cudaFree(dev_begin) );
	
	H_ERR(cudaFree(block_offset) );
	H_ERR(cudaFree(dev_count) );
	return NULL;	
}


