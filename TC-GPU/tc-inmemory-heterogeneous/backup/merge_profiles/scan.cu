//scan.cu
//#include "kernel.cu"
#include "comm.h"
#include "wtime.h"
#include "iostream"
#define max_thd 256
#define max_block 256 
#define thread_limit 256 
#define block_limit 1024 
using namespace std;

	
__global__ void naive_merge_kernel
(	vertex_t*	head,
	vertex_t*	adj,
	vertex_t*	adj_list,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,
	index_t*	count
)
{
	//phase 1, partition
	index_t tid = Ns + (threadIdx.x + blockIdx.x * blockDim.x)/32;
	int i = threadIdx.x%32;
	int p = threadIdx.x/32;
	long int mycount=0;
	__shared__ index_t local[max_thd];
	__shared__ vertex_t A_diag[33*8];
	__shared__ vertex_t B_diag[33*8];	
	while(tid<Ne){
		vertex_t A = head[tid];
		vertex_t B = adj[tid];
		index_t m = begin[A+1]-begin[A];//degree[A];
		index_t n = begin[B+1]-begin[B];//degree[B];

		vertex_t* a = &(adj_list[begin[A]]);
		vertex_t* b = &(adj_list[begin[B]]);
		
		if(i==0){
			A_diag[p*33+32]=m;
			B_diag[p*33+32]=n;
		}
		index_t index = (m+n)/32*i;
		vertex_t A_top, A_bottom, B_top, Ai, Bi;
		if(index>m){
			A_top = m;
			B_top = index-m;
		}
		else if(index<=m){
			A_top = index;
			B_top = 0;
		}
		if(index>n){
			A_bottom = index-n;
		}
		else if(index<=n){
			A_bottom = 0;
		}

		while(1){
			int offset=(A_top-A_bottom)/2;
			if(A_top==A_bottom){
				A_diag[p*33+i]=A_top;
				B_diag[p*33+i]=B_top;
				break;
			}

			Ai = A_top - offset;
			Bi = B_top + offset;
			if(offset<1){
				if(a[Ai-1]<b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi;
					break;
				}
				else if(a[Ai-1]>b[Bi]){
					A_diag[p*33+i]=Ai-1;
					B_diag[p*33+i]=Bi+1;
					break;
				}
				else if(a[Ai-1]==b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi+1;
					break;
				}
			}

			if(a[Ai]>b[Bi-1]){
				if(a[Ai-1]<b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi;
					break;
				}
				else if(a[Ai-1]>b[Bi]){
					A_top = Ai-1;
					B_top = Bi+1;
				}
				else if(a[Ai-1]==b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi+1;
					break;
				}
			}
			else if(a[Ai]<b[Bi-1]){
				A_bottom = Ai+1;
			}
			else if(a[Ai]==b[Bi-1]){
				A_diag[p*33+i]=Ai+1;
				B_diag[p*33+i]=Bi;
				break;
			}
		}

//		__syncthreads();

		vertex_t lowA  = A_diag[p*33+i];
		vertex_t lowB  = B_diag[p*33+i];
		vertex_t highA = A_diag[p*33+i+1];
		vertex_t highB = B_diag[p*33+i+1];
		vertex_t x,y;
		while(lowA<highA && lowB<highB){
			x=a[lowA];
			y=b[lowB];
			if(x<y){
				lowA++;
			}
			else if(x>y){
				lowB++;
			}
			else if(x==y){
				lowA++;
				lowB++;
				mycount++;
			}
		}
//		tid += blockDim.x * gridDim.x/32;
		tid += gridDim.x*blockDim.x/32;
		
//		__syncthreads();
	}
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
__global__ void warp_merge_kernel
(	vertex_t*	head,
	vertex_t*	adj,
	vertex_t*	adj_list,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,
	index_t*	count
)
{
	//phase 1, partition
	index_t tid = Ns + (threadIdx.x + blockIdx.x * blockDim.x)/32;
	int i = threadIdx.x%32;
	int p = threadIdx.x/32;
	long int mycount=0;
	__shared__ index_t local[max_thd];
	__shared__ vertex_t A_diag[33*8];
	__shared__ vertex_t B_diag[33*8];	
	while(tid<Ne){
		vertex_t A = head[tid];
		vertex_t B = adj[tid];
		index_t m = begin[A+1]-begin[A];//degree[A];
		index_t n = begin[B+1]-begin[B];//degree[B];

		vertex_t* a = &(adj_list[begin[A]]);
		vertex_t* b = &(adj_list[begin[B]]);
		
		if(i==0){
			A_diag[p*33+32]=m;
			B_diag[p*33+32]=n;
		}
		index_t index = (m+n)/32*i;
		vertex_t A_top, A_bottom, B_top, Ai, Bi;
		if(index>m){
			A_top = m;
			B_top = index-m;
		}
		else if(index<=m){
			A_top = index;
			B_top = 0;
		}
		if(index>n){
			A_bottom = index-n;
		}
		else if(index<=n){
			A_bottom = 0;
		}

		while(1){
			int offset=(A_top-A_bottom)/2;
			if(A_top==A_bottom){
				A_diag[p*33+i]=A_top;
				B_diag[p*33+i]=B_top;
				break;
			}

			Ai = A_top - offset;
			Bi = B_top + offset;
			if(offset<1){
				if(a[Ai-1]<b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi;
					break;
				}
				else if(a[Ai-1]>b[Bi]){
					A_diag[p*33+i]=Ai-1;
					B_diag[p*33+i]=Bi+1;
					break;
				}
				else if(a[Ai-1]==b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi+1;
					break;
				}
			}

			if(a[Ai]>b[Bi-1]){
				if(a[Ai-1]<b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi;
					break;
				}
				else if(a[Ai-1]>b[Bi]){
					A_top = Ai-1;
					B_top = Bi+1;
				}
				else if(a[Ai-1]==b[Bi]){
					A_diag[p*33+i]=Ai;
					B_diag[p*33+i]=Bi+1;
					break;
				}
			}
			else if(a[Ai]<b[Bi-1]){
				A_bottom = Ai+1;
			}
			else if(a[Ai]==b[Bi-1]){
				A_diag[p*33+i]=Ai+1;
				B_diag[p*33+i]=Bi;
				break;
			}
		}

//		__syncthreads();

		vertex_t lowA  = A_diag[p*33+i];
		vertex_t lowB  = B_diag[p*33+i];
		vertex_t highA = A_diag[p*33+i+1];
		vertex_t highB = B_diag[p*33+i+1];
		vertex_t x,y;
		while(lowA<highA && lowB<highB){
			x=a[lowA];
			y=b[lowB];
			if(x<y){
				lowA++;
			}
			else if(x>y){
				lowB++;
			}
			else if(x==y){
				lowA++;
				lowB++;
				mycount++;
			}
		}
//		tid += blockDim.x * gridDim.x/32;
		tid += gridDim.x*blockDim.x/32;
		
//		__syncthreads();
	}
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

__global__ void block_binary_kernel
(	vertex_t*	head,
	vertex_t*	adj,
	vertex_t*	adj_list,
	index_t*	begin,
	index_t	Ns,
	index_t	Ne,
	index_t*	count
)
{
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
	
		local[i]=a[i*m/max_thd];	
		__syncthreads();

	//search
		int j=i;
		while(j<n){
			vertex_t X = b[j];
			vertex_t Y;
			//phase 1: cache
			int bot = 0;
			int top = max_thd;
			int r;
			while(top>bot+1){
				r = (top+bot)/2;
				Y = local[r];
				if(X==Y){
					mycount++;
					bot = top + max_thd;
				}
				if(X<Y){
					top = r;
				}
				if(X>Y){
					bot = r;
				}
			}
			//phase 2
			bot = bot*m/max_thd;
			top = top*m/max_thd -1;
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
			j += max_thd;
		
		}
		tid += gridDim.x*blockDim.x/ max_thd;
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

void graph::scan(){

	cudaSetDevice(1);
	vertex_t*	dev_adj;
	vertex_t*	dev_head;
	index_t*	dev_begin;
	index_t*	dev_count;

	H_ERR(cudaMalloc(&dev_begin,  (vert_count+1)*sizeof(index_t)) );
	H_ERR(cudaMalloc(&dev_count,    max_block*sizeof(index_t)) );
	
	index_t* block_offset;
	H_ERR(cudaMalloc(&block_offset, max_block*sizeof(index_t)) );
	

	vertex_t*	src_head;
	vertex_t*	src_adj;
	
	H_ERR(cudaMalloc(&src_head, upperEdgeCount*sizeof(vertex_t)) );
	H_ERR(cudaMalloc(&src_adj,  upperEdgeCount*sizeof(vertex_t)) );
	
	H_ERR(cudaMemcpy(src_adj,    upperAdj, upperEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(src_head,   upperHead, upperEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(dev_begin,   upperBegin, (vert_count+1)*sizeof(index_t), cudaMemcpyHostToDevice) );

	dev_adj = src_adj;
	dev_head= src_head;
//		H_ERR(cudaMemcpy(src_degree, degree, vert_count*sizeof(index_t), cudaMemcpyHostToDevice) );
	
	//

double time1=wtime();


double time2=wtime();


	warp_merge_kernel<<<max_block,max_thd>>>
	(	dev_head,
		dev_adj,
		src_adj,
//			dev_degree,
		dev_begin,
		0,
		upperEdgeCount,
		dev_count
	);
	H_ERR(cudaDeviceSynchronize() );
	H_ERR(cudaDeviceSynchronize() );
	
	reduce_kernel2 <<<1,1>>>(dev_count);
	H_ERR(cudaDeviceSynchronize() );
	
	H_ERR(cudaMemcpy(&count[0], dev_count, sizeof(index_t), cudaMemcpyDeviceToHost));

double time4 = wtime();
	cout<<"total count = "<<count[0]<<endl;
	cout<<"GPU time = "<<time4-time2<<" seconds"<<endl;
	
	H_ERR(cudaFree(src_head) );
	H_ERR(cudaFree(src_adj) );

	H_ERR(cudaFree(dev_begin) );
	
	H_ERR(cudaFree(dev_count) );
}


