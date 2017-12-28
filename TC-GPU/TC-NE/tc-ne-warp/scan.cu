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
		tid += blockDim.x * gridDim.x/32;
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
//		count[blockIdx.x]+=val;
		count[blockIdx.x]=val;
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
//		if(n<=0){
//			continue;
//		}

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
	cudaSetDevice(GPU_id);
	H_ERR(cudaDeviceSynchronize() );

	vertex_t*	dev_adj;
	index_t*	dev_begin;
	index_t*	dev_count;

	index_t partEdgeCount = mygraph->partEdgeCount[i];
	vertex_t vert_count = mygraph->vert_count;
	vertex_t* partAdj = mygraph->partAdj[i];
	index_t* partBegin  = mygraph->partBegin[i];
	index_t* count    = mygraph->count;

	H_ERR(cudaMalloc(&dev_adj, partEdgeCount*sizeof(vertex_t)) );
	H_ERR(cudaMalloc(&dev_begin,  (vert_count+1)*sizeof(index_t)) );
	H_ERR(cudaMalloc(&dev_count,    max_block*sizeof(index_t)) );

		
	
	index_t* block_offset;
	H_ERR(cudaMalloc(&block_offset, max_block*sizeof(index_t)) );
	
	H_ERR(cudaMemcpy(dev_adj,    partAdj, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
	H_ERR(cudaMemcpy(dev_begin,  partBegin,  (vert_count+1)*sizeof(index_t),  cudaMemcpyHostToDevice) );


//---------------- go into loop-------------------
	double time2=wtime();
	for(int j=0; j<PART_NUM; j++){	

//---------------- check -------------------------
		if(mygraph->ds_status[i * GPU_NUM + j]>0){
			break;
		}
		mygraph->ds_status[i * GPU_NUM + j]=1;
		mygraph->ds_complete[i]++;
//---------------- run ---------------------------
		index_t totalEdgeCount = mygraph->partEdgeCount[j];
		vertex_t* 	head = mygraph->partHead[j];
		vertex_t* 	adj  = mygraph->partAdj[j];

		vertex_t*	src_head;
		vertex_t*	src_adj;
		H_ERR(cudaMalloc(&src_head, totalEdgeCount*sizeof(vertex_t)) );
		H_ERR(cudaMalloc(&src_adj,  totalEdgeCount*sizeof(vertex_t)) );
		
		H_ERR(cudaMemcpy(src_adj,    adj, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
		H_ERR(cudaMemcpy(src_head,   head, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );

	double time1=wtime();
		H_ERR(cudaDeviceSynchronize() );
		warp_merge_kernel<<<max_block,max_thd>>>
		(	src_head,
			src_adj,
			dev_adj,
			dev_begin,
			0,
			totalEdgeCount,
			dev_count
		);
		H_ERR(cudaDeviceSynchronize() );
		reduce_kernel2 <<<1,1>>>(dev_count);
		H_ERR(cudaDeviceSynchronize() );
		
		H_ERR(cudaMemcpy(&count[i], dev_count, sizeof(index_t), cudaMemcpyDeviceToHost));
		thd_count += count[i];
		
		H_ERR(cudaFree(src_head) );
		H_ERR(cudaFree(src_adj) );
//		cout<<"GPU "<<i<<" part "<<j<<"\n";

//---------------- write result ---------------------------
		mygraph->ds_count[i * GPU_NUM + j] = count[i];
	}
	
//---------------- work ending ---------------------------
//---------------- work stealing ---------------------------
	int check = 0;
	for(int k=0; k<GPU_NUM; k++){
		check += mygraph->ds_complete[k];
	}
	while(check<GPU_NUM*GPU_NUM){
//step 1: looking for the GPU with most remaining work
		int min=0;
		for(int k=GPU_NUM-1; k>=0; k--){
			if(mygraph->ds_complete[k]<mygraph->ds_complete[min]){
				min = k;
			}
			if(mygraph->ds_complete[k] == mygraph->ds_complete[min]){
				if(mygraph->ds_help[k] < mygraph->ds_help[min]){
					min = k;
				}
			}
		}
		mygraph->ds_help[min]++;

		if(mygraph->ds_complete[min]==GPU_NUM){
			for(int k=0; k<GPU_NUM; k++){
				check += mygraph->ds_complete[k];
			}
			continue;
		}

//step 2: check ds_complete array and start moving data
		H_ERR(cudaFree(dev_adj) );
//		H_ERR(cudaFree(dev_head) );

		index_t partEdgeCount = mygraph->partEdgeCount[min];
		partAdj = mygraph->partAdj[min];
//		partHead= mygraph->partHead[min];
		partBegin  = mygraph->partBegin[min];

		H_ERR(cudaMalloc(&dev_adj, partEdgeCount*sizeof(vertex_t)) );
//		H_ERR(cudaMalloc(&dev_head, partEdgeCount*sizeof(vertex_t)) );

		
		H_ERR(cudaMemcpy(dev_adj,    partAdj, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
//		H_ERR(cudaMemcpy(dev_head,   partHead, partEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
		H_ERR(cudaMemcpy(dev_begin,  partBegin,  (vert_count+1)*sizeof(index_t),  cudaMemcpyHostToDevice) );
		

//step 3: set flags and work
		int j = GPU_NUM - 1;
		while(j>=0){
			if(mygraph->ds_status[min * GPU_NUM + j]>0){
				j--;
				continue;
			}
			mygraph->ds_status[min * GPU_NUM + j] = 1;
			mygraph->ds_complete[min]++;
	//work		
			index_t totalEdgeCount = mygraph->partEdgeCount[j];
			vertex_t* 	head = mygraph->partHead[j];
			vertex_t* 	adj  = mygraph->partAdj[j];

			vertex_t*	src_head;
			vertex_t*	src_adj;
		
			H_ERR(cudaMalloc(&src_head, totalEdgeCount*sizeof(vertex_t)) );
			H_ERR(cudaMalloc(&src_adj,  totalEdgeCount*sizeof(vertex_t)) );
		
			H_ERR(cudaMemcpy(src_adj,    adj, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
			H_ERR(cudaMemcpy(src_head,   head, totalEdgeCount*sizeof(vertex_t), cudaMemcpyHostToDevice) );
			

	double time1=wtime();
			H_ERR(cudaDeviceSynchronize() );
			warp_merge_kernel<<<max_block,max_thd>>>
			(	src_head,
				src_adj,
				dev_adj,
				dev_begin,
				0,
				totalEdgeCount,
				dev_count
			);
			H_ERR(cudaDeviceSynchronize() );
			
			reduce_kernel2 <<<1,1>>>(dev_count);
			H_ERR(cudaDeviceSynchronize() );
			
			H_ERR(cudaMemcpy(&count[min], dev_count, sizeof(index_t), cudaMemcpyDeviceToHost));
//---------------- write result ---------------------------
			mygraph->ds_count[min * GPU_NUM + j] = count[min];
			thd_count += count[min];
			
			H_ERR(cudaFree(src_head) );
			H_ERR(cudaFree(src_adj) );
//			cout<<"steal GPU "<<i<<" GPU "<<min<<" part "<<j<<"\n";


		}
//step 4: write result
		for(int k=0; k<GPU_NUM; k++){
			check += mygraph->ds_complete[k];
		}
	}


//---------------- work stealing end -----------------------

	double time4 = wtime();
	count[i] = thd_count;
//	cout<<"gpu "<<i<<" binary count="<<count[i]<<"\n";
//	cout<<"time = "<<time4-time2<<" seconds"<<endl;
	H_ERR(cudaFree(dev_adj) );
//	H_ERR(cudaFree(dev_head) );
	H_ERR(cudaFree(dev_begin) );
	
	H_ERR(cudaFree(block_offset) );
	H_ERR(cudaFree(dev_count) );
	return NULL;	
}


