#include "comm.h"
#include "wtime.h"
#include "iostream"
#define max_thd 256 
#define max_block 256 

__global__ void tc_warp_kernel
(
	int*	adj_list,
	int*	head_list,
	int*	adj_card,
	int*	beg_pos,
	int*      count,
	int   vert_count,
	int   edge_count
){
	int tid = blockIdx.x * blockDim.x + threadIdx.x;//
	
	__shared__ int local[max_thd];	
	int mycount=0;
	int warp_id = tid/32;
	while(warp_id < edge_count){
		int U = head_list[warp_id];
		int V = adj_list[warp_id];
		int degreeU = adj_card[U];
		int degreeV = adj_card[V];
		int workId = threadIdx.x%32;
		int workload = degreeU * degreeV;
		//using a while loop to increase workId, to make a warp coordinate for one intersection
		while(workId < workload){
			int workU = workId/degreeV;
			int workV = workId%degreeV;
			int offsetU = beg_pos[U];
			int offsetV = beg_pos[V];
			if(adj_list[offsetU+workU] == adj_list[offsetV+workV]){
				mycount++;
			}
			workId+=32;
		}
		warp_id += blockDim.x * gridDim.x/32;
	}
	//reduce
	local[threadIdx.x] = mycount;
	__syncthreads();
	mycount=0;
	if(threadIdx.x==0){
		for(int i=0; i<max_thd; i++){
			mycount+= local[i];
		}
		count[blockIdx.x]=mycount;
	}
}

__global__ void tc_block_kernel
(
	int*	adj_list,
	int*	head_list,
	int*	adj_card,
	int*	beg_pos,
	int*      count,
	int   vert_count,
	int   edge_count
){
	int tid = blockIdx.x * blockDim.x + threadIdx.x;//
	__shared__ int local[max_thd];	
	int mycount=0;
	int grid_id = blockIdx.x;
	while(grid_id < edge_count){
		int U = head_list[grid_id];
		int V = adj_list[grid_id];
		int degreeU = adj_card[U];
		int degreeV = adj_card[V];
		int offsetU = beg_pos[U];
		int offsetV = beg_pos[V];
		int workId = threadIdx.x;
		int workload = degreeU * degreeV;
		//using a while loop to increase workId, to make a warp coordinate for one intersection
		while(workId < workload){
			int workU = workId/degreeV;
			int workV = workId%degreeV;
			if(adj_list[offsetU+workU] == adj_list[offsetV+workV]){
				mycount++;
			}
			workId+=blockDim.x;
		}
		grid_id += gridDim.x;
	}
	//reduce
	local[threadIdx.x] = mycount;
	__syncthreads();
	mycount=0;
	if(threadIdx.x==0){
		for(int i=0; i<max_thd; i++){
			mycount+= local[i];
		}
		count[blockIdx.x]=mycount;
	}
	__syncthreads();
}

__global__ void reduce_kernel(int* count)
{
	int val = 0;
	for(int i=0; i<256; i++){
		val += count[i];
	}
	count[0] = val;
}



__global__ void tc_kernel
(
	int*	adj_list,
	int*	adj_card,
	int*	beg_pos,
	int*      count,
	int   vert_count,
	int   edge_count
){
 	__shared__ int thd_count[max_thd];
	int v = blockIdx.x;//vertex v
	int M = vert_count;

	
	int mycount=0;
	//step
	while(v<M){
		int local_count=0;
		int i = threadIdx.x;
		int N = adj_card[v];
		int beg = beg_pos[v];
		while(i<N){
			int w = adj_list[beg + i];//vertex w
			int s1 = beg_pos[v];	//start 1
			int e1 = beg_pos[v] + adj_card[v];	//end 1
			int s2 = beg_pos[w];
			int e2 = beg_pos[w] + adj_card[w];
/*			while(s1<e1 && s2<e2){
				if(adj_list[s1]<adj_list[s2]){
					s1++;
				}
				else if(adj_list[s1]>adj_list[s2]){
					s2++;
				}
				else if(adj_list[s1]==adj_list[s2]){
					s1++;
					s2++;
					local_count++;
				}
			}
*/
			for(int j=s1; j<e1; j++){
				for(int k=s2; k<e2; k++){
					if(adj_list[j]==adj_list[k]){
						local_count++;
						mycount++;
//						break;
					}
				}
			}

			i += blockDim.x;
		}
		thd_count[threadIdx.x]=local_count;
		//sycn
		__syncthreads();
		if(threadIdx.x==0){
			local_count=0;
			int k=0;
			for(k=0;k<max_thd;k++){
				local_count += thd_count[k];
			}
			count[v] = local_count;
		}
		v += gridDim.x;
		__syncthreads();
	}
} 



//template <typename vertex_t, typename index_t>
void graph//<vertex_t, index_t>
:: triangle_count()
{
	
	int*	dev_adj_list;
	int*	dev_head_list;
	int*	dev_adj_card;
	int*	dev_beg_pos;
	int*	dev_count;
	int*	dev_count2;

	cudaMalloc(&dev_adj_list, upperEdgeCount*sizeof(int));
	cudaMalloc(&dev_head_list, upperEdgeCount*sizeof(int));
	cudaMalloc(&dev_adj_card, vert_count*sizeof(int));
	cudaMalloc(&dev_beg_pos,  vert_count*sizeof(int));
	
	cudaMalloc(&dev_count2,		vert_count*sizeof(int));
	cudaMalloc(&dev_count,		max_thd*sizeof(int));

	cudaMemcpy(dev_adj_list, upperAdj, upperEdgeCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_head_list, upperHead, upperEdgeCount*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_adj_card, upperDegree, vert_count*sizeof(int),  cudaMemcpyHostToDevice);
	cudaMemcpy(dev_beg_pos,  upperBegin,  vert_count*sizeof(int),  cudaMemcpyHostToDevice);
	
	//
	cudaDeviceSynchronize();
	tc_warp_kernel <<<max_block,max_thd>>>(
				dev_adj_list,
				dev_head_list,
				dev_adj_card,
				dev_beg_pos,
				dev_count,
				vert_count,
				upperEdgeCount
				);
	cudaDeviceSynchronize();
	reduce_kernel <<<1,1>>>(dev_count);
//	tc_kernel <<<max_block,max_thd>>>(
//				dev_adj_list,
//				dev_adj_card,
//				dev_beg_pos,
//				dev_count2,
//				vert_count,
//				edge_count
//				);
//	cudaMemcpy(valid, dev_count, vert_count*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(count, dev_count, sizeof(int), cudaMemcpyDeviceToHost);

//	cudaFree(dev_adj_list);
//	cudaFree(dev_head_list);
//	cudaFree(dev_adj_card);
//	cudaFree(dev_beg_pos);
//	cudaFree(dev_count);
	
}
