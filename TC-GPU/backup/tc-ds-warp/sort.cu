//sort.cu
#include "comm.h"
#include "graph.h"
#include "iostream"
#define T 56
using namespace std;


void printGraph(vertex_t vertCount, 
		vertex_t* head, 
		vertex_t* adj, 
		index_t* begin){
	for(vertex_t i=0; i<vertCount; i++){
		if(begin[i+1]>begin[i]){
			cout<<begin[i]<<" "<<begin[i+1]-begin[i]<<": ";
		}
//		for(int j=0; j<degree[i]; j++){
		for(vertex_t j=0; j<begin[i+1]-begin[i]; j++){
			cout<<head[begin[i]+j]<<"-"<<adj[begin[i]+j]<<" ";
		}
		if(begin[i+1]>begin[i]){
			cout<<"\n";
		}
	}
}

void quickSort(vertex_t * arr, index_t left, index_t right) {
      index_t i = left, j = right;
      vertex_t tmp;
      vertex_t pivot = arr[(left + right) / 2];
 
      /* partition */
      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      };
 
      /* recursion */
      if (left < j)
            quickSort(arr, left, j);
      if (i < right)
            quickSort(arr, i, right);
}


//sort and trim it to upper triangular
void graph::sort(){
	for(vertex_t i=0; i<vert_count; i++){
		index_t a=beg_pos[i];
		index_t b=beg_pos[i+1]-1;
		quickSort(adj_list,a,b);
	}
	
	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
	index_t k=0;
	for(vertex_t i=0; i<vert_count; i++){
		upperBegin[i+1]=upperBegin[i];//upperDegree[i]=0;
		index_t j=beg_pos[i];
		while(j<beg_pos[i+1]){
			if(adj_list[j]==adj_list[j+1]&&head_list[j]==head_list[j+1])
			{
				j++;
				continue;
			}
			if(head_list[j]<adj_list[j]){
				//upperHead[k] =adj_list[j];
				//upperAdj[k] =head_list[j];
				k++;
				upperBegin[i+1]++;//upperDegree[i]++;
			}
			j++;
		}
	}
	
	upperEdgeCount = k;
	upperAdj	= new vertex_t[upperEdgeCount];
	upperHead	= new vertex_t[upperEdgeCount];
//	upperDegree	= new index_t[vert_count];
//	int k=0;
	k=0;
	for(vertex_t i=0; i<vert_count; i++){
//		upperBegin[i+1]=0;//upperDegree[i]=0;
		index_t j=beg_pos[i];
		while(j<beg_pos[i+1]){
			if(adj_list[j]==adj_list[j+1]&&head_list[j]==head_list[j+1])
			{
				j++;
				continue;
			}
			if(head_list[j]<adj_list[j]){
				upperHead[k] =head_list[j];
				upperAdj[k] =adj_list[j];
				k++;
//				upperBegin[i+1]++;//upperDegree[i]++;
			}
			j++;
		}
	}
	
//	upperEdgeCount = k;
	cout<<"upper Edge Count= "<<upperEdgeCount<<"\n";
//	upperBegin[0] = 0;
//	for(int i=0; i<vert_count;i++){
//		upperBegin[i+1] += upperBegin[i];// + upperDegree[i-1];
//	}

}

void graph::reduce(){
	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
	index_t k=0;
	for(vertex_t i=0; i<vert_count; i++){
		upperBegin[i+1]=upperBegin[i];//upperDegree[i]=0;
		index_t j=beg_pos[i];
		while(j<beg_pos[i+1]){
			if(head_list[j]<adj_list[j]){
				k++;
				upperBegin[i+1]++;//upperDegree[i]++;
			}
			j++;
		}
	}
	
	upperEdgeCount = k;
	upperAdj	= new vertex_t[upperEdgeCount];
	upperHead	= new vertex_t[upperEdgeCount];
	k=0;
	for(vertex_t i=0; i<vert_count; i++){
		index_t j=beg_pos[i];
		while(j<beg_pos[i+1]){
			if(head_list[j]<adj_list[j]){
				upperHead[k] =head_list[j];
				upperAdj[k] =head_list[j];
				k++;
			}
			j++;
		}
	}
	
	cout<<"upper Edge Count= "<<upperEdgeCount<<"\n";
}


void graph::reverse_rank_by_degree(){
	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
	index_t k=0;
	for(vertex_t i=0; i<vert_count; i++){
		upperBegin[i+1]=upperBegin[i];//upperDegree[i]=0;
		index_t j=beg_pos[i];
			vertex_t h=head_list[j];
			index_t dh=beg_pos[h+1]-beg_pos[h];
		while(j<beg_pos[i+1]){
			vertex_t a=adj_list[j];
			index_t da=beg_pos[a+1]-beg_pos[a];
			if(dh>da||(dh==da && h<a)){
				k++;
				upperBegin[i+1]++;//upperDegree[i]++;
			}
			j++;
		}
	}
	
	upperEdgeCount = k;
	upperAdj	= new vertex_t[upperEdgeCount];
	upperHead	= new vertex_t[upperEdgeCount];
	k=0;
	for(vertex_t i=0; i<vert_count; i++){
		index_t j=beg_pos[i];
			vertex_t h=head_list[j];
			index_t dh=beg_pos[h+1]-beg_pos[h];
		while(j<beg_pos[i+1]){
			vertex_t a=adj_list[j];
			index_t da=beg_pos[a+1]-beg_pos[a];
			if(dh>da||(dh==da && h<a)){
				upperAdj[k] =adj_list[j];
				upperHead[k] =head_list[j];
				k++;
			}
			j++;
		}
	}
	
	cout<<"upper Edge Count= "<<upperEdgeCount<<"\n";
}

//rank-by-degree with trim
void graph::rank_by_degree(){
	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
//#pragma omp parallel for num_threads(56) schedule(static)
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(vertex_t i=0; i<vert_count; i++){
//		upperBegin[i+1]=upperBegin[i];//upperDegree[i]=0;
		upperBegin[i+1]=0;
		index_t j=beg_pos[i];
			vertex_t h=head_list[j];
			index_t dh=beg_pos[h+1]-beg_pos[h];
		while(j<beg_pos[i+1]){
			vertex_t a=adj_list[j];
			index_t da=beg_pos[a+1]-beg_pos[a];
			if(dh>da || (dh==da && h>a)){
				upperBegin[i+1]++;//upperDegree[i]++;
			}
			j++;
		}
	}
	
	for(vertex_t i=0; i<vert_count; i++){
		upperBegin[i+1] += upperBegin[i];//upperDegree[i]=0;
	}

	upperEdgeCount = upperBegin[vert_count];//k;
	upperAdj	= new vertex_t[upperEdgeCount];
	upperHead	= new vertex_t[upperEdgeCount];
//#pragma omp parallel for num_threads(56) schedule(static)
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(vertex_t i=0; i<vert_count; i++){
		index_t j=beg_pos[i];
		index_t jj=upperBegin[i];
		vertex_t h=head_list[j];
		index_t dh=beg_pos[h+1]-beg_pos[h];
		while(j<beg_pos[i+1]){
			vertex_t a=adj_list[j];
			index_t da=beg_pos[a+1]-beg_pos[a];
			if(dh>da || (dh==da && h>a)){
				upperAdj[jj] =adj_list[j];
				upperHead[jj] =head_list[j];
				jj++;//k++;
			}
			j++;
		}
	}
	
	cout<<"upper Edge Count= "<<upperEdgeCount<<"\n";
}
/*function to search the begin position to find proper place to cut adjacent list
 return is the smallest position in the data array that value is equal or larger then lookup x
 */
vertex_t BinarySearch(vertex_t x, vertex_t*A, vertex_t bot, vertex_t top){

//	for(int i=bot;i<=top;i++){
//		cout<<A[i]<<" ";
//	}
//	cout<<"\n";

	vertex_t r= (bot+top)/2;
//	int result;
	while(top>bot){
		if(x<A[r]){
			top = r;
		}
		else if(x>A[r]){
			bot = r+1;
		}
		else if(x==A[r]){
			break;
		}
		r = (bot+top)/2;
	}
	return r;
}
vertex_t BinarySearch(index_t x, index_t*A, vertex_t bot, vertex_t top){

//	for(int i=bot;i<=top;i++){
//		cout<<A[i]<<" ";
//	}
//	cout<<"\n";

	vertex_t r= (bot+top)/2;
//	int result;
	while(top>bot){
		if(x<A[r]){
			top = r;
		}
		else if(x>A[r]){
			bot = r+1;
		}
		else if(x==A[r]){
			break;
		}
		r = (bot+top)/2;
	}
	return r;
}

void graph::partition(){
	//step 1, evenly cut the upper CSR by using binary search in upperBegin
	partAdj  = new vertex_t*[PART_NUM];
	partHead = new vertex_t*[PART_NUM];
	partBegin  = new index_t*[PART_NUM];
	partEdgeCount = new index_t[PART_NUM];
	index_t offset[PART_NUM+1];
	offset[0] = 0;
	offset[PART_NUM] = upperEdgeCount;
	for(int i=1; i<PART_NUM; i++){
		index_t k=i*upperEdgeCount/PART_NUM;
		vertex_t index = BinarySearch(k, upperBegin, 0, vert_count-1);
		offset[i] = upperBegin[index];
	}
	vertex_t *destAdj = new vertex_t[upperEdgeCount];
	vertex_t *destHead= new vertex_t[upperEdgeCount];
	//we are going to use space on origin int pointer upperAdj and upperHead as output
	//thus we need two new array to store input data.
	
	//step 2, exchange two end points of each edge, get the lower CSC
	for(int i=0;i<PART_NUM;i++){
		
		vertex_t* tempHead = &upperAdj[offset[i]];
		vertex_t* tempAdj  = &upperHead[offset[i]];
		partHead[i] = &destHead[offset[i]];
		partAdj[i]  = &destAdj[offset[i]];

		partEdgeCount[i] = offset[i+1]-offset[i];

		partBegin[i] =  new index_t[vert_count+1];
		for(vertex_t j=0; j<vert_count+1; j++){
			partBegin[i][j]=0;
		}
		index_t *partOffset = new index_t[vert_count];
	//step 3, re-organize 1: go through new head list once to get new lowerDegree
		for(index_t j=0; j<partEdgeCount[i]; j++){
			//partDegree
			vertex_t head = tempHead[j];
			partBegin[i][head+1]++;

		}

	//step 4, re-organize 2: prefix, input is degree, output is begin position
		partBegin[i][0]=0;
		for(vertex_t j=0; j<vert_count; j++){
			partBegin[i][j+1] += partBegin[i][j];
			partOffset[j]=0;
		}
//		memset(partOffset, 0, vert_count*sizeof(index_t));

	//step 5, re-organize 3: go through again and moving data to transfer CSC to CSR
		for(index_t j=0; j<partEdgeCount[i]; j++){
			vertex_t head = tempHead[j];
			partHead[i][partBegin[i][head]+partOffset[head]] = tempHead[j];
			partAdj [i][partBegin[i][head]+partOffset[head]] = tempAdj [j];
			partOffset[head]++;
		}
//		printGraph(vert_count, partHead[i], partAdj[i], partBegin[i]);

	}


}

void graph::preproc(){
	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
	index_t k=0;
	index_t*inBegin = new index_t[vert_count+1];
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(vertex_t i=0; i<vert_count+1; i++){
		upperBegin[i]=0;
		inBegin[i]=0;

	}
//step 1: read round 1, to get the in-degree after orientation
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(vertex_t i=0; i<vert_count; i++){
//		upperBegin[i+1]=0;
		index_t j=beg_pos[i];
			vertex_t h=head_list[j];
			index_t dh=beg_pos[h+1]-beg_pos[h];
		while(j<beg_pos[i+1]){
			vertex_t a=adj_list[j];
			index_t da=beg_pos[a+1]-beg_pos[a];
			if(dh<da || (dh==da && h<a)){
//__sync_add_and_fetch(&k,1);
				upperBegin[i+1]++;//upperDegree[i]++;
				// to build rank-by-degree
//				__sync_fetch_and_add( &(inBegin[j+1]) , 1);
			}
			else{// if(dh>da || (dh==da && h>a)){
				inBegin[i+1]++;
			}
			j++;
		}
	}
	
	for(vertex_t i=0; i<vert_count; i++){
		upperBegin[i+1] += upperBegin[i];//upperDegree[i]=0;
		inBegin[i+1] += inBegin[i];//upperDegree[i]=0;
	}

	upperEdgeCount = upperBegin[vert_count];//k;
	upperAdj	= new vertex_t[upperEdgeCount];
	upperHead	= new vertex_t[upperEdgeCount];
//cout<<"test sycn_add_and_fetch k= "<<k<<endl;

//step 2: binary search in-degree for partition
	partAdj  = new vertex_t*[PART_NUM];
	partHead = new vertex_t*[PART_NUM];
	partBegin  = new index_t*[PART_NUM];
	partEdgeCount = new index_t[PART_NUM];
	index_t offset[PART_NUM+1];	// the vertex count begin from partition i
	vertex_t cutpoint[PART_NUM+1];	// the colum value begin from partition i
	offset[0] = 0;	
	offset[PART_NUM] = upperEdgeCount;
	cutpoint[0] = 0;
	cutpoint[PART_NUM] = vert_count+1;

//#pragma omp parallel for num_threads(PART_NUM) schedule(static)
	for(int i=1; i<PART_NUM; i++){
		index_t K=i*upperEdgeCount/PART_NUM;
		cout<<"K="<<K<<"\n";
		vertex_t index = BinarySearch(K, inBegin, 0, vert_count-1); // binary search
		cutpoint[i] = index; 			//used by each neigbhor list to find a place to cut	
		offset[i] = inBegin[index];
	}

//#pragma omp parallel for num_threads(PART_NUM) schedule(static)
	for(int i=0;i<PART_NUM;i++){
		partBegin[i] =  new index_t[vert_count+1];
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
		for(vertex_t n=0; n<vert_count+1; n++){
			partBegin[i][n]=0;
		}
	}


//step 3
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(vertex_t i=0; i<vert_count; i++){
		index_t j=beg_pos[i];
		index_t jj=upperBegin[i];
//		vertex_t h=head_list[j];
		vertex_t h=i;
		index_t dh=beg_pos[h+1]-beg_pos[h];
		//collect begin position for each partition
		vertex_t voffset[PART_NUM+1];
		voffset[0] = 0;
		voffset[PART_NUM]=dh;
		for(int n=0; n<PART_NUM; n++){
			voffset[n]=BinarySearch(cutpoint[n], &adj_list[j], 0, dh);
		}
		for(int n=0; n<PART_NUM; n++){
			for(int nn=voffset[n]; nn<voffset[n+1]; nn++ ){
				vertex_t a = adj_list[j+nn];
				index_t da=beg_pos[a+1]-beg_pos[a];
				if(dh<da || (dh==da && h<a)){
					partBegin[n][i+1]++;
					upperAdj[jj]  = adj_list[j+nn];
					upperHead[jj] = i;//head_list[j+nn];
					jj++;//k++;
				}
			}

		}
		
	}
	
//#pragma omp parallel for num_threads(PART_NUM) schedule(static)
	for(int i=0;i<PART_NUM;i++){
		for(vertex_t j=0; j<vert_count; j++){
			partBegin[i][j+1] += partBegin[i][j];//upperDegree[i]=0;
		}
	}

	for(int i=0;i<PART_NUM;i++){
		partEdgeCount[i] = partBegin[i][vert_count];		//set the edge number of each partition
		partAdj[i]   =  new vertex_t[partEdgeCount[i]];		//allocate space for each partition
	}

	
//step 4: moving partition data	
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(vertex_t i=0; i<vert_count; i++){
		index_t j=beg_pos[i];
		index_t jj=upperBegin[i];
		vertex_t h=head_list[j];
		index_t dh=beg_pos[h+1]-beg_pos[h];
		//collect begin position for each partition
		vertex_t voffset[5];
		voffset[0] = 0;
		voffset[PART_NUM]=dh;
		for(int n=0; n<PART_NUM; n++){
			voffset[n]=BinarySearch(cutpoint[n], &adj_list[j], 0, dh);
		}
		//build oriented graph
		for(int n=0; n<PART_NUM; n++){
			index_t nnn=partBegin[n][i];
			for(int nn=voffset[n]; nn<voffset[n+1]; nn++ ){
				vertex_t a = adj_list[j+nn];
				index_t da=beg_pos[a+1]-beg_pos[a];
				if(dh<da || (dh==da && h<a)){
					partAdj[n][nnn] = a;
					nnn++;	
				}
			}

		}

	}
	cout<<"upper Edge Count= "<<upperEdgeCount<<"\n";
}
