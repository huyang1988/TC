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

/*function to search the begin position to find proper place to cut adjacent list
 return is the largest position in the data array that value is equal or smaller then lookup x
 */
vertex_t BS(index_t x, index_t*A, vertex_t bot, vertex_t top){

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
		else if(x>=A[r]){
			bot = r+1;
		}
		else if(x==A[r]){
			while(x==A[r] && r<=top){
				r++;
			}
			r--;
			break;
		}
		r = (bot+top)/2;
	}
	return r;
}






void graph::vertical_partition(){

	//get the vertex range for partitions, step: normal ones, rem: last one.

	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
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
//	upperAdj	= new vertex_t[upperEdgeCount];
//	upperHead	= new vertex_t[upperEdgeCount];
	OrientedEdge	= new Edge[upperEdgeCount];
//cout<<"test sycn_add_and_fetch k= "<<k<<endl;

//step 2: binary search in-degree for partition
	partAdj  = new vertex_t*[PART_NUM];
//	partHead = new vertex_t*[PART_NUM];
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
cout<<"cut "<<inBegin[index]<<" at "<<index<<endl;
	}

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
//					upperAdj[jj]  = adj_list[j+nn];
//					upperHead[jj] = i;//head_list[j+nn];
					OrientedEdge[jj].B  = adj_list[j+nn];
					OrientedEdge[jj].A = i;//head_list[j+nn];
					jj++;//k++;
				}
			}

		}
		
	}
	
#pragma omp parallel for num_threads(PART_NUM) schedule(static)
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
		vertex_t h=head_list[j];
		index_t dh=beg_pos[h+1]-beg_pos[h];
		//collect begin position for each partition
		vertex_t voffset[PART_NUM+1];
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

//finished. initiate data structure for dynamic scheduling
	for(int j=0;j<PART_NUM;j++){

		stringstream ss;
		ss << "-";
		ss << j;
		string str = ss.str();
		
		string s_begin = "begin"+str;
		string s_adj = "adjacent"+str;

		char* begin_file = const_cast<char*>(s_begin.c_str());
		char* adj_file = const_cast<char*>(s_adj.c_str());
	

		FILE * a_File;
		a_File = fopen (adj_file, "wb");
		fwrite (partAdj[j] , sizeof(vertex_t), partEdgeCount[j], a_File);
		fclose (a_File);

		FILE * b_File;
		b_File = fopen (begin_file, "wb");
		fwrite (partBegin[j] , sizeof(index_t), vert_count+1, b_File);
		fclose (b_File);

	}
	string s_edge = "edge";		
	char* edge_file = const_cast<char*>(s_edge.c_str());
	FILE * a_File;
	a_File = fopen (edge_file, "wb");
	fwrite (OrientedEdge , sizeof(Edge), upperEdgeCount, a_File);
	fclose (a_File);
//
	
	

	int metadata[2];
	metadata[0] = PART_NUM;
	metadata[1] = 1;	
	string s_meta = "metadata";		

	char* meta_file = const_cast<char*>(s_meta.c_str());

	FILE * m_File;
	m_File = fopen (meta_file, "wb");
	fwrite (metadata , sizeof(int), 2 , m_File);
	fclose (m_File);
	
	/*
//verify
	for(int i=0; i<PART_NUM; i++){
		cout<<"part "<<i<<endl;
		cout<<"part size "<<partEdgeCount[i]<<endl;
//		for(vertex_t j = 0; j<5; j++){
		for(vertex_t j = 0; j<vert_count+1; j++){
			cout<<partBegin[i][j]<<" ";
		}
		cout<<endl;
//		for(index_t j = 0; j<10; j++){
		for(index_t j = 0; j<partEdgeCount[i]; j++){
			cout<<partAdj[i][j]<<" ";
		}
		cout<<endl;

	}
	*/
//	for(index_t i=0; i<upperEdgeCount; i++){
//		cout<<OrientedEdge[i].A<<"-"<<OrientedEdge[i].B<<" ";
//	}
//	cout<<endl;
	
}

