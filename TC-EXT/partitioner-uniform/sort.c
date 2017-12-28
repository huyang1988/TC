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

//-------------------------------------------------------------
/*
void graph::preproc(){
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
//	index_t offset[PART_NUM+1];	// the vertex count begin from partition i
	vertex_t cutpoint[PART_NUM+1];	// the colum value begin from partition i
//	offset[0] = 0;	
//	offset[PART_NUM] = upperEdgeCount;
	cutpoint[0] = 0;
	cutpoint[PART_NUM] = vert_count+1;

#pragma omp parallel for num_threads(PART_NUM) schedule(static)
	for(int i=1; i<PART_NUM; i++){
		index_t K=i*upperEdgeCount/PART_NUM;
		cout<<"K="<<K<<"\n";
		vertex_t index = BinarySearch(K, inBegin, 0, vert_count-1); // binary search
		cutpoint[i] = index; 			//used by each neigbhor list to find a place to cut	
//		offset[i] = inBegin[index];
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
					OrientedEdge[jj].A  = adj_list[j+nn];
					OrientedEdge[jj].B = i;//head_list[j+nn];
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
	ChunkNum = (upperEdgeCount-1)/BufferSize + 1;
cout<<"Chunk Number = "<<ChunkNum<<endl;
	
	ds_count 	= new index_t [PART_NUM*ChunkNum];
	ds_status	= new index_t [PART_NUM*ChunkNum];
	
	for(int i=0; i<DEV_NUM; i++){
		index_t t = ChunkNum/DEV_NUM;
		if( i<( ChunkNum % DEV_NUM) ) t++;
		ds_last[i] = t ;
cout<<"ds_last: "<<i<<" = "<<t<<endl;
	}

	
}

*/





void graph::preproc(){

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
//	index_t offset[PART_NUM+1];	// the vertex count begin from partition i
	vertex_t cutpoint[PART_NUM+1];	// the colum value begin from partition i
//	offset[0] = 0;	
//	offset[PART_NUM] = upperEdgeCount;
	cutpoint[0] = 0;
	cutpoint[PART_NUM] = vert_count+1;

#pragma omp parallel for num_threads(PART_NUM) schedule(static)
	for(int i=1; i<PART_NUM; i++){
		index_t K=i*upperEdgeCount/PART_NUM;
		cout<<"K="<<K<<"\n";
		vertex_t index = BinarySearch(K, inBegin, 0, vert_count-1); // binary search
		cutpoint[i] = index; 			//used by each neigbhor list to find a place to cut	
//		offset[i] = inBegin[index];
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
//
/*
//verify
	for(int i=0; i<PART_NUM; i++){
		cout<<"part "<<i<<endl;
		for(vertex_t j = 0; j<vert_count+1; j++){
			cout<<partBegin[i][j]<<" ";
		}
		cout<<endl;
		for(index_t j = 0; j<partEdgeCount[i]; j++){
			cout<<partAdj[i][j]<<" ";
		}
		cout<<endl;

	}
*/
	
}


void graph::edge_2d(){
	
	int N = _2DPART;
	
//	vertex_t step, rem;
	step = (vert_count - 1)/N + 1;
	cout<<"step ="<<step<<endl;
	rem = vert_count-(vert_count-1)/step*step;
	cout<<"rem ="<<rem<<endl;

	p2Edge = new Edge* [N*N];

	//step 1: count size of each partition, and allocate memory of them
	p2EdgeCount  = new index_t[N*N];
	for(vertex_t i=0; i<N*N; i++){
		p2EdgeCount[i]=0;
	}

	vertex_t u,v;
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(index_t i=0; i<upperEdgeCount; i++){
//__sync_add_and_fetch(&k,1);
		u = OrientedEdge[i].A;
		v = OrientedEdge[i].B;
		int j = u/step;
		int k = v/step;
		int l = j*N+k;//it is the partition number to write in
		__sync_fetch_and_add(&p2EdgeCount[l],1);
//cout<<OrientedEdge[i].A<<"	"<<OrientedEdge[i].B<<" : "<<j<<" "<<k<<endl;
//		p2EdgeCount[j*N+k]++;

	}
	for(int i = 0; i<N*N; i++){
		p2Edge[i] = new Edge[p2EdgeCount[i]];
	}

//	verify
/*
	cout<<"Oriented Edges: "<<endl;
	for(index_t i = 0; i<upperEdgeCount; i++){
		cout<<OrientedEdge[i].A<<"	"<<OrientedEdge[i].B<<endl;
	}
	for(int j = 0; j<N; j++){
		for(int k = 0; k<N; k++){
			cout<<"edge set "<<j<<"-"<<k<<" with count "<<p2EdgeCount[j*N+k]<<endl;
		}
	}
	for(int j = 0; j<N*N; j++){
		cout<<"edge set count "<<p2EdgeCount[j]<<endl;
	}
*/

	//step 2: write edges
	index_t* tempEdgeCount  = new index_t[N*N];
	for(vertex_t i=0; i<N*N; i++){
		tempEdgeCount[i]=0;
	}

	for(index_t i=0; i<upperEdgeCount; i++){
		u = OrientedEdge[i].A;
		v = OrientedEdge[i].B;
		int j = u/step;
		int k = v/step;
		int l = j*N+k;//it is the partition number to write in
		index_t n = __sync_fetch_and_add(&tempEdgeCount[l],1);

		p2Edge[l][n].A = u;
		p2Edge[l][n].B = v;
//cout<<"l : "<<l<<" n : "<<n<<endl;
//cout<<p2Edge[l][n].A<<"	"<<p2Edge[l][n].B<<endl;

	}

/*	
//	verify
	cout<<"Oriented Edges: "<<endl;
	for(int j = 0; j<N*N; j++){
		cout<<"edge set count "<<p2EdgeCount[j]<<endl;
		for(index_t k = 0; k<p2EdgeCount[j]; k++){
			cout<<p2Edge[j][k].A<<"	"<<p2Edge[j][k].B<<endl;
		}
	}

*/
	
	//step 3: partition the 1-D partitioned CSR
	//
	//
	metadata = new vertex_t[2+PART_NUM*N];
	metadata[0] = N;
	metadata[1] = PART_NUM;


	p2CSR 	= new index_t* [N];
	p2Adj 	= new vertex_t** [N];
	p2Begin	= new index_t** [N];
	for(int i=0; i<N; i++){
		p2CSR[i] = new index_t [PART_NUM];
		p2Adj[i] = new vertex_t* [PART_NUM];
		p2Begin[i] = new index_t* [PART_NUM];
		for(int j=0; j<PART_NUM; j++){

			//specify the starting vertex of this 2-D begin position
			metadata[2 + i * PART_NUM + j] = i*step;	

			if(i!=N-1){
				p2CSR[i][j]=partBegin[j][step*(i+1)]-partBegin[j][step*i];
				p2Begin[i][j]= new index_t [step+1];
			}
			else{
				p2CSR[i][j]=partBegin[j][step*i+rem]-partBegin[j][step*i];
				p2Begin[i][j]= new index_t [rem+1];

			}
			p2Adj[i][j] = &partAdj[j][partBegin[j][step*i]];
			
		}
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<PART_NUM; j++){
			index_t offset = partBegin[j][step*i];
			if(i!=N-1){
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
				for(index_t k=0; k<step+1; k++){
					p2Begin[i][j][k] = partBegin[j][step*i+k]-offset;
				}
			}
			else{
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
				for(index_t k=0; k<rem+1; k++){
					p2Begin[i][j][k] = partBegin[j][step*i+k]-offset;
				}

			}
			
		}
	}

}





void graph::write_back(){

	int N = _2DPART;
	string s_meta = "metadata";		

	char* meta_file = const_cast<char*>(s_meta.c_str());

	FILE * m_File;
	m_File = fopen (meta_file, "wb");
	fwrite (metadata , sizeof(vertex_t), PART_NUM*N+2 , m_File);
	fclose (m_File);

	//edge list partitions
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){

			stringstream ss;
			ss << i;
			ss << "-";
			ss << j;
			string str = ss.str();
			
			string s_edge = "edge"+str;		

			char* edge_file = const_cast<char*>(s_edge.c_str());

		
			cout<<edge_file<<endl;
			
			
			int t = N*i+j;

			FILE * a_File;
			a_File = fopen (edge_file, "wb");
			fwrite (p2Edge[t] , sizeof(Edge), p2EdgeCount[t], a_File);
			fclose (a_File);


		}
	}

	//csr partitions, CSR-i-j follows partition-of-row-column
	for(int i=0;i<N;i++){
		for(int j=0;j<PART_NUM;j++){

			stringstream ss;
			ss << i;
			ss << "-";
			ss << j;
			string str = ss.str();
			
			string s_begin = "begin"+str;
			string s_adj = "adjacent"+str;

			char* begin_file = const_cast<char*>(s_begin.c_str());
			char* adj_file = const_cast<char*>(s_adj.c_str());
		
			cout<<begin_file<<endl;
//			cout<<adj_file<<endl;
			
			
//			int t = N*i+j;

			FILE * a_File;
			a_File = fopen (adj_file, "wb");
			fwrite (p2Adj[i][j] , sizeof(vertex_t), p2CSR[i][j], a_File);
			fclose (a_File);

			if(i==N-1){
				FILE * b_File;
				b_File = fopen (begin_file, "wb");
				fwrite (p2Begin[i][j] , sizeof(index_t), rem+1, b_File);
				fclose (b_File);
			}		
			else{
				FILE * b_File;
				b_File = fopen (begin_file, "wb");
				fwrite (p2Begin[i][j] , sizeof(index_t), step+1, b_File);
				fclose (b_File);
			}	

		}
	}
	//verify
	/*
	for(int i=0; i<N; i++){
		for(vertex_t j = 0; j<PART_NUM; j++){
cout<<"begin "<<i<<"-"<<j<<endl;
			if(i!=N-1){
				for(index_t k=0; k<step+1; k++){
					cout<<p2Begin[i][j][k]<<" ";
				}
			}
			else{
				for(index_t k=0; k<rem+1; k++){
					cout<<p2Begin[i][j][k]<<" ";
				}

			}
			cout<<endl;
		}
		cout<<endl;

	}
	for(int i=0;i<PART_NUM*N+2;i++){
		cout<<metadata[i]<<" ";
	}
	cout<<endl;
	*/

}
