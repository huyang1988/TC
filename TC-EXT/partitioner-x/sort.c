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
//
//	/*
//verify
	for(int i=0; i<PART_NUM; i++){
		cout<<"part "<<i<<endl;
		cout<<"part size "<<partEdgeCount[i]<<endl;
		for(vertex_t j = 0; j<5; j++){
//		for(vertex_t j = 0; j<vert_count+1; j++){
			cout<<partBegin[i][j]<<" ";
		}
		cout<<endl;
		for(index_t j = 0; j<10; j++){
//		for(index_t j = 0; j<partEdgeCount[i]; j++){
			cout<<partAdj[i][j]<<" ";
		}
		cout<<endl;

	}
//	*/
	
}


void graph::further_partition(){
	//step 1: partition the 1-D partitioned CSR iteratively
	//
	//
	int N=0;	//horizontal partition number, accumulate on the fly
	
	long limit = 1+(upperEdgeCount+1)/PartitionSize;

cout<<"limit partition number "<<limit<<endl;

//	vertex_t* temp_meta = new vertex_t[limit+5];
	vertex_t* temp_meta = new vertex_t[1000];
	temp_meta[0]=0;

	
	
//	int i = 2;//search on vertical partition i
	vertex_t last[PART_NUM];
	for(int i=0; i<PART_NUM; i++)
		last[i] = 0;
	index_t	 last_size[PART_NUM];
	for(int i=0; i<PART_NUM; i++)
		last_size[i] = 0;
	vertex_t c,cc;
	index_t rem=upperEdgeCount;

	while(rem>0){
		cc = vert_count;
		for(int i=0; i<PART_NUM; i++){
//			c =  BinarySearch(PartitionSize+last_size[i], partBegin[i], last[i], vert_count);
			c =  BS(PartitionSize+last_size[i], partBegin[i], last[i], vert_count);
//			c =  BS(PartitionSize+last_size[i], partBegin[i], 0, vert_count);
cout<<"round "<<N<<" cut vert at "<<c<<" with size "<<partBegin[i][c]<<endl;
			if(c<cc){
				cc=c;
			}
		}

		vertex_t k = cc-1;
		if(cc==vert_count){
			k = cc;
		}
		temp_meta[N+1] = k;

		for(int i=0; i<PART_NUM; i++){
			cout<<"round "<<N<<" cut vert at "<<k<<" with size "<<partBegin[i][k]<<endl;
		}
		
		N++;

		rem = 0;
		for(int i=0; i<PART_NUM; i++){
			last[i] = cc;
			last_size[i] = partBegin[i][cc];
			rem += partEdgeCount[i] - last_size[i];
		}
cout<<"rem "<<rem<<" N "<<N<<endl;

	}

	//result
	cout<<"total round "<<N<<endl;









//	vertex_t 	
	//step 2: write metadata
	metadata = new vertex_t[2+N+1];
	metadata[0] = N;
	metadata[1] = PART_NUM;

	p2CSR 	= new index_t* [N];
	p2Adj 	= new vertex_t** [N];
	p2Begin	= new index_t** [N];
	
	for(int i=0; i<N+1; i++){		
		metadata[i+2] = temp_meta[i];
	}


	//step 3: partition CSR
//	for(int i=0; i<N+2; i++){
//		cout<<metadata[i]<<" ";
//	}
	cout<<endl;
	
	for(int i=0; i<N; i++){
		p2CSR[i] = new index_t [PART_NUM];
		p2Adj[i] = new vertex_t* [PART_NUM];
		p2Begin[i] = new index_t* [PART_NUM];
		for(int j=0; j<PART_NUM; j++){

			p2CSR[i][j]=partBegin[j][temp_meta[i+1]]-partBegin[j][temp_meta[i]];
			p2Adj[i][j] = &partAdj[j][partBegin[j][temp_meta[i]]];
			p2Begin[i][j]= new index_t [temp_meta[i+1]-temp_meta[i]+1];

//cout<<"edge number of partition "<<i<<"-"<<j<<" : "<<p2CSR[i][j]<<endl;
			
		}
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<PART_NUM; j++){
			index_t offset = partBegin[j][temp_meta[i]];
//#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
			for(vertex_t k=0; k<temp_meta[i+1]-temp_meta[i]+1; k++){
				p2Begin[i][j][k] = partBegin[j][temp_meta[i]+k]-offset;
			}

			
		}
	}
	//verify
	/*
	cout<<"metadata ";
	for(int i=0;i<N+2;i++){
		cout<<metadata[i]<<" ";
	}
	cout<<endl;
	for(int i=0; i<N; i++){
		for(vertex_t j = 0; j<PART_NUM; j++){
			cout<<"begin "<<i<<"-"<<j<<endl;
			for(index_t k=0; k<temp_meta[i+1]-temp_meta[i]+1; k++){
				cout<<p2Begin[i][j][k]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;

	}
	*/
	//csr partitions, CSR-i-j follows partition-of-row-column
	string s_meta = "metadata";		

	char* meta_file = const_cast<char*>(s_meta.c_str());

	FILE * m_File;
	m_File = fopen (meta_file, "wb");
	fwrite (metadata , sizeof(vertex_t), N+2 , m_File);
	fclose (m_File);
	
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
		
//			cout<<begin_file<<endl;
//			cout<<adj_file<<endl;
			
			
//			int t = N*i+j;

			FILE * a_File;
			a_File = fopen (adj_file, "wb");
			fwrite (p2Adj[i][j] , sizeof(vertex_t), p2CSR[i][j], a_File);
			fclose (a_File);

			FILE * b_File;
			b_File = fopen (begin_file, "wb");
			fwrite (p2Begin[i][j] , sizeof(index_t), temp_meta[i+1]-temp_meta[i]+1, b_File);
			fclose (b_File);

		}
	}
	
}



void graph::edge_2d(){

	vertex_t N = metadata[0];

	cout<<"N "<<N<<endl;
	cout<<"metadata ";
	for(int i=0;i<N+3;i++){
		cout<<metadata[i]<<" ";
	}
	cout<<endl;	

	p2Edge = new Edge* [N*N];

	//step 1: count size of each partition, and allocate memory of them
	p2EdgeCount  = new index_t[N*N];
	for(vertex_t i=0; i<N*N; i++){
		p2EdgeCount[i]=0;
	}

	vertex_t u,v;
//#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(index_t i=0; i<upperEdgeCount; i++){
		u = OrientedEdge[i].A;
		v = OrientedEdge[i].B;


		int j = 0;
		while(u>=metadata[3+j]){
			j++;
		}
		int k = 0;
		while(v>=metadata[3+k]){
			k++;
		}
//		int j = u/step;
//		int k = v/step;
		int l = j*N+k;//it is the partition number to write in



//		__sync_fetch_and_add(&p2EdgeCount[l],1);
//cout<<OrientedEdge[i].A<<"	"<<OrientedEdge[i].B<<" : "<<j<<" "<<k<<endl;
		p2EdgeCount[l]++;

	}

//allocate space
	for(int i = 0; i<N*N; i++){
		p2Edge[i] = new Edge[p2EdgeCount[i]];
	}


	//step 2: write edges
	index_t* tempEdgeCount  = new index_t[N*N];
	for(vertex_t i=0; i<N*N; i++){
		tempEdgeCount[i]=0;
	}

//#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(index_t i=0; i<upperEdgeCount; i++){
		u = OrientedEdge[i].A;
		v = OrientedEdge[i].B;
		int j = 0;
		while(u>=metadata[3+j]){
			j++;
		}
		int k = 0;
		while(v>=metadata[3+k]){
			k++;
		}
//		int j = u/step;
//		int k = v/step;
		int l = j*N+k;//it is the partition number to write in
//		index_t n = __sync_fetch_and_add(&tempEdgeCount[l],1);
		index_t n = tempEdgeCount[l];
		tempEdgeCount[l]++;

		p2Edge[l][n].A = u;
		p2Edge[l][n].B = v;
//cout<<OrientedEdge[i].A<<"	"<<OrientedEdge[i].B<<" : "<<j<<" "<<k<<endl;
//cout<<p2Edge[l][n].A<<"	"<<p2Edge[l][n].B<<" : "<<j<<" "<<k<<endl;

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
			int l = j*N+k;
			for(index_t t = 0; t<p2EdgeCount[l]; t++){
				cout<<p2Edge[l][t].A<<"	"<<p2Edge[j][k].B<<endl;
			}
		}
	}
	for(int j = 0; j<N*N; j++){
		cout<<"edge set count "<<p2EdgeCount[j]<<endl;
		for(index_t k = 0; k<p2EdgeCount[j]; k++){
			cout<<p2Edge[j][k].A<<"	"<<p2Edge[j][k].B<<endl;
		}
	}
	*/
//	verify
	for(int j = 0; j<N; j++){
		for(int k = 0; k<N; k++){
			int l = j*N+k;
			cout<<"edge set "<<j<<"-"<<k<<" count "<<p2EdgeCount[l]<<endl;
			for(index_t k = 0; k<5; k++){
//			for(index_t k = 0; k<p2EdgeCount[l]; k++){
				cout<<p2Edge[l][k].A<<"	"<<p2Edge[l][k].B<<endl;
			}
		}
	}


	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){

			stringstream ss;
			ss << i;
			ss << "-";
			ss << j;
			string str = ss.str();
			
			string s_edge = "edge"+str;		

			char* edge_file = const_cast<char*>(s_edge.c_str());

		
//			cout<<edge_file<<endl;
			
			
			int t = N*i+j;

			FILE * a_File;
			a_File = fopen (edge_file, "wb");
			fwrite (p2Edge[t] , sizeof(Edge), p2EdgeCount[t], a_File);
			fclose (a_File);


		}
	}

	

}





