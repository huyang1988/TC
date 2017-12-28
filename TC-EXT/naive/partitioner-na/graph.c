#include "graph.h"
#include "comm.h"
#include "wtime.h"
#include <fstream>
#include <omp.h>
#include <math.h> 
using namespace std;

//function BS is only used in further partition, it searches the position satisfy that A[n]+2n < x.
vertex_t BS(index_t x, index_t*A, vertex_t bot, vertex_t top){

//	for(int i=bot;i<=top;i++){
//		cout<<A[i]<<" ";
//	}
//	cout<<"\n";

	vertex_t r= (bot+top)/2;
//	int result;
	while(top>bot){

		index_t y = A[r] + r*2;
		if(x<y){
			top = r;
		}
		else if(x>=y){
			bot = r+1;
		}
		else if(x==y){
			while(x==y && r<=top){
				r++;
			}
			r--;
			break;
		}
		r = (bot+top)/2;
	}
	return r;
}

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

graph::graph(
	string jsonfile)//,
{




	cout<<"read from folder "<<jsonfile<<endl;
	
	string s_begin = jsonfile+"/begin.bin";
	string s_adj = jsonfile+"/adjacent.bin";
	string s_head = jsonfile+"/head.bin";
	string s_degree = jsonfile+"/degree.bin";

	char* begin_file = const_cast<char*>(s_begin.c_str());
	char* adj_file = const_cast<char*>(s_adj.c_str());

	vert_count = fsize(begin_file)/sizeof(index_t) - 1;
	edge_count = fsize(adj_file)/sizeof(vertex_t);

	cout<<"vert:"<< vert_count<<"  edge: "<<edge_count<<endl;
	
	PART_NUM = 1 + (int)sqrt((edge_count-1)/PartitionSize/2);
//	PART_NUM = (int)sqrt((graph_d->edge_count-1)/PartitionSize);
	cout<<"PART NUMBER = "<<PART_NUM<<endl;

//	edge_limit = PartitionSize - sizeof(index_t)/sizeof(vertex_t)*vert_count;
//	_2DPART   = (edge_count/2/PART_NUM + vert_count * 2-1)/PartitionSize + 1;


//	cout<<"edge limit = "<<edge_limit<<endl;
//if(edge_limit<=1000){
//	cout<<"Too few edge can be put in one CSR partition!"<<endl;
//	exit(1);
//}

//	PART_NUM = (edge_count/2 -1 )/edge_limit+1;
//	PART_NUM = (edge_count/2 -1 )/PartitionSize+1;
	
	cout<<"the partition number = "<<PART_NUM<<endl;

	FILE *pFile= fopen(adj_file,"rb");
	adj_list = (vertex_t *)malloc(fsize(adj_file));
	fread(adj_list,sizeof(vertex_t),edge_count,pFile);
	fclose(pFile);


	FILE *pFile3 = fopen(begin_file,"rb");
	beg_pos = (index_t *)malloc(fsize(begin_file));
	fread(beg_pos,sizeof(index_t),vert_count+1,pFile3);
	fclose(pFile3);
/*	
	//step 2: write metadata
	int* metadata = new int[2];
	metadata[0] = 1;
	metadata[1] = PART_NUM;
	string s_meta = "metadata";		

	char* meta_file = const_cast<char*>(s_meta.c_str());

	FILE * m_File;
	m_File = fopen (meta_file, "wb");
	fwrite (metadata , sizeof(int), 2 , m_File);
	fclose (m_File);
*/
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
			vertex_t h=i;
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
	vertex_t cutpoint[PART_NUM+1];	// the colum value begin from partition i
	cutpoint[0] = 0;
	cutpoint[PART_NUM] = vert_count+1;

//#pragma omp parallel for num_threads(PART_NUM) schedule(static)
	for(int i=1; i<PART_NUM; i++){
		index_t K=i*upperEdgeCount/PART_NUM;
cout<<"K="<<K<<"\n";
		vertex_t index = BinarySearch(K, inBegin, 0, vert_count-1); // binary search
		cutpoint[i] = index; 			//used by each neigbhor list to find a place to cut	
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
		vertex_t h=i;
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


void graph::further_partition(){
	
	//Partition each 1-D column into balanced parts
	//build new begin position arrays for each
	
	_2DPART = (upperEdgeCount/PART_NUM - 1) / PartitionSize + 1;
	cout<<"row-wise partition number = "<<_2DPART<<endl;

	index_t temp_meta[PART_NUM*_2DPART];//record the vertex id offset of each partition
	index_t temp_vert[PART_NUM*_2DPART];//record the vertex count of each partition, index k = (_2DPART*i+j) for partition(j,i)
		
	//step 2: binary search in-degree for partition
	for(int i=0; i<PART_NUM; i++){
		temp_meta[i*_2DPART] = 0;
		for(int j=1; j<_2DPART; j++){
			int k = _2DPART * i + j;
			index_t part_total = partEdgeCount[i]; 
			index_t K = j*part_total/_2DPART;
			vertex_t index = BinarySearch(K, partBegin[i], 0, vert_count-1); // binary search
			temp_meta[k] = index;
			temp_vert[k-1] = temp_meta[k]-temp_meta[k-1];	


		}
		temp_vert[(i+1)*_2DPART-1] = vert_count - temp_meta[(i+1)*_2DPART-1];
	}
	
	//step 2: write metadata
	int N = _2DPART;
	metadata = new vertex_t[2+PART_NUM*_2DPART];
	metadata[0] = _2DPART;
	metadata[1] = PART_NUM;

	p2CSR 	= new index_t* [N];
	p2Adj 	= new vertex_t** [N];
	p2Begin	= new index_t** [N];
	
	for(int i=0; i<N*PART_NUM; i++){		
		metadata[i+2] = temp_meta[i];
	}


	//step 3: partition CSR
	for(int i=0; i<N*PART_NUM+2; i++){
		cout<<metadata[i]<<" ";
	}
	cout<<endl;
	
	for(int i=0; i<N; i++){
		p2CSR[i] = new index_t [PART_NUM];
		p2Adj[i] = new vertex_t* [PART_NUM];
		p2Begin[i] = new index_t* [PART_NUM];
		for(int j=0; j<PART_NUM; j++){
			int t = j*_2DPART+i;

			p2CSR[i][j] = partBegin[j][temp_meta[t]+temp_vert[t]]-partBegin[j][temp_meta[t]];
			p2Adj[i][j] = &partAdj[j][partBegin[j][temp_meta[t]]];
			p2Begin[i][j]= new index_t [temp_vert[t]+1];

cout<<"edge number of partition "<<i<<"-"<<j<<" : "<<p2CSR[i][j]<<endl;
			
		}
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<PART_NUM; j++){
			int t = j*_2DPART+i;

			index_t offset = partBegin[j][temp_meta[t]];

//#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
			for(vertex_t k=0; k<temp_vert[t]+1; k++){
				p2Begin[i][j][k] = partBegin[j][temp_meta[t]+k]-offset;
			}

			
		}
	}
	//csr partitions, CSR-i-j follows partition-of-row-column
	string s_meta = "metadata";		

	char* meta_file = const_cast<char*>(s_meta.c_str());

	FILE * m_File;
	m_File = fopen (meta_file, "wb");
	fwrite (metadata , sizeof(vertex_t), N*PART_NUM+2 , m_File);
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
			fwrite (p2Begin[i][j] , sizeof(index_t), temp_vert[i]+1, b_File);
			fclose (b_File);

		}
	}
	
	//verify
	/*
	cout<<"meta ";
	for(int i=0;i<N*PART_NUM;i++){
		cout<<temp_meta[i]<<" ";
	}
	cout<<"vert ";
	for(int i=0;i<N*PART_NUM;i++){
		cout<<temp_vert[i]<<" ";
	}
	cout<<endl;
	for(int i=0; i<N; i++){
		for(vertex_t j = 0; j<PART_NUM; j++){
			int t = j*_2DPART+i;
			cout<<"begin "<<i<<"-"<<j<<endl;
			for(index_t k=0; k<temp_vert[t]+1; k++){
				cout<<p2Begin[i][j][k]<<" ";
			}
			cout<<endl;
			cout<<"adj "<<i<<"-"<<j<<endl;
			for(index_t k=0; k<p2CSR[i][j]; k++){
				cout<<p2Adj[i][j][k]<<" ";
			}
			cout<<endl;
		}
		cout<<endl;

	}
	*/


}






void graph::edge_2d(){
cout<<"PARTITION EDGE LIST"<<endl;
	vertex_t N = _2DPART;

	cout<<"N "<<N<<endl;

	p2Edge = new Edge** [PART_NUM];
	p2EdgeCount = new index_t* [PART_NUM];

	for(int k=0; k<PART_NUM; k++){
cout<<"	PARTITION EDGE LIST round"<<k<<endl;
		
		cout<<"metadata ";
		for(int i=0;i<N*PART_NUM+2;i++){
			cout<<metadata[i]<<" ";
		}
		cout<<endl;	

		p2Edge[k] = new Edge* [N*N];

		//step 1: count size of each partition, and allocate memory of them
		p2EdgeCount[k]  = new index_t[N*N];
		for(vertex_t i=0; i<N*N; i++){
			p2EdgeCount[k][i]=0;
		}

		vertex_t u,v;
	//#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
		for(index_t i=0; i<upperEdgeCount; i++){
			u = OrientedEdge[i].A;
			v = OrientedEdge[i].B;


			int m = 0;
			while(u>=metadata[3+m+_2DPART*k] && m<N-1){
				m++;
			}
			int n = 0;
			while(v>=metadata[3+n+_2DPART*k] && n<N-1){
				n++;
			}
	//		int j = u/step;
	//		int k = v/step;
			int l = m*N+n;//it is the partition number to write in



	//		__sync_fetch_and_add(&p2EdgeCount[l],1);
	//cout<<OrientedEdge[i].A<<"	"<<OrientedEdge[i].B<<" : "<<j<<" "<<k<<endl;
			p2EdgeCount[k][l]++;

		}

	//allocate space
		for(int i = 0; i<N*N; i++){
			p2Edge[k][i] = new Edge[p2EdgeCount[k][i]];
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
			int m = 0;
			while(u>=metadata[3+m+_2DPART*k] && m<N-1){
				m++;
			}
			int n = 0;
			while(v>=metadata[3+n+_2DPART*k] && n<N-1){
				n++;
			}
	//		int j = u/step;
	//		int k = v/step;
			int l = m*N+n;//it is the partition number to write in
	//		index_t n = __sync_fetch_and_add(&tempEdgeCount[l],1);
			index_t p = tempEdgeCount[l];
			tempEdgeCount[l]++;

			p2Edge[k][l][p].A = u;
			p2Edge[k][l][p].B = v;
	//cout<<OrientedEdge[i].A<<"	"<<OrientedEdge[i].B<<" : "<<j<<" "<<k<<endl;
	//cout<<p2Edge[l][n].A<<"	"<<p2Edge[l][n].B<<" : "<<j<<" "<<k<<endl;

		}

/*
 		for(int j = 0; j<N; j++){
			for(int p = 0; p<N; p++){
				int l = j*N+p;
//				cout<<"edge set "<<j<<"-"<<p<<" count "<<p2EdgeCount[k][l]<<endl;
				cout<<"edge set "<<j<<"-"<<p<<" count "<<p2EdgeCount[k][l]<<endl;
	//			for(index_t p = 0; p<5; p++){
				for(index_t p = 0; p<p2EdgeCount[k][l]; p++){
					cout<<p2Edge[k][l][p].A<<"	"<<p2Edge[k][l][p].B<<endl;
				}
			}
		}
*/

		for(int i=0;i<N;i++){
			for(int j=0;j<N;j++){

				stringstream ss;
				ss << i;
				ss << "-";
				ss << j;
				ss<<"-";
				ss<<k;
				string str = ss.str();
				
				string s_edge = "edge"+str;		

				char* edge_file = const_cast<char*>(s_edge.c_str());

			
	//			cout<<edge_file<<endl;
				
				
				int t = N*i+j;

				FILE * a_File;
				a_File = fopen (edge_file, "wb");
				fwrite (p2Edge[k][t] , sizeof(Edge), p2EdgeCount[k][t], a_File);
				fclose (a_File);


			}
		}
	}
	

}



