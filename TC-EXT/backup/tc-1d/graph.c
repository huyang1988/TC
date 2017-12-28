//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...


#include "graph.h"
#include "comm.h"
#include <fstream>
#include <omp.h>

#include "graph.h"
#define FILE_NOT_EXIST	1
#define FILE_EXIST	0
#define P	32
using namespace std;

graph::graph(
	string jsonfile)//,
{
	cout<<"read from folder "<<jsonfile<<endl;
	partEdgeCount 	= new index_t [PART_NUM];
	partAdj   	= new vertex_t * [PART_NUM];
	partBegin 	= new index_t * [PART_NUM];

	for(int i = 0; i<PART_NUM; i++){

		stringstream ss;
		ss<<i;
		string str = ss.str();

		string s_begin = jsonfile+"/begin-"+str;
		string s_adj = jsonfile+"/adjacent-"+str;

		char* begin_file = const_cast<char*>(s_begin.c_str());
		char* adj_file = const_cast<char*>(s_adj.c_str());

		vert_count 	 = fsize(begin_file)/sizeof(index_t)-1;	
cout<<"vert count = "<<vert_count<<endl;	
		partEdgeCount[i] = fsize(adj_file)/sizeof(vertex_t);
cout<<"part "<<i<<" edge count = "<<partEdgeCount[i]<<endl;
		
		
		FILE *pFile = fopen(adj_file,"rb");
		partAdj[i] = new vertex_t [partEdgeCount[i]];
		fread(partAdj[i], sizeof(vertex_t), partEdgeCount[i], pFile);
		fclose(pFile);

		FILE *pFile3 = fopen(begin_file,"rb");
		partBegin[i] = new index_t [vert_count+1];
		fread(partBegin[i], sizeof(index_t), vert_count+1 , pFile3);
		fclose(pFile3);
		
		


	}
	
	
	
	string s_edge = jsonfile+"/edge";
	char* edge_file = const_cast<char*>(s_edge.c_str());
	
		
	upperEdgeCount = fsize(edge_file)/sizeof(Edge);
cout<<"upperEdgeCount = "<<upperEdgeCount<<endl;
	OrientedEdge = new Edge[upperEdgeCount];
	FILE *pFile1 = fopen(edge_file,"rb");
	fread(OrientedEdge, sizeof(Edge), upperEdgeCount, pFile1);

	
	

	count = (index_t *)malloc(256*256*sizeof(index_t));
}




void graph::part_validation(){
//cout<<"validation "<<upperEdgeCount<<endl;
//
	long int global_count=0;
	for(int p=0; p<PART_NUM; p++){
	
		index_t mycount=0;
#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(dynamic,1024)
		for(index_t i=0; i<upperEdgeCount; i++){
			//get the edge A-B from global graph
			vertex_t A=OrientedEdge[i].A;
			vertex_t B=OrientedEdge[i].B;

			//degree of current partition
			index_t m=partBegin[p][A+1]-partBegin[p][A];
			index_t n=partBegin[p][B+1]-partBegin[p][B];

//cout<<"edge: "<<A<<"-"<<B<<endl;
	//cout<<"degree: "<<m<<" "<<n<<endl;

			vertex_t *a = &partAdj[p][partBegin[p][A]];
			vertex_t *b = &partAdj[p][partBegin[p][B]];

			vertex_t u1=0;
			vertex_t v1=0;
			while(u1<m && v1<n){
				vertex_t x=a[u1];
				vertex_t y=b[v1];
				if(x<y){
					u1++;
				}
				else if(x>y){
					v1++;
				}
				else if(x==y){
					u1++;
					v1++;
					mycount++;
//cout<<"find triangle "<<A<<"-"<<B<<"-"<<x<<endl;
				}
			}
		}
	
		global_count += mycount;
	}
	cout<<"merge version tc = "<<global_count<<endl;
}
