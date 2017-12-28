//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...


//#include "graph.h"
#include "comm.h"
#include "wtime.h"
#include <fstream>
#include <omp.h>

using namespace std;

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
	char* head_file = const_cast<char*>(s_head.c_str());
	char* degree_file = const_cast<char*>(s_degree.c_str());

	vert_count = fsize(begin_file)/sizeof(index_t) - 1;
	edge_count = fsize(head_file)/sizeof(vertex_t);

	cout<<"vert:"<< vert_count<<"  edge: "<<edge_count<<endl;

	PART_NUM = (edge_count/2 -1 )/PartitionSize+1;
	
	cout<<"the partition number = "<<PART_NUM<<endl;

	FILE *pFile= fopen(adj_file,"rb");
	adj_list = (vertex_t *)malloc(fsize(adj_file));
	fread(adj_list,sizeof(vertex_t),edge_count,pFile);
	fclose(pFile);

	FILE *pFile1= fopen(head_file,"rb");
	head_list = (vertex_t *)malloc(fsize(head_file));
	fread(head_list,sizeof(vertex_t),edge_count,pFile1);
	fclose(pFile1);

//	FILE *pFile2= fopen(degree_file,"rb");
//	adj_card = (index_t *)malloc(fsize(degree_file));
//	fread(adj_card,sizeof(index_t),vert_count,pFile2);
//	fclose(pFile2);

	FILE *pFile3 = fopen(begin_file,"rb");
	beg_pos = (index_t *)malloc(fsize(begin_file));
	fread(beg_pos,sizeof(index_t),vert_count+1,pFile3);
	fclose(pFile3);

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

	//cout<<"edge: "<<i<<" "<<U<<"-"<<V<<" ";
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
				}
			}
		}
	
		global_count += mycount;
	}
	cout<<"merge version tc = "<<global_count<<endl;
}
