//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...


#include "graph.h"
#include "comm.h"
#include <fstream>
#include <omp.h>

using namespace std;

graph::graph(
	string jsonfile)//,
{
	cout<<"read from folder "<<jsonfile<<endl;
	
	string s_begin = jsonfile+"/begin.bin";
	string s_adj = jsonfile+"/adjacent.bin";

	char* begin_file = const_cast<char*>(s_begin.c_str());
	char* adj_file = const_cast<char*>(s_adj.c_str());

	vert_count = fsize(begin_file)/sizeof(index_t) - 1;
	edge_count = fsize(adj_file)/sizeof(vertex_t);

	cout<<"vert:"<< vert_count<<"  edge: "<<edge_count<<endl;


	FILE *pFile= fopen(adj_file,"rb");
	adj_list = (vertex_t *)malloc(fsize(adj_file));
	fread(adj_list,sizeof(vertex_t),edge_count,pFile);
	fclose(pFile);



	FILE *pFile3 = fopen(begin_file,"rb");
	beg_pos = (index_t *)malloc(fsize(begin_file));
	fread(beg_pos,sizeof(index_t),vert_count+1,pFile3);
	fclose(pFile3);
}


void graph::process_graph()
{

	/*
	//Write to a file
	string filenames = "output.dat";
	//sprintf(filenames,"rmat_%lu_%lu.dat",log_numverts, edge_factor);
	std::ofstream outfile1(filenames.c_str());
	outfile1 <<"AdjacencyGraph" << endl;
	outfile1 <<vert_count << endl;
	outfile1 << upperBegin[vert_count] << endl;
	for(index_t i=0; i<vert_count; i++)
	{
		outfile1<< upperBegin[i] <<endl;
	}

	for(index_t i=0; i<vert_count; i++)
	{
		for(index_t j=upperBegin[i]; j<upperBegin[i+1]; j++)
		{   
			outfile1<<upperAdj[j]<< endl;
		}   
	}
	outfile1.close();
	*/

//	cout<<"size of vertex_t "<<sizeof(vertex_t)<<endl;
//	cout<<"size of index_t "<<sizeof(index_t)<<endl;
	FILE * a_File;
	a_File = fopen ("adjacent.bin", "wb");
	fwrite (upperAdj , sizeof(vertex_t), upperEdgeCount, a_File);
	fclose (a_File);
	
	FILE * h_File;
	h_File = fopen ("head.bin", "wb");
	fwrite (upperHead , sizeof(vertex_t), upperEdgeCount, h_File);
	fclose (h_File);

	FILE * b_File;
	b_File = fopen ("begin.bin", "wb");
	fwrite (upperBegin , sizeof(index_t), vert_count+1, b_File);
	fclose (b_File);

}

void graph::rank_by_degree(){
	upperBegin	= new index_t[vert_count+1];
	upperBegin[0]=0;
//#pragma omp parallel for num_threads(56) schedule(static)
#pragma omp parallel for num_threads(56) schedule(dynamic,1024)
	for(vertex_t i=0; i<vert_count; i++){
//		upperBegin[i+1]=upperBegin[i];//upperDegree[i]=0;
		upperBegin[i+1]=0;
		index_t dh=beg_pos[i+1]-beg_pos[i];
		index_t j=beg_pos[i];

		while(j<beg_pos[i+1]){
			vertex_t a=adj_list[j];
			index_t da=beg_pos[a+1]-beg_pos[a];
			if(dh<da || (dh==da && i<a)){
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
		index_t dh=beg_pos[i+1]-beg_pos[i];
		index_t jj=upperBegin[i];
		while(j<beg_pos[i+1]){
			vertex_t a=adj_list[j];
			index_t da=beg_pos[a+1]-beg_pos[a];
			if(dh<da || (dh==da && i<a)){
				upperAdj[jj] =adj_list[j];
				upperHead[jj]=i;
				jj++;//k++;
			}
			j++;
		}
	}
	
	cout<<"upper Edge Count= "<<upperEdgeCount<<"\n";
}
