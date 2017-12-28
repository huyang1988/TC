//Graph format: 
//Simplified json format: 
//src degree dest0 dest1 ...


//#include "graph.h"
#include "comm.h"
#include <fstream>
#include <omp.h>

#define FILE_NOT_EXIST	1
#define FILE_EXIST	0
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
	count = (index_t *)malloc(256*256*sizeof(index_t));
//	valid = (int *)malloc(vert_count*sizeof(int));
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
	
	for(int i=0;i<EDGE_PART;i++){
		for(int j=0;j<EDGE_PART;j++){

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
			
			
			int t = 4*i+j;

			FILE * a_File;
			a_File = fopen (adj_file, "wb");
			fwrite (partAdj[t] , sizeof(vertex_t), partEdgeCount[t], a_File);
			fclose (a_File);

			if(i==EDGE_PART-1){
				FILE * b_File;
				b_File = fopen (begin_file, "wb");
				fwrite (partBegin[t] , sizeof(index_t), rem+1, b_File);
				fclose (b_File);
			}		
			else{
				FILE * b_File;
				b_File = fopen (begin_file, "wb");
				fwrite (partBegin[t] , sizeof(index_t), step+1, b_File);
				fclose (b_File);
			}	

		}
	}

/*
	cout<<"size of vertex_t "<<sizeof(vertex_t)<<endl;
	cout<<"size of index_t "<<sizeof(index_t)<<endl;
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
*/
}

