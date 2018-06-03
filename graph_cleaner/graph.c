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

void graph::validation(){
	index_t mycount=0;
//cout<<"validation "<<upperEdgeCount<<endl;
	for(index_t i=0; i<upperEdgeCount; i++){
		vertex_t U=upperHead[i];
		vertex_t V=upperAdj[i];
		index_t m=upperBegin[U+1]-upperBegin[U];
		index_t n=upperBegin[V+1]-upperBegin[V];

//cout<<"edge: "<<i<<" "<<U<<"-"<<V<<" ";
//cout<<"degree: "<<m<<" "<<n<<endl;

		vertex_t *u = &upperAdj[upperBegin[U]];
		vertex_t *v = &upperAdj[upperBegin[V]];

		vertex_t u1=0;
		vertex_t v1=0;
		while(u1<m && v1<n){
			vertex_t x=u[u1];
			vertex_t y=v[v1];
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
	cout<<"validation version tc = "<<mycount<<endl;
}

void graph::bsvalidation(){
	index_t mycount=0;
	for(index_t i=0; i<upperEdgeCount; i++){
		vertex_t A=upperHead[i];
		vertex_t B=upperAdj[i];
		index_t m=upperBegin[A+1]-upperBegin[A];
		index_t n=upperBegin[B+1]-upperBegin[B];
	
		if(m<n){
			A=A^B;
			B=A^B;
			A=A^B;
			m=m^n;
			n=m^n;
			m=m^n;
		}


		vertex_t *a = &upperAdj[upperBegin[A]];
		vertex_t *b = &upperAdj[upperBegin[B]];
//		cout<<b[0];

		index_t bot;
		index_t top;
		index_t r;

		for(size_t j=0; j<n; j++){
			vertex_t X = b[j];
			vertex_t Y;
			bot = 0;
			top = m-1;
			
			while(top>=bot){
				r = (top+bot)/2;
				Y = a[r];
				if(X==Y){
					mycount++;
					break;
				}
				else if(X<Y){
					top = r-1;
				}
				else if(X>Y){
					bot = r+1;
				}
			}
		}

	}
	cout<<"validation version tc = "<<mycount<<endl;
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

}

