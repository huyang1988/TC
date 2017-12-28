#include "worker.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include "comm.h"
#include <fstream>
#include <omp.h>


using namespace std;

worker::worker(string s){


	time_io = 0;
	time_compute = 0;

	csr1 =new CSR(s);
	csr2 =new CSR(s);
	buffer = new Edge_chunk(s);	
	count=0;
	edge_error=0;
	
	
	string s_meta = s + "/metadata";
	char* meta_file = const_cast<char*>(s_meta.c_str());
	size_t meta_size = fsize(meta_file)/sizeof(vertex_t);
	FILE *mFile= fopen(meta_file,"rb");
	vertex_t *meta = (vertex_t *)malloc(fsize(meta_file));
	fread(meta,sizeof(vertex_t),meta_size,mFile);
	fclose(mFile);
	N1 = meta[0];
	N2 = meta[1];

	offset = &meta[2];
//	/*
	cout<<"N1 "<<N1<<" N2 "<<N2<<endl;
//	for(size_t i=0;i<meta_size;i++){
//		cout<<meta[i]<<" ";
//	}
//	cout<<endl;
	for(size_t i=0;i<meta_size-2;i++){
		cout<<offset[i]<<" ";
	}
	cout<<endl;
//	*/
}

void worker::work(){
	// outer loop start from k! very important..
	for(int k = 0; k<N2; k++){
		for(int i = 0; i<N1; i++){
			double t0 = wtime();			
			
			csr1->CSR_load(i,k);
			csr1->start = offset[i];
			csr1->end = csr1->start + csr1->beg_size -1;
			double t1 = wtime();	
			time_io  += t1-t0; 
			
			
			for(int j=0; j<N1; j++){
				double t0 = wtime();			
				csr2->CSR_load(j,k);
				csr2->start = offset[j];
				csr2->end = csr2->start + csr2->beg_size -1;
				
				double t1 = wtime();	
				time_io  += t1-t0; 
				proc(i,j,k);
				csr2->CSR_free();
			}
			csr1->CSR_free();
		}
	}
	cout<<"total number = "<<count<<endl;
	cout<<"edge error = "<<edge_error<<endl;
//	return;
}

void worker::proc(int i, int j, int k){

//	cout<<"subproblem "<<i<<"-"<<j<<"-"<<k<<endl;
//	vertex_t offsetA = offset[i];
//	vertex_t offsetB = offset[j];
	vertex_t offsetA = offset[i];
	vertex_t offsetB = offset[j];


	long sub_edge_error=0;

	buffer->Edge_open(i,j);
	for(int c = 0; c<buffer->chunk_num; c++){
		double t0 = wtime();			
		buffer->Edge_reload();
		double t1 = wtime();	
time_io  += t1-t0; 
//cout<<"chunk "<<c<<endl;
		index_t mycount=0;
		long local_edge_error=0;
		double t2 = wtime();
#pragma omp parallel for num_threads(12) reduction(+:mycount) schedule(dynamic,512)
		for(vertex_t t =0; t<buffer->chunk_size; t++){
//cout<<buffer->edge[t].A<<" -- "<<buffer->edge[t].B<<endl;
			//tc
			vertex_t Ao = buffer->edge[t].A;
			vertex_t Bo = buffer->edge[t].B;
			
			if(Ao<csr1->start || Ao>csr1->end || Bo<csr2->start || Bo>csr2->end){
				local_edge_error++;
//				cout<<"edge not in partition"<<endl;
				continue;
			}

			//get the offset of CSR partitions
//			vertex_t offsetA = offset[i];
//			vertex_t offsetB = offset[j];
			
			vertex_t A = Ao-offsetA;
			vertex_t B = Bo-offsetB;	
//cout<<A<<" - "<<B<<endl;
			
			
			index_t m=csr1->begin[A+1] - csr1->begin[A];
			index_t n=csr2->begin[B+1] - csr2->begin[B];
//cout<<"the two degree of this edge = "<<m<<" "<<n<<endl;


			vertex_t *a = &csr1->adj[csr1->begin[A]];
			vertex_t *b = &csr2->adj[csr2->begin[B]];

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
//cout<<"find triangle "<<Ao<<" "<<Bo<<" "<<x<<endl;
				}
			}
		}
double t3 = wtime();	
time_compute  += t3-t2; 
		count += mycount;
		edge_error += local_edge_error;
		sub_edge_error += local_edge_error;
	}
	buffer-> Edge_free();
//	cout<<"subproblem "<<i<<"-"<<j<<"-"<<k<<"  ";
//	cout<<"subproblem edge error = "<<sub_edge_error<<endl;
//	cout<<"csr1 starts from "<<csr1->start<<" and end with "<<csr1->end<<endl;
//	cout<<"csr2 starts from "<<csr2->start<<" and end with "<<csr2->end<<endl;
		
	
}
