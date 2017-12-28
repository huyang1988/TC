#include "worker.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include "comm.h"
#include <fstream>
#include <omp.h>
#include <math.h>


using namespace std;

worker::worker(string s){
	
//	gdata = new GPU_data[GPU_NUM];
//	for(int i=0; i<GPU_NUM; i++){
//		gdata[i].id = i;
//		gdata[i].EdgeBuffer = new Edge* [2];
//	}


	time_io = 0;
	time_cp = 0;
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

	subproblem = N1*(N1+1)*N2/2;

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
/*
	for(int k = 0; k<N2; k++){
		for(int i = 0; i<N1; i++){
			
			
			for(int j=0; j<N1; j++){
			
				proc(i,j,k);

			}
		}
	}
*/
	for(int i=0; i<subproblem; i++){
		proc(i);
	}
	cout<<"total number = "<<count<<endl;
	cout<<"edge error = "<<edge_error<<endl;
//	return;
}

void worker::proc(int t){
	if(t<0 || t>=subproblem){
		cout<<"invalid subproblem id "<<t<<endl;
		return;
	}

	int k = t/(N1*(N1+1)/2);
	int p = t%(N1*(N1+1)/2);
	int i = (int)(-0.5+sqrt(2*p+1));
	int j = p - i*(i+1)/2;

	
	proc(i,j,k);

}


void worker::proc(int i, int j, int k){
	count=0;
	double t0, t1;

	if(i!=j){	
	
		t0 = wtime();
		csr1->CSR_load(i,k);
		csr1->start = offset[i];
		csr1->end = csr1->start + csr1->beg_size -1;
					
		csr2->CSR_load(j,k);
		csr2->start = offset[j];
		csr2->end = csr2->start + csr2->beg_size -1;
					
		if(csr1->adj_size<=0 || csr2->adj_size<=0){
			return;
		}
		t1 = wtime();	
		time_io  += t1-t0; 
		
		t0 = wtime();		
		init_gpu();
		t1 = wtime();		
		time_cp  += t1-t0; 
			

		//proc for (i,j,k)
		buffer->Edge_open(i,j);
		if(buffer->edge_size>0){
			for(int c = 0; c<buffer->chunk_num; c++){
				t0 = wtime();			
				buffer->Edge_reload();
				t1 = wtime();	
				time_io  += t1-t0; 
				
				t0 = wtime();			
				call_gpu();
				t1 = wtime();			
				time_compute  += t1-t0; 
				count += gdata.host_count;

			}
			buffer-> Edge_free();
		}	
		//proc for (j,i,k)
		buffer->Edge_open(j,i);
		if(buffer->edge_size>0){
			for(int c = 0; c<buffer->chunk_num; c++){
				t0 = wtime();			
				buffer->Edge_reload();
				t1 = wtime();	
				time_io  += t1-t0; 
				
				t0 = wtime();			
				call_gpu_reverse();
				t1 = wtime();	
				
				time_compute  += t1-t0; 
				count += gdata.host_count;

			}
			buffer-> Edge_free();
		}	


		//free				
		csr1->CSR_free();
		csr2->CSR_free();

		free_gpu();
	}
	else{
		csr1->CSR_load(i,k);
		csr1->start = offset[i];
		csr1->end = csr1->start + csr1->beg_size -1;
	
		t0 = wtime();	
		init_gpu_internal();
		t1 = wtime();
		time_cp += t1-t0;

		buffer->Edge_open(i,j);
		if(buffer->edge_size>0){
			for(int c = 0; c<buffer->chunk_num; c++){
				t0 = wtime();			
				buffer->Edge_reload();
				t1 = wtime();	
				time_io  += t1-t0; 
				
				t0 = wtime();			
				call_gpu_internal();
				t1 = wtime();			
				
				time_compute  += t1-t0; 
				count += gdata.host_count;

			}
			buffer-> Edge_free();
		}	
		csr1->CSR_free();
		free_gpu_internal();
	}
		
	
}


