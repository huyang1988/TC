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
	
//	gdata = new GPU_data[GPU_NUM];
//	for(int i=0; i<GPU_NUM; i++){
//		gdata[i].id = i;
//		gdata[i].EdgeBuffer = new Edge* [2];
//	}


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

	subproblem = N1*N1*N2;

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
	int i = t%N1;
	t = t/N1;
	int j = t%N1;
	t = t/N1;
	proc(i,j,t);

	cout<<"task "<<t<<" with count "<<gdata.host_count<<endl;
	cout<<"accumulative count "<<count<<endl;

}


void worker::proc(int i, int j, int k){

	double t0 = wtime();			
	
	csr1->CSR_load(i,k);
	csr1->start = offset[i];
	csr1->end = csr1->start + csr1->beg_size -1;
				
	csr2->CSR_load(j,k);
	csr2->start = offset[j];
	csr2->end = csr2->start + csr2->beg_size -1;
				
	double t1 = wtime();	
	time_io  += t1-t0; 
	
	buffer->Edge_open(i,j);
				
	if(csr1->adj_size<=0 || csr2->adj_size<=0 || buffer->edge_size<=0){
		return;
	}	
	
	
	init_gpu();

//	cout<<"subproblem "<<i<<"-"<<j<<"-"<<k<<endl;
//	vertex_t offsetA = offset[i];
//	vertex_t offsetB = offset[j];



	for(int c = 0; c<buffer->chunk_num; c++){
		double t0 = wtime();			
		buffer->Edge_reload();
		double t1 = wtime();	
	time_io  += t1-t0; 
//cout<<"chunk "<<c<<endl;
//		index_t mycount=0;
		double t2 = wtime();
		
		//input data: csr1, csr2, buffer
		//input parameters: buffer->chunksize, offset(A and B for csr1 and csr2)
		//output: return triangle count of this chunk to mycount
			
		/*fault tolerance
			if(Ao<csr1->start || Ao>csr1->end || Bo<csr2->start || Bo>csr2->end){
				local_edge_error++;
				cout<<"edge not in partition"<<endl;
				continue;
			}
		*/

		call_gpu();
		
		double t3 = wtime();	
		time_compute  += t3-t2; 
		count += gdata.host_count;
//		count = gdata.host_count;//when using MPI, it does not keep a local count anymore, instead, it send count of current task back to the master.

	}


	buffer-> Edge_free();
			
	csr1->CSR_free();
	csr2->CSR_free();

	free_gpu();
		
	
}


