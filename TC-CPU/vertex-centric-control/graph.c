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
	heapAdj = new vertex_t[edge_count];
}

void graph::validation(){
	index_t mycount=0;
//cout<<"validation "<<upperEdgeCount<<endl;
#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(dynamic,1024)
	for(index_t i=0; i<upperEdgeCount; i++){
		vertex_t A=upperHead[i];
		vertex_t B=upperAdj[i];
		index_t m=upperBegin[A+1]-upperBegin[A];
		index_t n=upperBegin[B+1]-upperBegin[B];
cout<<"launch intersection for edge "<<A<<" "<<B<<endl;

//cout<<"edge: "<<i<<" "<<U<<"-"<<V<<" ";
cout<<"degree: "<<m<<" "<<n<<endl;

		vertex_t *a = &upperAdj[upperBegin[A]];
		vertex_t *b = &upperAdj[upperBegin[B]];

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
cout<<"A-B-C : "<<A<<" "<<B<<" "<<x<<endl;
			}
		}
	}
	cout<<"merge version tc = "<<mycount<<endl;
}

void graph::bsvalidation(){
	index_t mycount=0;
#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(dynamic,512)
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

		for(index_t j=0; j<n; j++){
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
	cout<<"bs validation version tc = "<<mycount<<endl;
}

	

index_t bsp(vertex_t *a, index_t m, vertex_t x){
	
	vertex_t *base =a;
	index_t n=m;
	while (n > 1) {
		const index_t half = n / 2;
		base = (base[half] < x) ? &base[half] : base;
		n -= half;
	}
	//return value is the first position larger than or equal to x.
	return (*base < x) + base - a;
}

// intersection function for m and n larger than 256?
index_t intersec256(vertex_t * a, vertex_t * b, index_t m, index_t n){
	if(m<n){
		vertex_t * tempv;
		index_t    temp;
		tempv=a;
		a=b;
		b=tempv;
		temp=m;
		m=n;
		n=temp;
	}
	index_t mycount = 0;
	//cache 31 vertices 
	vertex_t cache[P+1];	//cache content in array A
	index_t  index[P+1];	//cache index in array A
	index_t  part[P+1];	//partition index in array B
	for(int i=0;i<P+1;i++){
		index[i] = m * i/(P);
		cache[i] = a[index[i]];
		part[i]=0;
	}
	part[P]=n;

	
	//phase 1: search-partitioning
	index_t step =P/2;
	while(step>0){
		index_t j=0;
		index_t k=step;
		while(k<P){
			part[k]=bsp(&b[part[k-step]],part[k+step]-part[k-step],cache[k])+part[k-step];
			j++;
			k=(2*j+1)*step;
		}
		step /=2;
	} 
		
	//phase 2: search in each partition
	for(int i=0;i<P;i++){
		for(int j=part[i]; j<part[i+1]; j++){//last one is not include
			vertex_t X = b[j];
		
			index_t max=index[i+1]-index[i];
			vertex_t * base = &a[index[i]];
			vertex_t * temp = base;
			while (max > 1) {
				const index_t half = max/2;
				base = (base[half] < X) ? &base[half] : base;
				max -= half;
			}
			if(temp[(*base < X) + base - temp]==X && base<(temp+m-1)){
				mycount++;
			}
			
		}
	}

	return mycount;
}

index_t intersec1(vertex_t *a, vertex_t *b, index_t m, index_t n){
	if(m<n){
		vertex_t * tempv;
		index_t    temp;
		tempv=a;
		a=b;
		b=tempv;
		temp=m;
		m=n;
		n=temp;
	}
	index_t mycount=0;
	index_t bot;
	index_t top;
	index_t r;

	for(index_t j=0; j<n; j++){
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
	return mycount;
}

index_t intersecmerge(vertex_t *a, vertex_t *b, index_t m, index_t n){
	index_t mycount=0;
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
	return mycount;
}


void graph::part_validation(){
//cout<<"validation "<<upperEdgeCount<<endl;
//
	long int global_count=0;
	for(int p=0; p<PART_NUM; p++){
	
		index_t mycount=0;
//#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(dynamic,1024)
#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(static,1024)
		for(index_t i=0; i<upperEdgeCount; i++){
			//get the edge A-B from global graph
			vertex_t A=upperHead[i];
			vertex_t B=upperAdj[i];

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
cout<<"A-B-C : "<<A<<" "<<B<<" "<<x<<endl;
				}
			}
		}
	
		global_count += mycount;
	}
	cout<<"merge version tc = "<<global_count<<endl;
}


void graph::vertex_centric_cpu(){
	index_t mycount=0;
cout<<"vertex number: "<<vert_count<<endl;
//cout<<"validation "<<upperEdgeCount<<endl;
#pragma omp parallel for num_threads(56) reduction(+:mycount) schedule(dynamic,1024)
	for(index_t i=0; i<vert_count; i++){
		vertex_t A=i;
		vertex_t *a = &upperAdj[upperBegin[A]];
		index_t m=upperBegin[A+1]-upperBegin[A];

		for(index_t j=upperBegin[A]; j<upperBegin[A+1]; j++){
			
		
			vertex_t B=upperAdj[j];
			index_t n=upperBegin[B+1]-upperBegin[B];
			vertex_t *b = &upperAdj[upperBegin[B]];
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
	}
	cout<<"node-iterator version tc = "<<mycount<<endl;
}
