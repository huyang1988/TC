#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

using namespace std;

typedef struct packed_edge {
	  int v0;
	  int v1;
} packed_edge;

inline off_t fsize(const char *filename) {
    struct stat st; 
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1; 
}

		
main(int argc, char** argv){
	int fd,i;
	char* ss;
	char* ss_head;

	size_t file_size = fsize(argv[1]);
//	int N = atoi(argv[2]);
//	cout<<N<<endl;
			
	long int vert_count;
//	long int** bin = (long int**)malloc(N*sizeof(long int*));//ptr to the N new files


	fd=open( argv[1],O_CREAT|O_RDWR,00666 );



	ss_head = (char*)mmap(NULL,file_size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);

	size_t head_offset=0;
	while(ss_head[head_offset]=='#'){
		while(ss_head[head_offset]!='\n'){
			head_offset++;
		}
		head_offset++;
	}
	ss = &ss_head[head_offset];
	file_size -= head_offset;
//	char temp[50];
//	int a;

	size_t curr=0;
	size_t next=0;

	//step 1. vert_count,edge_count,
	size_t edge_count=0;
	while(next<file_size){
		if(ss[next]=='\n'){
			edge_count++;
		}
		next++;
	}
//cout<<"step1"<<endl;
//	edge_count++;

	cout<<"edge count: "<<edge_count<<endl;
	next=0;
	//step 2. each file size
	//step 3. write to bin file
	int fd1 = open( "output.bin",O_CREAT|O_RDWR,00666 );
	ftruncate(fd1, edge_count*sizeof(packed_edge));
	packed_edge* edge = (packed_edge*)mmap(NULL,edge_count*sizeof(packed_edge),PROT_READ|PROT_WRITE,MAP_SHARED,fd1,0);
	size_t offset=0;

//cout<<"step2"<<endl;
	while(offset<edge_count){
		char* sss=ss+curr;
		edge[offset].v0 = atoi(sss);
//		cout<<edge[offset].v0<<"  ";
		while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
			next++;
		}
		while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
			next++;
		}
		curr = next;

		char* sss1=ss+curr;
		edge[offset].v1 = atoi(sss1);
//		cout<<edge[offset].v1<<endl;
		while((ss[next]!=' ')&&(ss[next]!='\n')&&(ss[next]!='\t')){
			next++;
		}
		while((ss[next]==' ')||(ss[next]=='\n')||(ss[next]=='\t')){
			next++;
		}
		curr = next;

		offset++;
//cout<<"step2.5"<<endl;
	}


	munmap( ss,sizeof(char)*file_size );
	munmap( edge,sizeof(packed_edge)*edge_count );
}
