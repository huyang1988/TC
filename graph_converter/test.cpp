#include <iostream>
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
	  long int v0;
	  long int v1;
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

	size_t file_size = fsize(argv[1]);
//	int N = atoi(argv[2]);
//	cout<<N<<endl;
	
	fd=open( argv[1],O_CREAT|O_RDWR,00777 );
cout<<"file_size  "<<file_size<<endl;


	ss = (char*)mmap(NULL,file_size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);

	char* temp;
long	int a;

	size_t curr=0;
	size_t next=0;

	while(next<file_size){
		if(ss[next]==' '){
			//space
			while(ss[next]==' '||ss[next]=='\n'){
				next++;
			}
			char* sss=ss+curr;
	//		strncpy(temp,sss,next-curr);
			a=atol(sss);
			cout<<"atoi"<<a<<endl;
			curr = next;
      		

		}
		else if(ss[next]=='\n'){
			//change line
			
			while(ss[next]==' '||ss[next]=='\n'){
				next++;
			}
			next++;
			char* sss=ss+curr;
		//	strncpy(temp,sss,next-curr);
			a=atol(sss);
			cout<<"atoi"<<a<<endl;
			curr = next;
		}
		else{
			next++;
		}
//cout<<"step2.5 "<<next<<endl;
	}

//	for(i = 0; i<file_size; i++){
//		cout<<*(p_map+i);
//	}
//	cout<<endl;
	
//	cout<<sline<<endl;
//	std::string ss;
//	getline(sline,ss);
//	int a;
//	sline>>a;
//	cout<<a<<endl;

	munmap( ss,sizeof(char)*file_size );
}
