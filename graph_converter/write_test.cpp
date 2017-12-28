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
inline off_t fsize(const char *filename) {
    struct stat st; 
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1; 
}

main(int argc, char** argv){
	int fd,i;
	char* ss;

//	size_t file_size = fsize(argv[1]);
	fd=open( argv[1],O_CREAT|O_RDWR,00777 );
	if(ftruncate(fd, 100000*sizeof(char))){
		cout<<"ftruncate err\n";
	}
	ss = (char*)mmap(NULL,100000*sizeof(char),PROT_READ|PROT_WRITE,MAP_SHARED,fd,0);

	char* temp;
	int a;

	size_t curr=0;
	size_t next=0;
	for(int i=0;i<99;i++){
		ss[i]='4';
	}
	ss[100]='\0';

	munmap( ss,sizeof(char)*100 );
//	close(fd);
}
