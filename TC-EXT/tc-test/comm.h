#ifndef	COMM_HEADER
#define	COMM_HEADER
//#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//--------------------------
//typedef 	unsigned long int 	index_t;
//typedef		unsigned int		vertex_t;

typedef 	long int 	index_t;
typedef		int		vertex_t;
//--------------------------------

//#define BufferSize	8388608//16777216 //16M edges, or 64MB*2 for edge buffer
#define BufferSize	67108864//16777216 //16M edges, or 64MB*2 for edge buffer



inline off_t fsize(const char *filename) {
	struct stat st;
	if (stat(filename, &st) == 0){
		return st.st_size;
	}
	return -1;
}

typedef struct EDGE{
	vertex_t A;
	vertex_t B;
} Edge;



#endif
