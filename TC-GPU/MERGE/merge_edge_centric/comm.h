#ifndef	COMM_HEADER
#define	COMM_HEADER
#include <iostream>
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

int GPU_NUM=1; 
int PART_NUM=2; 

inline off_t fsize(const char *filename) {
	struct stat st;
	if (stat(filename, &st) == 0){
		return st.st_size;
	}
	return -1;
}


static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", \
        cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define H_ERR( err ) \
  (HandleError( err, __FILE__, __LINE__ ))



//////////////////////////////////////////////////
//SCALE*THDS_NUMS*sizeof(int) should be 
//limited by the size of the shared memory
/////////////////////////////////////////////////
//////////////////////////////////////////////////
#define	THDS_NUM			256	
#define	BLKS_NUM			256	

#define	V_NON_INC			-1
#define V_INI_HUB			-2

#define VALIDATE_TIMES		1
#define NUM_SRC		 		1	
enum ex_q_t
{
  SML_Q,
  MID_Q,
  LRG_Q,
  NONE
};

typedef struct Gdata{
	cudaStream_t	*stream;
} *GPU_data;


#define	VIS				0x02
#define UNVIS			0x00
#define FRT				0x01
#define	SET_VIS(a)		((a)=0x02)

#define	SET_FRT(a)		((a)=0x01)

#define	IS_FRT(a)		((a)==0x01)
#define	IS_VIS(a)		((a)==0x02)
#define	IS_UNVIS(a)		((a)==0x00)

//----------------------------------
//GLOBAL VARIABLES
//---------------------------------
//--------------------------------
#define INFTY			255	
#endif

#ifndef EXTERN
#define EXTERN

#define HUB_SZ			1536	
#define HUB_BU_SZ		1920//should be 1.25 of HUB_SZ
								//since there is no status
								//array in __shared__ mem
#define HUB_CRITERIA	0	

#define	Q_CARD	3

#endif
