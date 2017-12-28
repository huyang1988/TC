#include <iostream>
#include <stdio.h>
#include "wtime.h"
#include <errno.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>


typedef long index_t;
typedef long vertex_t;
typedef float path_t;
typedef long depth_t;
#define INFTY -1	

template<typename T, typename IT>
void print(T *data, IT count)
{
	IT i=0;
	for(i=0;i<count;++i)
		std::cout<<i<<" "<<data[i]<<"\n";
	std::cout<<"\n";
}

inline off_t fsize(const char *filename) {
    struct stat st; 
    if (stat(filename, &st) == 0)
        return st.st_size;
    return -1; 
}

depth_t bfs_mmap(
	index_t root,
	index_t *csr,
	index_t *beg_pos,
	depth_t *sa,
	path_t *sp_count,
	index_t vert_count
)
{
	sp_count[root] = 1;
	sa[root] = 0;
	depth_t level = 0;
	while(true)
	{
		double ltm= wtime();
		index_t front_count = 0;
		for(vertex_t vert_id = 0; vert_id<vert_count; vert_id++)
		{
			if(sa[vert_id] == level)
			{
				index_t my_beg = beg_pos[vert_id];
				index_t my_end = beg_pos[vert_id + 1];

				for(; my_beg < my_end; my_beg ++)
				{
					vertex_t nebr=csr[my_beg];
					if(sa[nebr] == INFTY)
					{
						sa[nebr] = level+1;
						front_count++;
					}
					
					if(sa[nebr] == level+1)
						sp_count[nebr]+=sp_count[vert_id];
				}
			}
		}

		if(front_count == 0) break;
	//	std::cout<<"Level "<<(int) level<<": "<<front_count<<" "
	//						<<wtime() - ltm<<"\n";
		level ++;
	}
//print<path_t, index_t>(sp_count, vert_count);
//exit(-1);
	return level+1;
}


void compute_bc(
	path_t *bc,
	path_t *tmp_bc,
	path_t *sp_count,
	index_t *csr,
	index_t *beg_pos,
	depth_t *sa,
	index_t root,
	index_t vert_count,
	depth_t depth_count
){
	for(depth_t i=depth_count;i>=0;--i)
	{
		for(vertex_t vert_id = 0; vert_id < vert_count; vert_id++)
		{
			if(sa[vert_id] == i)
			{
				depth_t my_depth=sa[vert_id];
				index_t my_beg = beg_pos[vert_id];
				index_t my_end = beg_pos[vert_id + 1];
				path_t tbc = 0;
				path_t tpath=sp_count[vert_id];
				for(; my_beg < my_end; my_beg ++)
				{
					vertex_t nebr=csr[my_beg];
					if(my_depth+1 == sa[nebr])
						tbc+=tpath*(1+tmp_bc[nebr])/sp_count[nebr];
				}
				tmp_bc[vert_id] = tbc;
			}
		}
	}
	tmp_bc[root]=0;
	for(vertex_t vert_id = 0; vert_id < vert_count; vert_id++)
		bc[vert_id] += tmp_bc[vert_id];
}

int main( int args, char **argv)
{
	std::cout<<"Input: ./exe /path/to/beg /path/to/csr \n";
	if(args != 3)
	{
		std::cout<<"Input format wrong\n";
		exit(-1);
	}
	const char *beg_filename = argv[1];
	const char *csr_filename = argv[2];
	size_t align=512;
	struct stat csr_st;
	index_t vert_count=fsize(beg_filename)/sizeof(vertex_t)-1;
		
	//first step: test async direct io
	//generate a frontier queue
	index_t* beg_pos;
	depth_t* sa;
	path_t* bc;
	path_t* sp_count;
	path_t* tmp_bc;

	if(	posix_memalign((void **)&beg_pos, align, sizeof(index_t)*vert_count+align) ||
			posix_memalign((void **)&sa, 			align, sizeof(depth_t)*vert_count) ||
			posix_memalign((void **)&bc, 			align, sizeof(path_t)*vert_count) ||
			posix_memalign((void **)&sp_count, align, sizeof(path_t)*vert_count)||
			posix_memalign((void **)&tmp_bc, align, sizeof(path_t)*vert_count))
	{
		perror("posix_memalign");
		return -1; 
	}

	int fd_beg = open(beg_filename, O_RDONLY | O_DIRECT);
	if(fd_beg==-1)
	{
		fprintf(stdout,"Wrong open %s\n",beg_filename);
		exit(-1);
	}
	
	//To investigate: alignment read
	index_t aligned_beg = sizeof(vertex_t)*(vert_count+1);
	if(aligned_beg & (align - 1)) 
		aligned_beg += align - (aligned_beg & (align - 1));
	
	if(read(fd_beg, beg_pos,aligned_beg) != (vert_count+1) * sizeof(index_t))
	{
		perror("read wrong\n");
		exit(-1);
	}
	close(fd_beg);
	int fd = open(csr_filename, O_RDONLY);
	if(fd==-1)
	{
		fprintf(stdout,"Wrong open %s\n",csr_filename);
		exit(-1);
	}
	
	fstat(fd, &csr_st);
	index_t *csr = (index_t*) mmap(0, csr_st.st_size, 
			PROT_READ, MAP_PRIVATE, fd, 0);
	perror("mmap");
	
	index_t report=1;
	memset(bc, 0, sizeof(path_t)*vert_count);
	double btm= wtime();
	for(index_t i=0;i<512;i++)
	{
		if(report <= i)
		{
			std::cout<<i<<" finished with "<<wtime()-btm<<" seconds\n";
			report <<=1;
		}

		memset(sa, INFTY, sizeof(depth_t)*vert_count);
		memset(sp_count, 0, sizeof(path_t)*vert_count);
		memset(tmp_bc, 0, sizeof(path_t)*vert_count);
		depth_t tdepth=bfs_mmap(i, csr, beg_pos, sa, sp_count, vert_count);
		compute_bc(bc, tmp_bc, sp_count, csr, beg_pos, sa, i, vert_count, tdepth);
	}
	std::cout<<"Total time: "<<wtime()-btm<<" seconds\n";
	std::cout<<"BC value: \n";
	print<path_t, index_t>(bc, vert_count);
	return 0;
}
