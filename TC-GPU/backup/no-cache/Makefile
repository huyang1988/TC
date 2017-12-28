EXE=tc

CUFLAGS= -O3 --compiler-options -Wall -Xptxas -v -lpthread
#CUFLAGS= -O3 -Wall -lpthread -openmp
#CUFLAGS+=-O3 --compiler-options -Wall -Xptxas -v -lpthread
CUCC= "$(shell which nvcc)"

#CUFLAGS= -arch=sm_35  ${COMMFLAGS}#-Xptxas -dlcm=cg#disable l1 cache
#CUFLAGS+= -ccbin=g++ -Xcompiler -fopenmp
CUFLAGS+=  -Xcompiler -fopenmp
#CUFLAGS+= -O0 -G -g -lpthread

ifeq ($(enable_monitor), 1)
#	CUFLAGS+= -DENABLE_MONITORING
endif

ifeq ($(enable_check), 1)
#	CUFLAGS+= -DENABLE_CHECKING
endif


OBJS=   main.o 
DEPS= 	Makefile \
		graph.c\
		sort.cu\
		scan.cu\
		comm.h \
#		expander.cuh \
		kernel.cu\
		inspector.cuh \
		graph.cuh \
		bfs_gpu_opt.cuh \
		wtime.h \
		validate.h \
		scan.cuh \
		allocator.cuh 

%.o:%.cu $(DEPS)
	${CUCC} -c  ${CUFLAGS} $< -o $@

${EXE}:${OBJS}
	${CUCC} ${OBJS} $(CUFLAGS) -o ${EXE}

clean:
	rm -rf *.o ${EXE}
