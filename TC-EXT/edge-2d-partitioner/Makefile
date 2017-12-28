EXE=tc

COMMFLAGS= -O3 -compiler-options -Wall -Xptxas -v -lpthread
#CUFLAGS+=-O3 --compiler-options -Wall -Xptxas -v -lpthread
CUCC= "$(shell which nvcc)"

#CUFLAGS= -arch=sm_35  ${COMMFLAGS}#-Xptxas -dlcm=cg#disable l1 cache
#CUFLAGS+= -ccbin=g++ -Xcompiler -fopenmp
#CUFLAGS+= -O0 -G -g -lpthread

ifeq ($(enable_monitor), 1)
#	CUFLAGS+= -DENABLE_MONITORING
endif

ifeq ($(enable_check), 1)
#	CUFLAGS+= -DENABLE_CHECKING
endif


OBJS=  	main.o 
DEPS= 	Makefile \
		graph.c\
		sort.c\
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
test:
	./tc ~/rmat_21_128.dat.json.sorted &>>log/log
	nvprof --metrics warp_execution_efficiency ./tc ~/rmat_21_128.dat.json.sorted &>>log/log1
	nvprof --metrics gld_transactions_per_request ./tc ~/rmat_21_128.dat.json.sorted &>>log/log2
	nvprof --metrics global_cache_replay_overhead ./tc ~/rmat_21_128.dat.json.sorted &>>log/log3
	nvprof --metrics shared_replay_overhead ./tc ~/rmat_21_128.dat.json.sorted &>>log/log4
	nvprof --metrics shared_load_transactions ./tc ~/rmat_21_128.dat.json.sorted &>>log/log5
test1:
	nvprof --metrics l1_cache_global_hit_rate ./tc ~/rmat_21_128.dat.json.sorted &>>log/log6
test2:
	nvprof --metrics l1_cache_local_hit_rate ./tc ~/rmat_21_128.dat.json.sorted &>>log/log6
