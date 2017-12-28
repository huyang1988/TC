EXE=tc


COMMFLAGS= -O3 --compiler-options -Wall -Xptxas -v -lpthread
CUCC= "$(shell which nvcc)"
CC= "$(shell which mpicxx)"

CUFLAGS= -arch=sm_35  ${COMMFLAGS}#-Xptxas -dlcm=cg#disable l1 cache
CUFLAGS+= -ccbin=g++ -Xcompiler -fopenmp
#CUFLAGS+= -O0 -G -g -lpthread
MPCFLAGS  = -Wall -I"$(shell dirname $(CUCC))/../include" -L"$(shell dirname $(CUCC))/../lib64" -lcudart -fopenmp



OBJS=  	main.o \
		scan.o\
		wtime.o\
		graph.o\
		worker.o\

DEPS= 	Makefile \
		graph.h\
		worker.h\
		comm.h\
		wtime.h\
#		cuUtil.cuh\
		graph.cuh \

%.o:%.c $(DEPS)
	${CC} -c  ${MPCFLAGS} $< -o $@

%.o:%.cu $(DEPS)
	${CUCC} -c  ${CUFLAGS} $< -o $@

${EXE}:${OBJS}
	${CC} ${OBJS} $(MPCFLAGS) -o ${EXE}

clean:
	rm -rf *.o ${EXE}



test1:
	mpirun -N 1 -n 1 ./tc ~/data/orkut/
