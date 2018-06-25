EXE=rank_by_degree.bin


COMMFLAGS= -O3 --compiler-options -Wall -Xptxas -v -lpthread
CC= "$(shell which g++)"

MPCFLAGS  = -Wall -fopenmp -O3



OBJS=  	main.o \
		graph.o\
#		sort.o\

DEPS= 	Makefile \
		graph.h\
		comm.h\
#		graph.cuh \

%.o:%.c $(DEPS)
	${CC} -c  ${MPCFLAGS} $< -o $@

#%.o:%.cu $(DEPS)
#	${CUCC} -c  ${CUFLAGS} $< -o $@

${EXE}:${OBJS}
	${CC} ${OBJS} $(MPCFLAGS) -o ${EXE}

clean:
	rm -rf *.o ${EXE}



test1:
	mpirun -N 1 -n 1 ./tc ~/data/orkut/
