EXE=tc


CC= "$(shell which g++)"
FLAGS = -O3 -Wall -fopenmp -lpthread
#FLAGS = -O0 -g -fopenmp -lpthread
#CUFLAGS= -O3 --compiler-options -Wall -Xptxas -v -lpthread
#CUCC= "$(shell which nvcc)"

#CUFLAGS+=  -Xcompiler -fopenmp
#CUFLAGS+= -O0 -G -g -lpthread

ifeq ($(enable_monitor), 1)
#	CUFLAGS+= -DENABLE_MONITORING
endif

ifeq ($(enable_check), 1)
#	CUFLAGS+= -DENABLE_CHECKING
endif


OBJS=   main.o \
	graph.o\
	edge.o\
	wtime.o\

DEPS= 	Makefile \
		graph.h\
		edge.h\
		wtime.h\
		comm.h \


%.o:%.c $(DEPS)
	${CC} -c ${FLAGS} $< -o $@

${EXE}:${OBJS}
	${CC} ${OBJS} $(FLAGS) -o ${EXE}

clean:
	rm -rf *.o ${EXE}
