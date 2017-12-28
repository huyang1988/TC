EXE=partition-x


CC= "$(shell which g++)"
FLAGS = -O3 -Wall -fopenmp -lpthread
#FLAGS = -O0 -g



OBJS=   main.o \
	graph.o\
	wtime.o\

DEPS= 	Makefile \
		graph.h\
		wtime.h\
		comm.h \


%.o:%.c $(DEPS)
	${CC} -c ${FLAGS} $< -o $@

${EXE}:${OBJS}
	${CC} ${OBJS} $(FLAGS) -o ${EXE}

clean:
	rm -rf *.o ${EXE}
