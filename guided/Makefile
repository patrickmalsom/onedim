CC=gcc
CFLAGS=-Wall -Warray-bounds -fopenmp -march=native -O2 -g
LFLAGS=-fopenmp -lm -lgsl -lgslcblas

all: guidedGenPaths

guidedGenPaths.o: guidedGenPaths.c
	${CC} ${CFLAGS} ${EXTFLAGS} -c guidedGenPaths.c

guidedGenPaths: guidedGenPaths.o 
	${CC} ${LFLAGS} guidedGenPaths.o -o guidedGenPaths.out
 
clean: 
	rm -f guidedGenPaths.o guidedGenPaths.out


