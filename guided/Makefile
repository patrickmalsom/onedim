CC=gcc
CFLAGS=-Wall -Warray-bounds -fPIC -fopenmp -march=native -O2
LFLAGS=-shared -fopenmp -lm -lgsl -lgslcblas -Wl,-soname,GenPaths.so

all: GenPaths

GenPaths.o: GenPaths.c
	${CC} ${CFLAGS} ${EXTFLAGS} -c GenPaths.c

GenPaths: GenPaths.o 
	${CC} ${LFLAGS} -o GenPaths.so GenPaths.o
 
clean: 
	rm -f GenPaths.so GenPaths.o

