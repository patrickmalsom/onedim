CC=gcc
#CFLAGS=-Wall -Warray-bounds -Wl,-soname,guided -fPIC -fopenmp -march=native -O2 -g 
CFLAGS=-Wall -Warray-bounds -fPIC -fopenmp -march=native -O2
#LFLAGS=-fopenmp -lm -lgsl -lgslcblas -shared 
LFLAGS=-shared -fopenmp -lm -lgsl -lgslcblas -Wl,-soname,guidedGenPaths.so

# this might have worked...
#gcc -Wall -fPIC -c guidedGenPaths.c 
# gcc -shared -lm -lgsl -Wl,-soname,guidedGenPaths.so -glscblas -o guidedGenPaths.so guidedGenPaths.o

all: guidedGenPaths

guidedGenPaths.o: guidedGenPaths.c
	${CC} ${CFLAGS} ${EXTFLAGS} -c guidedGenPaths.c

guidedGenPaths: guidedGenPaths.o 
	${CC} ${LFLAGS} -o guidedGenPaths.so guidedGenPaths.o
 
clean: 
	rm -f guidedGenPaths.so guidedGenPaths.o

