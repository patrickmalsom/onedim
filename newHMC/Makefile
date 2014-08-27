CC=gcc
CFLAGS=-Wall -Warray-bounds -fPIC -I/opt/local/include  -march=native -O2
LFLAGS=-shared -L/opt/local/lib  -lm -lgsl -lgslcblas -Wl,-install_name,NewMALA-func.so
#CFLAGS=-Wall -Warray-bounds -fPIC -fopenmp -march=native -O2
#LFLAGS=-shared -fopenmp -lm -lgsl -lgslcblas -Wl,-soname,NewMALA-func.so

all: NewMALA-func

NewMALA-func.o: NewMALA-func.c
	${CC} ${CFLAGS} ${EXTFLAGS} -c NewMALA-func.c

NewMALA-func: NewMALA-func.o 
	${CC} ${LFLAGS} -o NewMALA-func.so NewMALA-func.o
 
clean: 
	rm -f NewMALA-func.so NewMALA-func.o

