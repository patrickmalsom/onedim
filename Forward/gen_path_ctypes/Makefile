# Makefile to compile the .so files used by ctypes in the HMC python code

CC=gcc

CFLAGS=-Wall -Warray-bounds -fPIC -O2 -g
LFLAGS=-lgsl -lgslcblas -lm
INCLUDEDIRS= -I/usr/include -I./gsl-build/include 
LIBDIRS= -L/lib64 -L./gsl-build/lib 
#INCLUDEDIRS= -I/usr/include
#LIBDIRS= -L/lib64

all: forward_B

forward_B.o: forward_B.c ; $(CC) -static $(INCLUDEDIRS) $(CFLAGS) -c forward_B.c
forward_B : forward_B.o; $(CC) -static  forward_B.o -o forward_B.out $(LFLAGS) $(LIBDIRS)

clean: ; rm forward_B.o forward_B.out
