# Makefile to compile the .so files used by ctypes in the HMC python code

CC=gcc

CFLAGS=-shared -Wall -Warray-bounds -fPIC -march=native -O2 -g
LFLAGS=-lm -lgsl -lgslcblas
INCLUDEDIRS= -I/usr/include
LIBDIRS= -L/usr/lib64

all: forward_B

forward_B.o: forward_B.c ; $(CC) $(INCLUDEDIRS) $(CFLAGS) -c forward_B.c
forward_B : forward_B.o; $(CC) $(LIBDIRS) $(LFLAGS) forward_B.o -o forward_B.out

clean: ; rm forward_B.o forward_B.out
