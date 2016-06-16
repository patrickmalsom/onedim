# Makefile to compile the .so files used by ctypes in the HMC python code

CC=gcc

CFLAGS=-shared -Wall -Warray-bounds -fPIC -march=native -O2 -g
LFLAGS=-lm -lgsl -lgslcblas
INCLUDEDIRS= -I/usr/include
LIBDIRS= -L/usr/lib64

all: calc_move

calc_move.o: calc_move.c ; $(CC) $(INCLUDEDIRS) $(CFLAGS) -c calc_move.c
calc_move : calc_move.o; $(CC) $(LIBDIRS) $(LFLAGS) calc_move.o -o calc_move.out

clean: ; rm calc_move.o calc_move.out
