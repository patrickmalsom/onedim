# Makefile to compile the .so files used by ctypes in the HMC python code

CC=gcc

CFLAGS=-shared -Wall -Warray-bounds -fPIC -march=native -O2 -g
LFLAGS=-lm
INCLUDEDIRS= -I/usr/include
LIBDIRS= -L/usr/lib64

all: calc_move.so
%.o: %.c ; $(CC) $(INCLUDEDIRS) $(CFLAGS) -c calc_move.c -o calc_move.o
%.so: %.o ; $(CC) -shared $(LIBDIRS) $(LFLAGS) -o calc_move.so calc_move.o

clean: ; rm calc_move.so
