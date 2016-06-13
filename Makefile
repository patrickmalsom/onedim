# Makefile to compile the .so files used by ctypes in the HMC python code
# the potential definitions must be in the potential_defns directory
# make will compile all files in that directory

# compiler
CC=gcc

# find all .c files in the potential_defns directory
SRCS = $(wildcard potential_defns/*.c)
# list of the .c files with no extensions
PROGS = $(patsubst %.c,%.so,$(SRCS))

# test the OS type (mac or linux supported)
UNAME := $(shell uname)

# Linux build flags (with OpenMP enabled)
ifeq ($(UNAME), Linux) 
CFLAGS=-Wall -Warray-bounds -fPIC -fopenmp -march=native -O2 -g
#CFLAGS=-Wall -Warray-bounds -fPIC -march=native -O2 -g
LFLAGS=-fopenmp -lm
INCLUDEDIRS= -I/usr/include
LIBDIRS= -L/usr/lib64
endif

# Mac/Darwin build flags (no parallelization)
ifeq ($(UNAME), Darwin)
CFLAGS=-Wall -Warray-bounds -fPIC -g
LFLAGS=-lm 
INCLUDEDIRS= -I/opt/local/include
LIBDIRS= -L/opt/local/lib
endif

# make all .o files and then all .so files
# this will loop over all files in PROGS
all: $(PROGS)

%.o: %.c ; $(CC) $(INCLUDEDIRS) $(CFLAGS) -c $< -o $@
%.so: %.o ; $(CC) -shared $(LIBDIRS) $(LFLAGS) -o $@ $<

# clean all .o and .so files from potential_defns directory
clean: ; rm -f potential_defns/*.o potential_defns/*.so
