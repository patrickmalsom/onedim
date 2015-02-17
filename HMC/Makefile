
# define compiler
CC=gcc

# get the source files from the potentials directory
# This will compile all files in that directory
VPATH = potential_defns
SRCS = $(wildcard potential_defns/*.c)
PROGS = $(patsubst %.c,%.so,$(SRCS))

UNAME := $(shell uname)

ifeq ($(UNAME), Linux) # linux build flags
CFLAGS=-Wall -Warray-bounds -fPIC -fopenmp -march=native -O2 -g
LFLAGS=-fopenmp -lm
INCLUDEDIRS= -I/usr/include
LIBDIRS= -L/usr/lib64
endif

ifeq ($(UNAME), Darwin) # mac build flags
CFLAGS=-Wall -Warray-bounds -fPIC -g
LFLAGS=-lm -Wl,-install_name,$(CODENAME).so
INCLUDEDIRS= -I/opt/local/include
LIBDIRS= -L/opt/local/lib
endif

all: $(PROGS)

%.o: %.c ; $(CC) $(INCLUDEDIRS) $(CFLAGS) -c $< -o $@
%.so: %.o ; $(CC) -shared $(LIBDIRS) $(LFLAGS) -o $@ $<

clean: ; rm -f potential_defns/*.o potential_defns/*.so
