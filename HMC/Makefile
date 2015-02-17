CC=gcc

# Name of the code to be compiled
CODENAME=onedimHMC

SOURCE=$(CODENAME).c
OBJECT=$(CODENAME).o
SHAREDLIB=$(CODENAME).so

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

all: $(CODENAME)

$(OBJECT): $(SOURCE) ; $(CC) $(INCLUDEDIRS) $(CFLAGS) -c $(SOURCE)

$(CODENAME): $(OBJECT) ; $(CC) -shared $(LIBDIRS) $(LFLAGS) -o $(SHAREDLIB) $(OBJECT) 

clean: ; rm -f $(OBJECT) $(SHAREDLIB)
