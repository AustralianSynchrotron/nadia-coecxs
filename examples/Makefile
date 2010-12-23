#include ../make.inc

LD_RUN_PATH=$(PWD)/../lib

EXAMPLE_SRC=real_example.c simulation_example.c FCDI_example.c FCDI_WF_example.c Angelas.c

EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

#LIBS := -L$(CXS_ROOT)/lib/ -lCOECXS -L$(HDF_ROOT)/lib/ -lmfhdf -lfftw3 -ltiff -lstdc++
#LIBS := -ldf -lmfhdf -lfftw3 -ltiff ../lib/libCOECXS.a # -lstdc++
#INC := -I../include 

LIB_A = ../lib/libCOECXS.a

all: $(EXAMPLE_EXEC)

%.exe: %.c
	g++ -g -O2 -I../include $< -o $@  -lprofiler -ltiff -lmfhdf -lm -lz -ljpeg -ldf -lfftw3 -lstdc++  -L../lib -lCOECXS

clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe
