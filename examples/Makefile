
LD_RUN_PATH=$(PWD)/../lib

EXAMPLE_SRC=real_example.c simulation_example.c FCDI_example.c FCDI_WF_example.c
EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

LIB_A = ../lib/libCOECXS.a

all: $(EXAMPLE_EXEC)

%.exe: %.c
	g++ -g -O2 -I../include $< -o $@  -ltiff -lmfhdf -lm -lz -ljpeg -ldf -lfftw3 -lstdc++  -L../lib -lCOECXS

clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe
