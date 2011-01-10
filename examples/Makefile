export LD_RUN_PATH=;$(PWD)/../lib

EXAMPLE_SRC=PlanarCDI_example.c PlanarCDI_simulation_example.c	\
FresnelCDI_example.c FresnelCDI_WF_example.c
EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

LIB_A = ../lib/libCOECXS.a

all: $(EXAMPLE_EXEC)

%.exe: %.c
	g++ -g -O2 -I../include $< -o $@ \
	-L../lib/ -lCOECXS \
	 -ltiff -lmfhdf -lz -ljpeg -ldf -lm -lfftw3 -lstdc++ 


clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe
