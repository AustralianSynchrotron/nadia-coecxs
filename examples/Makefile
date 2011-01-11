export LD_RUN_PATH=:/usr/lib64/hdf:/home/ndavidson/cxs/cxs_software/lib

EXAMPLE_SRC=PlanarCDI_example.c PlanarCDI_simulation_example.c	\
FresnelCDI_example.c FresnelCDI_WF_example.c
EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

all: $(EXAMPLE_EXEC)

%.exe: %.c
	g++  -fPIC -I/home/ndavidson/cxs/cxs_software/include $< -o $@ \
	-L/home/ndavidson/cxs/cxs_software/lib -lCOECXS \
	-L/usr/lib64/hdf  -ltiff -lmfhdf -lz -ljpeg -ldf -lm -lfftw3 -lstdc++ 


clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe
