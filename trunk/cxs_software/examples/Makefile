export LD_RUN_PATH=:/home/tdjempire/CXS/DIXE/lib

EXAMPLE_SRC=PlanarCDI_example.c PlanarCDI_simulation_example.c  \
	    FresnelCDI_example.c FresnelCDI_WF_example.c \
	    FresnelCDI_simulation_example.c ComplexConstraints_example.c \
	    PhaseDiverse_example.c PartialCDI_example.c PartialCDI_simulation_example.c\

#PlanarCDI_example.c PlanarCDI_simulation_example.c	\
  FresnelCDI_example.c FresnelCDI_WF_example.c \
  FresnelCDI_simulation_example.c ComplexConstraints_example.c \
  PhaseDiverse_example.c PartialCDI_example.c PartialFresnelCDI_WF_example.c\
  PartialFresnelCDI_W_example.c
EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

all: $(EXAMPLE_EXEC)

%.exe: %.c
	g++ -fPIC -O3 -ffast-math  -DHAVE_HDF_MFHDF_H -I/home/tdjempire/CXS/DIXE/include $< -o $@ \
	  -L/home/tdjempire/CXS/DIXE/lib -lCOECXS \
	   -ltiff -lmfhdf -ldf -lz -ljpeg -llapack -lm -lfftw3f -lstdc++ 


clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe
