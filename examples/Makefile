include ../make.inc

EXAMPLE_SRC=real_example.c simulation_example.c simple_example.c
EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

#LIBS := -L$(CXS_ROOT)/lib/ -lCOECXS -L$(HDF_ROOT)/lib/ -lmfhdf -lfftw3 -ltiff -lstdc++
LIBS := -ldf -lmfhdf -L../lib/ -lCOECXS -lfftw3 -ltiff # -lstdc++
INC := -I../include 

all: $(EXAMPLE_EXEC)

%.exe: %.c
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@ $(LIBS)

clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe