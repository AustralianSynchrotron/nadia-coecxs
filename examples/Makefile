include ../make.inc

EXAMPLE_SRC=real_example.c simulation_example.c simple_example.c
EXAMPLE_EXEC=$(EXAMPLE_SRC:.c=.exe)

#LIBS := -lmfhdf -lfftw3 -ltiff -lprofiler ../lib/libCOECXS.so 
LIBS := -lmfhdf -lfftw3 -ltiff ../lib/libCOECXS.so 
#LIBS := ../lib/libCOECXS.so 
INC := -I../include 

all: $(EXAMPLE_EXEC)

%.exe: %.c
	$(CXX) $(CXXFLAGS) $(LIBS) $(INC) $< -o $@

clean:                                                     
	rm -f *~

clobber: clean
	rm -f *.exe