#export LD_RUN_PATH=:/home/ndavidson/cxs/cxs_software/lib

IDL_SRC=IDL_interface.c
IDL_O=IDL_interface.o
LIB=libIDLCOECXS.so

all: $(LIB)

$(LIB): $(IDL_O)
	gcc -shared -g -O2 -o $(LIB) $(IDL_O) \
	-L/home/ndavidson/cxs/cxs_software/lib/ -lCOECXS \
	 -ltiff -lmfhdf -ldf -lz -ljpeg -lm -lfftw3 -lstdc++ 

$(IDL_O): $(IDL_SRC)
	g++ -g -O2 -c \
	-I/home/ndavidson/cxs/cxs_software/include \
	$(IDL_SRC)

clean:                                                     
	rm -f *~ $(IDL_O)

clobber: clean
	rm -f $(LIB)
