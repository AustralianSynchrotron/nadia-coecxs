export LD_RUN_PATH=:/usr/lib64/hdf:/data/nadia/cxs_software_dev/cxs_software/lib

IDL_SRC=IDL_interface.c
IDL_O=IDL_interface.o
LIB=libIDLCOECXS.so

all: $(LIB)

$(LIB): $(IDL_O)
	gcc -shared -Wl,-soname,libname.so.1 -o $(LIB) $(IDL_O) \
	-L/data/nadia/cxs_software_dev/cxs_software/lib -lCOECXS \
	-L/usr/lib64/hdf  -ltiff -lmfhdf -ldf -lz -ljpeg -lm -lfftw3 \
	-lstdc++  


%.o: %.c
	g++ -Wall -fPIC -c \
	-I/data/nadia/cxs_software_dev/cxs_software/include \
	-I/usr/local/itt/idl/idl80/external/include/ \
	$(IDL_SRC)

clean:                                                     
	rm -f *~ $(IDL_O)

clobber: clean
	rm -f $(LIB)
