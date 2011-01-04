#include ../make.inc

SOURCE_FILES_CXX=Complex_2D.c++ FourierT.c++ \
	PlanarCDI.c++ Config.c++ FresnelCDI_WF.c++ FresnelCDI.c++
SOURCE_FILES_C=io_hdf.c io_ppm.c io_tiff.c io_dbin.c \
	io_cplx.c
OBJECT_FILES=$(SOURCE_FILES_CXX:.c++=.o) $(SOURCE_FILES_C:.c=.o)
HEADER_FILES=$(SOURCE_FILES_CXX:.c++=.h) io.h Double_2D.h #$(SOURCE_FILES_C:.c=.h)

LIB_SO=../lib/libCOECXS.so
LIB_A=../lib/libCOECXS.a

#External libraries to include
#LIBS := -lfftw3 -ltiff -L$(HDF_ROOT)/lib/ -ldf -lmfhdf -lstdc++
#INC := -I/usr/include/hdf/ -I$(HDF_ROOT)/include/ 

all: $(LIB_A) $(LIB_SO) install_headers #$(EXEC_TOOLS)

$(LIB_A): $(OBJECT_FILES) *.h
	ar cr $(LIB_A)$  $(OBJECT_FILES)

$(LIB_SO): $(OBJECT_FILES) *.h
	g++ -O3 -shared -o $(LIB_SO) $(OBJECT_FILES) # -ltiff -lmfhdf -lm -lz -ljpeg -ldf -lfftw3 -lstdc++ 

#	gcc -shared -o ../lib/cygCOECXS.dll \
#	-Wl,--out-implib=../lib/libCOECXS.dll.a \
#	-Wl,--export-all-symbols \
#	-Wl,--enable-auto-import \
#	-Wl,--whole-archive $(OBJECT_FILES) \
#	-Wl,--no-whole-archive $(LIBS)



install_headers: $(HEADER_FILES)
	cp $(HEADER_FILES) ../include

%.o: %.c++
	g++ -g -O3 -c $< -o $@


%.o: %.c
	g++ -g -O3 -c $< -o $@


clean:                                                     
	rm -f *.o *~ 

clobber: clean
	rm -f *.exe 
