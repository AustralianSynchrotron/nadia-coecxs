include ../make.inc

SOURCE_TOOLS=hdf2ppm.c tiff2ppm.c planar_CDI_reconstruction.c
EXEC_TOOLS=$(SOURCE_TOOLS:.c=.exe)

LIB_SO=../lib/libCOECXS.so
LIB_A=../lib/libCOECXS.a

#External libraries to include
#LIBS := -lfftw3 -ltiff -L$(HDF_ROOT)/lib/ -ldf -lmfhdf -lstdc++
#INC := -I/usr/include/hdf/ -I$(HDF_ROOT)/include/ -I../include/

all: $(LIB_A) $(LIB_SO) $(EXEC_TOOLS) 

$(LIB_A): ../src/*.c++ ../src/*.h ../src/*.c
	make -C ../src

$(LIB_SO): ../src/*.c++ ../src/*.h ../src/*.c
	make -C ../src

%.exe: %.c ../src/*.c++ ../src/*.h ../src/*.c
	g++ -g -O2 -I../include $< -o $@  -ltiff -lmfhdf -ldf -lfftw3 -lstdc++  $(LIB_A) 

#-L../lib/ -lCOECXS

clean:                                                     
	rm -f *.o *~ 

clobber: clean
	rm -f *.exe 
