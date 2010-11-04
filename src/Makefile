include ../make.inc

SOURCE_FILES_CXX=Complex_2D.c++ FourierT.c++ Projection.c++ Config.c++
SOURCE_FILES_C=io.c
OBJECT_FILES=$(SOURCE_FILES_CXX:.c++=.o) $(SOURCE_FILES_C:.c=.o)
HEADER_FILES=$(SOURCE_FILES_CXX:.c++=.h) $(SOURCE_FILES_C:.c=.h)

SOURCE_TOOLS=hdf2ppm.c
EXEC_TOOLS=$(SOURCE_TOOLS:.c=.exe)

LIB_SO=../lib/libCOECXS.so
LIB_A=../lib/libCOECXS.a

#External libraries to include
LIBS := -lfftw3 -ltiff -L$(HDF_ROOT)/lib/ -ldf -lmfhdf -lstdc++
INC := -I/usr/include/hdf/ -I$(HDF_ROOT)/include/ 

all: $(LIB_A) $(LIB_SO) install_headers #$(EXEC_TOOLS)

$(LIB_A): $(OBJECT_FILES) *.h
	ar cr $(LIB_A)$  $(OBJECT_FILES)

$(LIB_SO): $(OBJECT_FILES) *.h
	$(LD) $(LDFLAGS) $(SOFLAGS) -o $(LIB_SO) $(OBJECT_FILES) $(LIBS)
#	gcc -shared -o ../lib/cygCOECXS.dll \
#	-Wl,--out-implib=../lib/libCOECXS.dll.a \
#	-Wl,--export-all-symbols \
#	-Wl,--enable-auto-import \
#	-Wl,--whole-archive $(OBJECT_FILES) \
#	-Wl,--no-whole-archive $(LIBS)



install_headers: $(HEADER_FILES)
	cp $(HEADER_FILES) ../include

%.exe: %.c
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@ $(LIBS) -L../lib/ -lCOECXS


%.o: %.c++
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@


%.o: %.c
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@


clean:                                                     
	rm -f *.o *~ 

clobber: clean
	rm -f *.exe 