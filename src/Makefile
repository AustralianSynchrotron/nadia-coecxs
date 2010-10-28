include ../make.inc

SOURCE_FILES_CXX=Complex_2D.c++ FourierT.c++ Projection.c++
SOURCE_FILES_C=io.c
OBJECT_FILES=$(SOURCE_FILES_CXX:.c++=.o) $(SOURCE_FILES_C:.c=.o)
HEADER_FILES=$(SOURCE_FILES_CXX:.c++=.h) $(SOURCE_FILES_C:.c=.h)

SOURCE_TOOLS=hdf2ppm.c
EXEC_TOOLS=$(SOURCE_TOOLS:.c=.exe)


LIB_SO=../lib/libCOECXS.so
LIB_A=../lib/libCOECXS.a

#External libraries to include
LIBS := -lmfhdf -lfftw3 -ltiff -lprofiler 
INC := -I/usr/include/hdf/ 


all: $(LIB_SO) $(LIB_A) $(EXEC_TOOLS) install_headers

$(LIB_A): $(OBJECT_FILES) *.h
	ar cr $(LIB_A)$  $(OBJECT_FILES)

$(LIB_SO): $(OBJECT_FILES) *.h
	$(LD) $(LDFLAGS) $(SOFLAGS) -o $(LIB_SO) $(OBJECT_FILES)

install_headers: $(HEADER_FILES)
	cp $(HEADER_FILES) ../include

%.exe: %.c
	$(CXX) $(CXXFLAGS) $(LIBS) $(INC) $< -o $@


%.o: %.c++
	$(CXX) $(CXXFLAGS) $(LIBS) $(INC) -c $< -o $@


%.o: %.c
	$(CXX) $(CXXFLAGS) $(LIBS) $(INC) -c $< -o $@


clean:                                                     
	rm -f *.o *~ 

clobber: clean
	rm -f *.exe 