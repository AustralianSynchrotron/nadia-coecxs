

CXX      = g++
CXXFLAGS = -O3 -g -Wall -fno-exceptions -fPIC $(ROOTCFLAGS)
#CXXFLAGS = -O4 -g -Wall -fno-exceptions -funroll-loops -fPIC $(ROOTCFLAGS)

LD       = gcc
LDFLAGS  = #-fno-stack-protector -fno-stack-protector-all
SOFLAGS  = -shared -Wl


#External libraries to include
#LIBS := -lmfhdf -lfftw3 -ltiff -lprofiler 
#INC := -I/usr/include/hdf/ 

