CXS_ROOT = $(HOME)/cxs/cxs_software
HDF_ROOT = $(HOME)/hdf4.2.5-cygwin-szip

CXX      = g++
CXXFLAGS = -O3 -g -Wall -fno-exceptions # -mno-cygwin


LD       = gcc
LDFLAGS  = #-fno-stack-protector -fno-stack-protector-all
SOFLAGS  = -shared -Wl


#External libraries to include
#LIBS := -lmfhdf -lfftw3 -ltiff -lprofiler 
#INC := -I/usr/include/hdf/ 

