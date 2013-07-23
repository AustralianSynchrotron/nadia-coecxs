cimport numpy as np
cimport cio as cio
cimport double2d
#from double2d cimport Double_2D,PyDouble2D
from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool

cdef class PyDouble2D:
    #cdef Double_2D *thisptr
    #def __cinit__(self):
    #    self.thisptr= new Double_2D()
    #    if self.thisptr == NULL:
    #        raise MemoryError('Not enough memory.')
    def __cinit__(self, int x_size=0, int y_size=0):
        if (x_size==0 and y_size==0):
           self.thisptr=new Double_2D()
        else:
           self.thisptr = new Double_2D(x_size,y_size)
        if self.thisptr == NULL:
            raise MemoryError('Not enough memory.')
    def __dealloc__(self):
        del self.thisptr
    def getMax(self):
        return self.thisptr.get_max()
    def getShape(self):
        return (self.thisptr.get_size_x(),self.thisptr.get_size_y())
    def get(self,int i, int j):
        return self.thisptr.get(i,j)
    def getArray(self):
        shape = self.getShape()
        arr=np.empty(shape)
        for i in range(0, shape[0]):
            for j in range(0,shape[1]):
                arr[i][j]=self.thisptr.get(i,j)
        return arr
    def read_tiff(self,string filename):
        return cio.read_tiff(filename,deref(self.thisptr) ) 
    def read_dbin(self,string filename,int nx, int ny):
        return cio.read_dbin(filename,nx,ny,deref(self.thisptr))
    def read_spec(self,string filename):
        return cio.read_spec(filename, deref(self.thisptr))
    def read_ppm(self,string filename):
        return cio.read_ppm(filename,deref(self.thisptr))
    def read_hdf4(self, string filename):
        return cio.read_hdf4(filename,deref(self.thisptr))
    def write_ppm(self,string filename,bool log_scale=False,double min=0,double max=0):
        if min==0 and max==0:
           if log_scale:
              return cio.write_ppm(filename, deref(self.thisptr),log_scale)
           else:
              return cio.write_ppm(filename,deref(self.thisptr))
        else:
           return cio.write_ppm(filename, deref(self.thisptr),log_scale,min,max)
    #def write_spec(self, string filename,int nx=0, int ny=0):
    #    if nx==0 or ny==0:
    #       return cio.write_spec(filename, deref(self.thisptr))
   #     else:
   #        return cio.write_spec(filename, deref(self.thisptr), nx, ny)
    def write_tiff(self,string filename,bool log_scale=False, double min=0,double max=0):
        if min==0 and max==0:
           if log_scale:
              return cio.write_tiff(filename, deref(self.thisptr),log_scale)
           else:
              return cio.write_tiff(filename,deref(self.thisptr))
        else:
           return cio.write_tiff(filename, deref(self.thisptr),log_scale,min,max)
    def write_dbin(self, string filename):
        return cio.write_dbin(filename, deref(self.thisptr))
    def getSizeX(self):
        return self.thisptr.get_size_x()
    def getSizeY(self):
        return self.thisptr.get_size_y()    