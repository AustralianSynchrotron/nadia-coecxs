cimport complex2d
cimport cio
cimport double2d
from complex2d cimport Complex_2D, MAG, IMAG
from double2d cimport Double_2D, PyDouble2D
from cython.operator cimport dereference as deref
from libcpp.string cimport string


cdef class PyComplex2D:
    def __cinit__(self, int nx, int ny):
        self.thisptr= new Complex_2D(nx,ny)
    def copyconstruct(self, PyComplex2D other):
        self.thisptr=new complex2d.Complex_2D(<Complex_2D>deref(other.thisptr))
    def __dealloc__(self):
        del self.thisptr
    def getSizeX(self):
        return self.thisptr.get_size_x()
    def getSizeY(self):
        return self.thisptr.get_size_y()
    def setValue(self,int x, int y ,int type, double value):
        self.thisptr.set_value(x,y,type,value)
    def setRead(self,int x,int y, double value):
        self.thisptr.set_real(x,y,value)
    def setImag(self,int x, int y, double value):
        self.thisptr.set_imag(x,y,value)
    def get2dMAG(self):
        result = PyDouble2D(self.getSizeX(),self.getSizeY()) 
        self.thisptr.get_2d(MAG,(deref(result.thisptr)))
        return result
    def get2dPHASE(self):
        result = PyDouble2D(self.getSizeX(),self.getSizeY()) 
        self.thisptr.get_2d(PHASE,(deref(result.thisptr)))
        return result
    def read_cplx(self, string filename):
        return cio.read_cplx(filename, deref(self.thisptr))   
    def write_cplx(self, string filename):
        return cio.write_cplx(filename, deref(self.thisptr))