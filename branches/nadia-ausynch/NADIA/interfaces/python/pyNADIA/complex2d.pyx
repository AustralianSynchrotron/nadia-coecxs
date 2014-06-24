"""!
This package provides python wrappers to the Complex_2D class. 
"""
cimport cio
from complex2d cimport Complex_2D, MAG, IMAG
from double2d cimport Double_2D, PyDouble2D
from cython.operator cimport dereference as deref
from libcpp.string cimport string
import numpy as np
cimport numpy as np

cdef class PyComplex2D:
    """!A class wrapping the Complex_2D C class. 
    
    Contains a pointer to a Complex_2D class and provides python only versions
    of the public methods of the C class. For further details on each method, see the corresponding method in Complex_2D.
    """
    
    def __cinit__(self, int nx=0, int ny=0):
        if nx==0 and ny ==0:
            self.thisptr = new Complex_2D()
        else:
            self.thisptr = new Complex_2D(nx, ny)
        self.numpy_flag=0

#    def set_array_ptr(self,np.ndarray[np.complex128_t, ndim=2, mode="c"] input):
    def set_array_ptr(self,double complex[:,:] input):
       # cdef np.ndarray[np.complex128_t, ndim=2, mode="c"] input_c = input
        self.numpy_flag=1
        self.thisptr.set_array_ptr( <void *> &input[0,0],input.shape[0],input.shape[1])

    def copyconstruct(self, PyComplex2D other):
        self.thisptr = new Complex_2D(< Complex_2D > deref(other.thisptr))

    def __dealloc__(self):
        if self.numpy_flag==1:
            self.thisptr.unset_sizes()
        if self.thisptr:
            del self.thisptr

    def getSizeX(self):
        """!Get the number of pixels in the x co ordinate.
        Corresponds to Complex_2D.get_size_x()
        @return An int, the number of pixels in the x coordinate.
        """
        return self.thisptr.get_size_x()

    def getSizeY(self):
        """!Get the number of pixels in the y co-ordinate.
        Corresponds to Complex_2D.get_size_y()
        @return An int, the number of pixels in the y coordinate.
        """
        return self.thisptr.get_size_y()
    
    def getReal(self, x, y):
        """!Get the real part of the value at location x,y
        @param x
        @param y
        """
        return self.thisptr.get_real(x,y)
    def getImag(self,x,y):
        """!Get the imaginary part of the value at location x,y
        @param x
        @param y
        """
        return self.thisptr.get_imag(x,y)
    def getMag(self,x,y):
        """!Get the magnitude of the value at location x,y
        @param x
        @param y
        """
        return self.thisptr.get_mag(x,y)
    def getPhase(self, x, y):
        """!Get the phase of the value at location x,y
        @param x
        @param y
        """
        return self.thisptr.get_phase(x,y)
    def getValue(self,x,y,type):
        """!Get the value at location x,y where the type is one of real, imaginary, magnitude or phase 
        @param x
        @param y
        @param type REAL,IMAG,MAG,PHASE
        """
        return self.thisptr.get_value(x,y,type)

    def setValue(self, x, y , type, value):
        """!Set a value at the x,y coordinate, type can be REAL,IMAG,MAG,PHASE. See the examples for how to use. Setting 
        magnitude will preserve the phase.
        Coresponds to Complex_2D.set_value()
        
        @param x An int, the x coordinate
        @param y An int, the y coordinate
        @param type REAL,IMAG,MAG or PHASE depending on which part of the complex value to set
        @param value A double, the value to set
        """
        self.thisptr.set_value(x, y, type, value)

    def setReal(self, x, y, value):
        """!Sets the real part of the complex at coordinate x,y.
        
        @param x An int, the x-coordinate to set.
        @param y An int, the y coordinate.
        @param value A double, the value to set.
        """
        self.thisptr.set_real(x, y, value)

    def setImag(self, x, y, value):
        """!Sets the imaginary part of the complex at coordinate x,y.
        
        @param x An int, the x-coordinate to set.
        @param y An int, the y coordinate.
        @param value A double, the value to set.
        """
        self.thisptr.set_imag(x, y, value)

    def get2dMAG(self):
        """!Returns a PyDouble2D object containing the value of the magnitude at each location.
        This method has no direct equivalent in Complex_2D
        @return A PyDouble2D containing the magnitude of the complex value at each point.
        """
        result = PyDouble2D(self.getSizeX(), self.getSizeY()) 
        self.thisptr.get_2d(MAG, (deref(result.thisptr)))
        return result

    def get2dPHASE(self):
        """!Returns a PyDouble2D object containing the value of the phase at each location.
        This method has no direct equivalent in Complex_2D
        @return A PyDouble2D containing the phase part of the value at each point.
        """
        result = PyDouble2D(self.getSizeX(), self.getSizeY()) 
        self.thisptr.get_2d(PHASE, (deref(result.thisptr)))
        return result

    def get2dREAL(self):
        """!Returns a PyDouble2D object containing the value of the real part at each location.
        This method has no direct equivalent in Complex_2D
        @return A PyDouble2D containing the real part of the complex value at each point.
        """
        result = PyDouble2D(self.getSizeX(), self.getSizeY()) 
        self.thisptr.get_2d(REAL, (deref(result.thisptr)))
        return result

    def get2dIMAG(self):
        """!Returns a PyDouble2D object containing the value of the imaginary part at each location.
        This method has no direct equivalent in Complex_2D
        @return A PyDouble2D containing the imaginary part of the complex value at each point.
        """
        result = PyDouble2D(self.getSizeX(), self.getSizeY()) 
        self.thisptr.get_2d(IMAG, (deref(result.thisptr)))
        return result

    def scale(self,scale_factor):
        self.thisptr.scale(scale_factor)
    def add(self, PyComplex2D c2, scale=1):
        self.thisptr.add(deref(c2.thisptr),scale)
    def multiply(self,PyComplex2D c2, scale=1):
        self.thisptr.multiply(deref(c2.thisptr),scale)
    def getNorm(self):
        return self.thisptr.get_norm()
    def invert(self,scale=False):
        self.thisptr.invert(scale)
    def conjugate(self):
        self.thisptr.conjugate()
    def performForwardFFT(self):
        self.thisptr.perform_forward_fft()
    def performBackwardFFT(self):
        """!Do the forward fourier transform on this Complex_2D object
        """
        self.thisptr.perform_backward_fft()
    def performForwardFFTReal(self):
        """!Do the forward fourier transform returning 
        """
        input = PyDouble2D(self.nx,self.ny)
        self.thisptr.perform_forward_fft_real(deref(input.thisptr))
        return input
    def performBackwardFFTReal(self):
        """!Do the backwards fourier transform ...
        """
        result=PyDouble2D(self.nx,self.ny)
        self.thisptr.perform_backward_fft_real(deref(result.thisptr))
        return result
    
    def getPadded(self,x_add,y_add):
        """!Create the same complex with some padding. The padding is filled with 0s.
        @param x_add The number of points to add on each side in horizontal direction
        @param y_add The number of points to add on each side in the vertical direction
        @return A PyComplex2D resized to the padded
        """
        padded=PyComplex2D(self.nx+2*x_add, self.ny+2*y_add)
        cdef Complex_2D * cpad = new Complex_2D(self.thisptr.get_padded(x_add,y_add))
        padded.thisptr=cpad
        return padded
    
    def getUnpadded(self,x_add,y_add):
        """!Return the complex without the padding
        @param x_add The number of padding points to remove from each end in the horizontal direction
        @param y_add The number of padding points to remove from each end in the vertical direction
        """
        unpadded=PyComplex2D(self.nx-2*x_add, self.ny-2*y_add)
        cdef Complex_2D * cupad = new Complex_2D(self.thisptr.get_unpadded(x_add,y_add))
        unpadded.thisptr=cupad
        return unpadded
    
    def read_cplx(self, filename):
        """!Read a cplx file containing an array of complex values into this PyComplex2D object.
        Corresponds to io.read_cplx
        @param filename A string, the name of the file to write to.
        """
        return cio.read_cplx(filename, deref(self.thisptr))

    def write_cplx(self, filename):
        """!Write the array of complex values from this PyComplex2D object to a cplx file.
        Corresponds to io.write_cplx
        @param filename A string, the name of the file to write to.
        """
        return cio.write_cplx(filename, deref(self.thisptr))
        


