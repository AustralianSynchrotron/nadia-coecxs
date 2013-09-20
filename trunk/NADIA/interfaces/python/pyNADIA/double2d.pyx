"""!Python interface to the Double_2D class.
"""

cimport cio as cio
from double2d cimport Double_2D
from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool


cdef class PyDouble2D:
    """!A class wrapping the Double_2D python class. 
    
    PyDouble2D contains a pointer to a Double_2D class and provides a pure python interface
    to the public methods of the C class. The read/write io functions which operate on Double_2D have been incorporated as methods
    into this class. See Double_2D documentation for further details on methods.
    
    """
    
    def __cinit__(self, int x_size=0, int y_size=0):
        """!The constructor."""
        if (x_size == 0 and y_size == 0):
            self.thisptr = new Double_2D()
        else:
            self.thisptr = new Double_2D(x_size, y_size)
        if self.thisptr == NULL:
            raise MemoryError('Not enough memory.')
    
    def __dealloc__(self):
        del self.thisptr
    
    def allocateMemory(self,x_size,y_size):
        """! Allocate memory for the underlying Double_2D
        @param x_size The size in the horizontal direction
        @param y_size The size in the vertical direction
        """
        self.thisptr.allocate_memory(x_size,y_size)
    
    def getMax(self):
        """!Returns the maximum value in the 2D array.
        Corresponds to Double_2D.get_max()
        @return The maximum value in the array.
        """
        return self.thisptr.get_max()
    
    def getMin(self):
        """!Returns the minimum value in the 2D array.
        Corresponds to Double_2D.get_min()
        @return The minimum value in the array.
        """
        return self.thisptr.get_min()
    
    def getSum(self):
        """!Returns the sum of all values in the 2D array.
        Corresponds to Double_2D.get_sum()
        @return The sum of the values in the array.
        """
        return self.thisptr.get_sum()
    
    def getAbsSum(self):
        """!Returns the sum of the absolute values of the 2D array elements.
        Corresponds to Double_2D.get_abs_sum()
        @return The sum of the absolute values in the array.
        """
        return self.thisptr.get_abs_sum()
    
    def getShape(self):
        """!Returns a tuple containing the dimensions of the array (nx,ny).
        
        @return a tuple, with the (x,y) sizes of the array (int,int)
        """
        return (self.thisptr.get_size_x(), self.thisptr.get_size_y())
    
    def get(self, i, j):
        """!Returns the value at position i,j.
        
        @param i (int) The position in the horizontal direction
        @param j (int) The position in the vertical direction
        """
        return self.thisptr.get(i, j)
    
    def set(self, i, j,value):
        """!Returns the value at position i,j.
        
        @param i (int) The position in the horizontal direction
        @param j (int) The position in the vertical direction
        @param value (real) The value to set
        """
        self.thisptr.set(i, j,value)
    
    def add(self, PyDouble2D other, norm=1.0):
        """!Add another PyDouble2D to this one
        
        @param other (PyDouble2D) The other PyDouble2D to add to this one
        @param norm (double) Normalise other array by this factor
        """
        self.thisptr.add(deref(other.thisptr),norm)
        
    def scale(self, scale_factor):
        """!Scale all items in the array in-place
        @param scale_factor
        """
        self.thisptr.scale(scale_factor)
    def sqrt(self):
        """!Take the square root of all items in the array in-place
        """
        self.thisptr.sq_root()
    def square(self):
        """!Take the square of all items in the array in-place
        """
        self.thisptr.square()
   # def getArray(self):
   #     """! Return the 2D array as a numpy array.
        
    #    @return Returns a numpy array of shape (nx,ny) containing the data.
    #    """
    #    shape = self.getShape()
    #    arr = np.empty(shape)
    #    for i in range(0, shape[0]):
    #        for j in range(0, shape[1]):
     #           arr[i][j] = self.thisptr.get(i, j)
   #     cdef * 
    #    return arr
    
    def read_tiff(self, filename):
        """! Read data from a tiff file into a PyDouble2D object.
        
        @param filename The name of the tiff file to read from.
        """
        return cio.read_tiff(filename, deref(self.thisptr))
    
    def read_dbin(self, filename, nx, ny):
        """! Read data from a dbin file into a PyDouble2D object of size nx,ny.
        
        @param filename The name of the file to read from
        @param nx (int) Number of samplings in the horizontal direction
        @param ny (int) Number of samplings in the vertical direction
        """
        return cio.read_dbin(filename, nx, ny, deref(self.thisptr))
    
    def read_spec(self, filename):
        """! Read data from a spec file into a PyDouble2D object.
        @param filename A string, the name of the file to read from.
        """
        return cio.read_spec(filename, deref(self.thisptr))
    
    def read_ppm(self, filename):
        """! Read data from a ppm file into a PyDouble2D object.
        
        @param filename A string, the name of the file to read from.
        """
        return cio.read_ppm(filename, deref(self.thisptr))
    
    def read_hdf4(self, filename):
        """! Read data from a hdf4 file into a PyDouble2D object.
        
        @param filename A string, the name of the file to read from
        """
        return cio.read_hdf4(filename, deref(self.thisptr))
    
    def write_ppm(self, filename, log_scale=False, min=0, max=0):
        """! Write data from the PyDouble2D object to a ppm file.
        
        @param filename  The name of the file to write to
        @param log_scale A bool, write values on a log scale.
        @param min
        @param max
        """
        if min == 0 and max == 0:
            if log_scale:
                return cio.write_ppm(filename, deref(self.thisptr), log_scale)
            else:
                return cio.write_ppm(filename, deref(self.thisptr))
        else:
            return cio.write_ppm(filename, deref(self.thisptr), log_scale, min, max)

    def write_tiff(self, filename, log_scale=False, min=0, max=0):
        """! Write data from the PyDouble2D object to a tiff file.
        
        @param filename  A string, the name of the file to write to
        @param log_scale A boolean, display on a log scale.
        @param min A double.
        @param max A double.
        """
        if min == 0 and max == 0:
            if log_scale:
                return cio.write_tiff(filename, deref(self.thisptr), log_scale)
            else:
                return cio.write_tiff(filename, deref(self.thisptr))
        else:
            return cio.write_tiff(filename, deref(self.thisptr), log_scale, min, max)
        
    def write_dbin(self, filename):
        """! Write data from the PyDouble2D object to a dbin file.
        
        @param filename A string, the file name to write to.
        """
        return cio.write_dbin(filename, deref(self.thisptr))

    def write_image(self,file_name, log_scale=False, min=0, max=0):
        """! Write data from the PyDouble2D object to an image file.
        @param file_name
        @param log_scale
        @param min default:0
        @param max default:0
        """
        if min == 0 and max == 0:
            if log_scale:
                cio.write_image(file_name, deref(self.thisptr), log_scale)
            else:
                cio.write_image(file_name, deref(self.thisptr))
        else:
            cio.write_image(file_name, deref(self.thisptr), log_scale, min, max)
    
    def getSizeX(self):
        """! Return size of the PyDouble2D array in the x dimension.
        
        @return The number of horizontal points
        """

        return self.thisptr.get_size_x()

    def getSizeY(self):
        """! Return size of the PyDouble2D array in the y dimension.
        
        @return The number of positions in the vertical direction
        """
        return self.thisptr.get_size_y()
