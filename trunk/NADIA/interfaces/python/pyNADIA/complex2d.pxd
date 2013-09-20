cimport cython
from double2d cimport Double_2D
from libcpp cimport bool 

cdef extern from "Complex_2D.h":
    cdef int MAG "MAG"
    cdef int PHASE "PHASE"
    cdef int SUCCESS "SUCCESS"
    cdef int FAILURE "FAILURE"
    cdef int REAL "REAL"
    cdef int IMAG "IMAG"
    cdef int MAG_SQ "MAG_SQ"
    cdef cppclass ComplexR_2D[T]:
        #Complex_2D() except+
        ComplexR_2D(int,int) except+
        ComplexR_2D(const ComplexR_2D[T]& )
        #Complex_2D& operator=(const Complex_2D& rhs)
        void set_value(int x, int y, int type, T value)
        void set_real(int x, int y, T value)
        void set_imag(int x, int y, T value)
        T get_real(int x,int y)
        T get_imag(int x, int y)
        T get_mag(int x, int y)
        T get_phase(int x, int y)
        T get_value(int x, int y, int type)
        int get_size_x()
        int get_size_y()
        void get_2d(int type, Double_2D & result)
        void scale(T scale_factor)
        void add(ComplexR_2D[T] & c2)
        void add(ComplexR_2D[T] & c2, T scale)
        void multiply(ComplexR_2D[T] & c2)
        void multiply(ComplexR_2D[T] & c2, T scale)
        T get_norm()
        ComplexR_2D * clone()
        void copy(const ComplexR_2D[T] & c)
        void invert()
        void invert(bool scale)
        void conjugate()
        void perform_forward_fft()
        void perform_backward_fft()
        void perform_backward_fft_real(Double_2D & result)
        void perform_forward_fft_real(Double_2D & input)
        void set_fftw_type(int type)
        ComplexR_2D get_padded(int x_add, int y_add)
        ComplexR_2D get_unpadded(int x_add, int y_add)
IF DOUBLE_PRECISION !='1':
    ctypedef ComplexR_2D[float] Complex_2D
ELSE:
    ctypedef ComplexR_2D[double] Complex_2D
        
cdef class PyComplex2D:
    cdef Complex_2D *thisptr
    cdef object __weakref__
    