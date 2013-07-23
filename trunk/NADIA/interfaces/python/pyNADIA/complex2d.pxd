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
    cdef cppclass Complex_2D:
        #Complex_2D() except+
        Complex_2D(int,int) except+
        Complex_2D(const Complex_2D& )
        #Complex_2D& operator=(const Complex_2D& rhs)
        void set_value(int x, int y, int type, double value)
        void set_real(int x, int y, double value)
        void set_imag(int x, int y, double value)
        double get_imag(int x, int y)
        double get_mag(int x, int y)
        double get_phase(int x, int y)
        double get_value(int x, int y, int type)
        int get_size_x()
        int get_size_y()
        void get_2d(int type, Double_2D & result)
        void scale(double scale_factor)
        void add(Complex_2D & c2)
        void add(Complex_2D & c2, double scale)
        void multiply(Complex_2D & c2)
        void multiply(Complex_2D & c2, double scale)
        double get_norm()
        Complex_2D * clone()
        void copy(const Complex_2D & c)
        void invert()
        void invert(bool scale)
        void conjugate()
        void perform_forward_fft()
        void perform_backward_fft()
        void perform_backward_fft_real(Double_2D & result)
        void perform_forward_fft_real(Double_2D & input)
        void set_fftw_type(int type)
        Complex_2D get_padded(int x_add, int y_add)
        Complex_2D get_unpadded(int x_add, int y_add)
    
    
cdef class PyComplex2D:
    cdef Complex_2D *thisptr
    