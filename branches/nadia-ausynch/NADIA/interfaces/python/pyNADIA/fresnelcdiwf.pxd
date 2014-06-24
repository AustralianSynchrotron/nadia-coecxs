cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D
from basecdi cimport BaseCDI
from libcpp.string cimport string

cdef extern from "FresnelCDI_WF.h":
    cdef cppclass FresnelCDI_WF(BaseCDI):
        FresnelCDI_WF(Complex_2D & c, double beam_wavelength, \
                        double zone_focal_length, double focal_detector_length,\
                        double pixel_size) except+
        FresnelCDI_WF(Complex_2D & c, double beam_wavelength, \
                        double zone_focal_length, double focal_detector_length,\
                        double pixel_size, int n_best) except+
        void initialise_estimate(int seed) except+
        void propagate_to_detector(Complex_2D & c)except+
        void propagate_from_detector(Complex_2D & c) except+
        void scale_intensity(Complex_2D & c) except+
        int iterate() 
        void set_support(double z_factor) except+
        void set_support(double z_factor, double size) except+
        void multiply_factors(Complex_2D & c, int direction) except+

cdef class PyFresnelCDIWF:
    cdef FresnelCDI_WF *thisptr
    cdef object __weakref__
    cdef int nx,ny