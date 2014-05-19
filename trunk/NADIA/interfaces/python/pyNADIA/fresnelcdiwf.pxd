cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D
from basecdi cimport BaseCDI
from libcpp.string cimport string

cdef extern from "FresnelCDI_WF.h":
    cdef cppclass FresnelCDI_WF(BaseCDI):
        FresnelCDI_WF(Complex_2D & complex, double beam_wavelength, \
                        double zone_focal_length, double focal_detector_length,\
                        double pixel_size) except+
        FresnelCDI_WF(Complex_2D & complex, double beam_wavelength, \
                        double zone_focal_length, double focal_detector_length,\
                        double pixel_size, unsigned int n_best) except+
        void initialise_estimate(int seed)
        void propagate_to_detector(Complex_2D & c)
        void propagate_from_detector(Complex_2D & c)
        void scale_intensity(Complex_2D & c)
        int iterate()
        void set_support(double z_factor)
        void set_support(double z_factor, double size)
        void multiply_factors(Complex_2D & c, int direction)

cdef class PyFresnelCDIWF:
    cdef FresnelCDI_WF *thisptr
    cdef int nx,ny