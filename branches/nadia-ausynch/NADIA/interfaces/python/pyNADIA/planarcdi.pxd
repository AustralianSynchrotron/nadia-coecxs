cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D
from basecdi cimport BaseCDI
from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "PlanarCDI.h":
     cdef cppclass PlanarCDI(BaseCDI):
          PlanarCDI(Complex_2D & complex,unsigned int n_best)
          void get_intensity_autocorrelation(Double_2D & autoc)
          void initialise_estimate(int seed)
          void propagate_to_detector(Complex_2D & c)
          void propagate_from_detector(Complex_2D & c)
          void scale_intensity(Complex_2D & c)
          int iterate()
          void apply_support(Complex_2D & c)
          void set_support(const Double_2D & object_support)
          void set_support(const Double_2D & object_support, bool soften)
cdef class PyPlanarCDI:
     cdef PlanarCDI *thisptr
     cdef int nx,ny