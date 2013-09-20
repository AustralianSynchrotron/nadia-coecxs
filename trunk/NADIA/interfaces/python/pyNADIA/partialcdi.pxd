cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D
from basecdi cimport BaseCDI
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "PartialCDI.h":
     cdef cppclass PartialCDI(BaseCDI):
          PartialCDI()
          PartialCDI(Complex_2D & initial_guess, double lcx, double lcy, \
                                double pxsize,double pysize,double energy, \
                                double zsd, int nleg, int nmode, unsigned int n_best)
          void initialise_estimate()
          void initialise_estimate(int seed)
          int iterate()
          Complex_2D get_transmission()
          void set_transmission(Complex_2D & new_transmission)
          void set_threshold(double new_threshold)
          Double_2D propagate_modes_to_detector()
          Complex_2D get_mode(int mode)
          void set_support(const Double_2D & object_support)
          void set_support(const Double_2D & object_support, bool soften)
        
cdef class PyPartialCDI:
     cdef PartialCDI *thisptr
     cdef int nx
     cdef int ny