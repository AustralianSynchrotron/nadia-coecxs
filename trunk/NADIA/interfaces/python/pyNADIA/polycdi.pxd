cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D
from basecdi cimport BaseCDI
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "PolyCDI.h":
     cdef cppclass PolyCDI(BaseCDI):
          PolyCDI(Complex_2D & initial_guess)
          PolyCDI(Complex_2D & initial_guess,double beta, unsigned int n_best, bool parallel)
          void initialise_estimate()
          void initialise_estimate(int seed)
          void scale_intensity(Complex_2D & c)
          void expand_wl(Complex_2D & c)
          Double_2D sum_intensity(vector[Complex_2D] & c)
          int iterate()
          void update_transmission()
          void set_iterations_per_cycle(int iterations)
          Complex_2D get_transmission()
          void set_transmission(Complex_2D & new_transmission)
          void set_spectrum(Double_2D spectrum_in)
          void set_spectrum(string file_name)
          Double_2D get_intensity()
          Double_2D propagate_modes_to_detector()
          Complex_2D get_mode(int mode)
          void propagate_to_detector(Complex_2D & c)
          void propagate_from_detector(Complex_2D & c)
        
cdef class PyPolyCDI:
     cdef PolyCDI *thisptr
     cdef int nx,ny