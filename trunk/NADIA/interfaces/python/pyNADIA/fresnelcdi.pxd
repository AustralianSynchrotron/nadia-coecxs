cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D
from basecdi cimport BaseCDI
from libcpp.string cimport string

cdef extern from "FresnelCDI.h":
     cdef cppclass FresnelCDI(BaseCDI):
         #FresnelCDI() except+
         FresnelCDI(Complex_2D & initial_guess) except+
         FresnelCDI(Complex_2D & initial_guess, int n_best) except+
         FresnelCDI(Complex_2D & initial_guess,\
	     Complex_2D & white_field,\
	     double beam_wavelength, double focal_detector_length,\
	     double focal_sample_length, double pixel_size)
         FresnelCDI(Complex_2D & initial_guess,\
	     Complex_2D & white_field,\
	     double beam_wavelength,double focal_detector_length,\
	     double focal_sample_length, double pixel_size, double normalisation)
         FresnelCDI(Complex_2D & initial_guess,\
	     Complex_2D & white_field,\
	     double beam_wavelength, double focal_detector_length,\
	     double focal_sample_length, double pixel_size, double normalisation,int n_best)
         void initialise_estimate() 
         void initialise_estimate(int seed)
         void auto_set_norm()
         int iterate()
         void get_transmission_function(Complex_2D & result)
         void get_transmission_function(Complex_2D & result, Complex_2D * esw)
         void apply_support(Complex_2D & complex)
         void set_transmission_function(Complex_2D & transmission)
         void set_transmission_function(Complex_2D & transmission, Complex_2D *esw)
         void scale_intensity(Complex_2D & c)
         void propagate_from_detector(Complex_2D & c)
         void propagate_to_detector(Complex_2D & c)
         void set_normalisation(double normalisation)
         void set_experimental_parameters(double beam_wavelength, \
              double focal_detector_length, double focal_sample_length, double pixel_size)
         void set_norm(double new_normalisation)
         const Complex_2D & get_illumination_at_sample()
         Complex_2D illumination
         Complex_2D * illumination_at_sample
         Complex_2D * transmission
         
cdef class PyFresnelCDI:
     cdef FresnelCDI *thisptr            
     cdef int nx,ny