"""!
This package provides python wrappers to the FresnelCDI class.

More package details.
"""
from fresnelcdi cimport FresnelCDI
from double2d cimport Double_2D, PyDouble2D
from transmissionconstraint cimport PyTransmissionConstraint
from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool
from complex2d cimport Complex_2D, PyComplex2D


cdef class PyFresnelCDI:
    """! PyFresnelCDI class.
    
    A python wrapper for the FresnelCDI class.
    """

    def __cinit__(self, PyComplex2D object_estimate, PyComplex2D wf, double wavelength, double ftd, double fts, double pixel_size, double norm):
        self.thisptr = new FresnelCDI(deref(object_estimate.thisptr), deref(wf.thisptr), wavelength, ftd, fts, pixel_size, norm)
        self.nx = object_estimate.thisptr.get_size_x()
        self.ny = object_estimate.thisptr.get_size_y()
    
    def __dealloc__(self):
        del self.thisptr
    
    def iterate(self):
        """!Perform one iteration."""
        self.thisptr.iterate()

    def setSupport(self, PyDouble2D support):
        """!Set the support.
        @param support A PyDouble2D object containing the support.
        """
        cdef Double_2D csupport = deref(support.thisptr)
        self.thisptr.set_support(csupport, False)
    
    def getTransmissionFunction(self):
        """! Returns the transmission function.
        @returns a PyComplex2D containing the transmission function.
        """
        trans = PyComplex2D(self.nx, self.ny) 
        self.thisptr.get_transmission_function(deref(trans.thisptr))
        return trans
    
    def applyShrinkwrap(self, x, y):
        """! Apply the shrinkwrap algorithm.
        @param x 
        @param y
        """
        self.thisptr.apply_shrinkwrap(x, y)
    
    def getIlluminationAtSample(self):  
        """!Return the illumination at the sample plane.
        @return A PyComplex2D containing the illumination at the sample plane.
        """
        result = PyComplex2D(self.nx, self.ny)
        result.thisptr = self.thisptr.get_illumination_at_sample().clone()
        return result 
    
    def getError(self):
        """!Returns the current error 
        @return A double, the current value of the error.
        """
        return self.thisptr.get_error()
    
    def initialiseEstimate(self, n=0):
        """! Initialise the object esw at the sample plane.
        @param n An optional parameter, the seed of the random number generator 
        """
        self.thisptr.initialise_estimate(n)
    
    def setIntensity(self, PyDouble2D intensity):
        """!Set the intensity
        @param intensity A PyDouble2D containing the intensity 
        """
        self.thisptr.set_intensity(deref(intensity.thisptr))
    
    def setAlgorithm(self, algorithm):
        """!Set the algorithm used for the iterations.
        @param algorithm A string
        """
        self.thisptr.set_algorithm(self.thisptr.getAlgFromName(algorithm))
    
    def getBestResult(self, error, index=0):
        """! Get the current best reconstructed image
        @param error 
        @param index 
        @return a PyComplex2D containing the best reconstruction data
        """
        result = PyComplex2D(self.nx, self.ny) 
        result.thisptr = self.thisptr.get_best_result(error, index)
        return result
    
    def setBeamStop(self, PyDouble2D beamstop):
        """! Set the beam stop
        @param beamstop A PyDouble2D representing the beamstop data.
        """
        self.thisptr.set_beam_stop(deref(beamstop.thisptr))
    
    def getSupport(self):
        """! Get the current support.
        @return a PyDouble2D containing the support
        """
        support = PyDouble2D()
        cdef Double_2D * tmp = new Double_2D(self.thisptr.get_support())
        support.thisptr = tmp
        return support
    
    def resetBest(self):
        """! Reset the best reconstruction results 
        """
        self.thisptr.reset_best()
    
    def set_fftw_type(self, type):
        """! Set the FFTW type.
        
        @param type FFTW_ESTIMATE, FFTW_MEASURE or FFTW_PATIENT
        """
        self.thisptr.set_fftw_type(type)
    
    def autoSetNorm(self):
        self.thisptr.auto_set_norm()
        
    def applySupport(self,PyComplex2D c):
        """!Apply the support constraint)
        @param c A PyComplex2D to apply the support constraint to.
        @return The complex field with support constraint applied (PyComplex2D)
        """
        self.thisptr.apply_support(deref(c.thisptr))
        return c
    def setTransmissionFunction(self,PyComplex2D transmission, PyComplex2D esw=None):
        if esw is None:
            self.thisptr.set_transmission_function(deref(transmission.thisptr))
        else:
            esw=PyComplex2D(self.nx,self.ny)
            self.thisptr.set_transmission_function(deref(transmission.thisptr),esw.thisptr)

    def scaleIntensity(self,PyComplex2D c):
        """!Intensity is scaled to match the data, error also is updated.
        @param c (PyComplex2D) A complex field to apply the intensity scaling on.
        @return The scaled complex field (PyComplex2D)
        """
        self.thisptr.scale_intensity(deref(c.thisptr))
        return c
    def propagateToDetector(self,PyComplex2D c):
        """!Propagate a complex field to the detector plane
        @param c (PyComplex2D) A complex field to propagate to the detector plane using an fft
        @return the transformed complex field (PyComplex2D)
        """
        self.thisptr.propagate_to_detector(deref(c.thisptr))
    def propagateFromDetector(self,PyComplex2D c):
        """!Propagate a complex field to the sample plane
        @param c (PyComplex2D) A complex field to propagate to the sample plane using an fft
        @return the transformed complex field (PyComplex2D)
        """
        self.thisptr.propagate_from_detector(deref(c.thisptr))
        
    def setNormalisation(self,normalisation):
        self.thisptr.set_normalisation(normalisation)
    def setExperimentalParameters(self, beam_wavelength, focal_detector_length,focal_sample_length,pixel_size):
        self.thisptr.set_experimental_parameters(beam_wavelength, focal_detector_length, focal_sample_length, pixel_size)
        
    def setNorm(self, new_normalisation):
        self.thisptr.set_norm(new_normalisation)
    def setComplexConstraint(self, PyTransmissionConstraint tc):
        self.thisptr.set_complex_constraint(deref(tc.thisptr))
 