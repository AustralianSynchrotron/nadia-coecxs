"""!
Package for wrapping the FresnelCDI_WF class.
"""
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from fresnelcdiwf cimport FresnelCDI_WF
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint

from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref


cdef class PyFresnelCDIWF:
    """! PyFresnelCDIWF class.
    PyFresnelCDIWF is the python wrapper for the FresenelCDI_WF class. In general there is a Python method corresponding to each public method
    in the C++ class. This documentation should detail the parameters and return values for the python methods. 
    For more details on the purpose of each function see the documentation in FresnelCDI_WF (or BaseCDI for inherited methods).
    """
    def __cinit__(self, PyComplex2D object_estimate, double beam_wavelength, double zone_focal_length, double focal_detector_length, \
                    double pixel_size,int n_best=1):
        self.thisptr = new FresnelCDI_WF(deref(object_estimate.thisptr), beam_wavelength, zone_focal_length, focal_detector_length,\
                                         pixel_size, n_best)
        self.nx = object_estimate.thisptr.get_size_x()
        self.ny = object_estimate.thisptr.get_size_y()
    def __dealloc__(self):
        del self.thisptr
    
    def initialiseEstimate(self, seed=0):
        """!Initialise the estimated diffraction pattern.
        Calls FresnelCDI_WF.initialise_estimate()
        @param seed (int) An optional seed parameter, default=0
        """
        self.thisptr.initialise_estimate(seed)
    def iterate(self):
        """!Perform one iteration, as implemented in PartialCDI.iterate()

        """
        return self.thisptr.iterate()
    def getBestResult(self, double & error, int index=0):
        """! Get the current best result of the reconstruction.
        @return the best result (PyComplex2D).
        """
        result = PyComplex2D(self.nx,self.ny)
        result.thisptr=self.thisptr.get_best_result(error,index)
        return result
    def setSupport(self, z_factor,size=1.01):
        """! Set the support
        @param z_factor (double)
        @param size (double) default=1.01
        """
        self.thisptr.set_support(z_factor,size)
    def setIntensity(self, PyDouble2D intensity):
        """! Set the detector diffraction image (the square of the amplitude of the wavefield at the detector).
        @param intensity (PyDouble2D) The intensity at the detector
        """
        self.thisptr.set_intensity(deref(intensity.thisptr))
    def setBeamStop(self, PyDouble2D beamstop):
        """!Set the beam-stop position in the detector plane.        
        @param beamstop (PyDouble2D) The region of the beam stop in the detector plane.
        """
        self.thisptr.set_beam_stop(deref(beamstop.thisptr))
    def getError(self):
        """!Get the error.
        Returns the difference between the estimated diffraction adn the actual diffraction pattern.
        @return the error metric (double)
        """
        return self.thisptr.get_error()
    def getSupport(self):
        """! Get the support
        @return The object support. (PyDouble2D)
        """
        result = PyDouble2D()
    def resetBest(self):
        self.thisptr.reset_best()
    def set_fftw_type(self, int type):
        self.thisptr.set_fftw_type(type)
    def multiplyFactors(self, PyComplex2D c, direction):
        self.thisptr.multiply_factors(deref(c.thisptr), direction)
  