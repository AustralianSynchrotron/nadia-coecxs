"""!
Package for wrapping the PartialCDI class.
"""
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from partialcdi cimport PartialCDI
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint

from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref


cdef class PyPartialCDI:
    """! PyPartialCDI class.
    PyPartialCDI is the python wrapper for the PartialCDI class. In general there is a Python method corresponding to each public method
    in the C++ class. This documentation should detail the parameters and return values for the python methods. 
    For more details on the purpose of each function see the documentation in PartialCDI (or BaseCDI for inherited methods).
    """
    def __cinit__(self, PyComplex2D object_estimate, double lcx, double lcy, double pxsize, double pysize, double energy, double zsd, int nleg, int nmode, unsigned int n_best):
        self.thisptr = new PartialCDI(deref(object_estimate.thisptr), lcx, lcy, pxsize, pysize, energy, zsd, nleg, nmode, n_best)
        self.nx = object_estimate.thisptr.get_size_x()
        self.ny = object_estimate.thisptr.get_size_y()
    def __dealloc__(self):
        del self.thisptr
    
    def initialiseEstimate(self, seed=0):
        """!Initialise the estimated diffraction pattern.
        Calls PartialCDI.initialise_estimate()
        @param seed (int) An optional seed parameter, default=0
        """
        self.thisptr.initialise_estimate(seed)
    def iterate(self):
        """!Perform one iteration, as implemented in PartialCDI.iterate()

        """
        return self.thisptr.iterate()
    def getTransmission(self):
        """! Get the 'global' sample function.
        Corresponds to PartialCDI.get_transmission()
        @return A PyComplex2D the 'global' sample function.
        """
        trans = PyComplex2D(self.nx, self.ny)
        cdef Complex_2D * tmp = new Complex_2D(self.thisptr.get_transmission())
        trans.thisptr = tmp
        return trans
    def setTransmission(self, PyComplex2D newTrans):
        """! Set the 'global' sample function.
        
        Useful for starting a reconstruction using results of a previous one. Calls PartialCDI.set_transmission()
        @param newTrans A PyComplex2D the new 'global' sample function.
        """
        self.thisptr.set_transmission(deref(newTrans.thisptr))
    def setThreshold(self, double new_threshold):
        """!Set the minimum contribution of a mode as a proportion of the dominant mode for it to be included in the reconstruction.
        
        Corresponds to PartialCDI.set_threshold()
        @param new_threshold (double)
        """
        self.thisptr.set_threshold(new_threshold)
    def getMode(self, int modenum):
        """! Returns a given mode.
        @param modenum (int) The mode number.
        @return the mode (PyComplex2D)
        """
        mode = PyComplex2D(self.nx, self.ny)
        cdef Complex_2D * tmp = new Complex_2D(self.thisptr.get_mode(modenum))
        mode.thisptr = tmp
        return mode
    def getBestResult(self, double & error, int index=0):
        """! Get the current best result of the reconstruction.
        @return the best result (PyComplex2D).
        """
        result = PyComplex2D(self.nx,self.ny)
        result.thisptr=self.thisptr.get_best_result(error,index)
        return result
    def setSupport(self, PyDouble2D support, bool soften=False):
        """! Set the support
        @param support (PyDouble2D) The support function
        @param soften (bool) Convolve the edge of the support with a 3 pixel wide Gaussian. Default=False
        """
        self.thisptr.set_support(deref(support.thisptr), soften)
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
    def setAlgorithm(self, int alg):
        """!Set the algorithm. 
        By default HIO is used. See pyNADIA.utils for the the constants used.
        @param arg (int) The algorithm
        """
        self.thisptr.set_algorithm(alg)
    def resetBest(self):
        self.thisptr.reset_best()
    def set_fftw_type(self, int type):
        self.thisptr.set_fftw_type(type)
