"""!@package pyNADIA.planarcdi 
Python interface to the PlanarCDI class.

"""
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from planarcdi cimport PlanarCDI
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint
from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef class PyPlanarCDI:
    """!PyPlanarCDI is a wrapper class for the PlanarCDI class. 
    Wrapper around the PlanarCDI class and provides Python only access to its methods.
    """
    def __cinit__(self,PyComplex2D object_estimate, int n_best):
        self.thisptr=new PlanarCDI(deref(object_estimate.thisptr),n_best)
        self.nx=object_estimate.thisptr.get_size_x()
        self.ny=object_estimate.thisptr.get_size_y()
      
    def __dealloc__(self):
        del self.thisptr
    def getIntensityAutocorrelation(self):
        """!Returns the autocorrelation of the data in a PyDouble2D object.
           
        @return The autocorrelation as a PyDouble2D  
           
        """
        autoc = PyDouble2D(self.nx,self.ny)
        self.thisptr.get_intensity_autocorrelation(deref(autoc.thisptr))
        return autoc
    
    # all these are inherited from the baseCDI class
    def iterate(self):
        """!Performs one iteration of propagating the image to and from the detector plane.
        
        Performs one iteration of propagating the image to and from the detector plane.
           
        """
        self.thisptr.iterate()
         
    def setRelaxationParameter(self, relax_parameter):
        """!Sets the relaxation parameter.
        
        @param relax_parameter The relaxation parameter
             
        """
        self.thisptr.set_relaxation_parameter(relax_parameter)
         
    def getBestResult(self, error, index=0):
        """!Returns a PyComplex2D holding the current best result.
        
        @param error 
        @param index Index of the best result
           
        @return  A PyComplex2D containing the current best result.
        """
        result = PyComplex2D(self.nx,self.ny) 
        result.thisptr=self.thisptr.get_best_result(error,index)
        return result
    def initialiseEstimate(self, seed=0):
        """! Initialises the reconstructed image, with an optional random seed.
            
        Initialises the reconstructed image, with an optional random seed.
         
        @param seed 
        """
        self.thisptr.initialise_estimate(seed)
    def setSupport(self, PyDouble2D support, soften=False):
        """!Set the support.
        @param support A PyDouble2D object containing the support.
        """
        self.thisptr.set_support(deref(support.thisptr),soften)
    def setIntensity(self,PyDouble2D intensity):
        """!Set the intensity
        @param intensity A PyDouble2D containing the intensity 
        """
        self.thisptr.set_intensity(deref(intensity.thisptr))
    def setBeamStop(self,PyDouble2D beamstop):
        """! Set the beam stop
        @param beamstop A PyDouble2D representing the beamstop data.
        """
        self.thisptr.set_beam_stop(deref(beamstop.thisptr))
    def getSizeX(self):
        """!Get the number of pixels in the x co ordinate.
        
        @return An int, the number of pixels in the x coordinate.
        """
        return self.nx
    def getSizeY(self):
        """!Get the number of pixels in the y co ordinate.
        
        @return An int, the number of pixels in the y coordinate.
        """
        return self.ny
    def getExitSurfaceWave(self):
        """!Get the exit surface wave
        @returns The exit surface wave (PyComplex2D)
        """
        esw = PyComplex2D(self.nx,self.ny)
        cdef Complex_2D * tmp = new Complex_2D(self.thisptr.get_exit_surface_wave())
        esw.thisptr=tmp
        return esw
    def setAlgorithm(self, algorithm):
        """!Set the algorithm to use
         @param algorithm A string
        """
        self.thisptr.set_algorithm(self.thisptr.getAlgFromName(algorithm))
    def setCustomAlgorithm(self, ms):
        """! Set a custom algorithm
        @param ms A list of 10 doubles to specify a custom algorithm. See
        set_custom_algorithm() documentation for details
        """
        if ms.length()!=10:
            print "argument must be an array of 10 doubles"
            exit 
        self.thisptr.set_custom_algorithm(ms[0],ms[1],ms[2],ms[3],ms[4], ms[5],ms[6],ms[7],ms[8],ms[9])
    def printAlgorithm(self):
        """!Print the algorith used
        """
        self.thisptr.print_algorithm()
    def getError(self):
        """!Returns the current error 
        @return The current value of the error metric (double).
        """
        return self.thisptr.get_error()
    def applyShrinkwrap(self,gauss_width=1.5,threshold=0.1):
        """! Apply shrinkwrap algorithm.
        
        @param gauss_width (default=1.5)
        @param threshold (default=0.1)
        
        """
        self.thisptr.apply_shrinkwrap(gauss_width,threshold)
    def getSupport(self):
        """! Get the current support.
        @return a PyDouble2D containing the support
        """
        support = PyDouble2D()
        cdef Double_2D * tmp = new Double_2D(self.thisptr.get_support())
        support.thisptr=tmp
        return support
    def applySupport(self,PyComplex2D c):
        """!Apply the support constraint)
        @param c A PyComplex2D to apply the support constraint to.
        @return The complex field with support constraint applied (PyComplex2D)
        """
        self.thisptr.apply_support(deref(c.thisptr))
        return c
    def projectIntensity(self, PyComplex2D c):
        """!Apply the intensity constraint
        @param c A PyComplex2D to apply to the intensity constraint on
        @return A complex field with intensity constraint applied. (PyComplex2D) 
        """
        self.thisptr.project_intensity(deref(c.thisptr))
        return c
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
    def setFFTWType(self,type):
        """Set the fft type
        @param type FFTW_ESTIMATE, FFTW_MEASURE or FFTW_PATIENT
        """
        self.thisptr.set_fftw_type(type)
    def setComplexConstraint(self,PyTransmissionConstraint ptc):
        """! Sets a complex Constraint.
        @param PyTransmissionConstraint ptc
        """
        self.thisptr.set_complex_constraint(deref(ptc.thisptr))


