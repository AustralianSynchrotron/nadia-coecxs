"""!
This package provides python wrappers to the PhaseDiverseCDI class. 
"""
cimport phasediversecdi
from basecdi cimport BaseCDI, PyBaseCDI
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from phasediversecdi cimport PhaseDiverseCDI
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref
from sets import Set


cdef class PyPhaseDiverseCDI:
    """! PhaseDiverseCDI wrapper.
    
    The PyPhaseDiverseCDI class wraps the C++ PhaseDiverseCDI class 
    and provides a pure python interface to it.
    
    """
       
    def __cinit__(self, beta=1.0, gamma=1.0, parallel=False, granularity=1):
        self.thisptr = new PhaseDiverseCDI(beta, gamma, parallel, granularity)
        self.cdiobjs=set()
        
    def __dealloc__(self):
        del self.thisptr
        
    def addNewPosition(self, local, x=0.0, y=1.0, alpha=1.0):
        """! Add new position
        @param local An instance of PyFresnelCDI or PyPlanarCDI
        @param x 
        @param y
        @param alpha
        """
        basecdi=<PyBaseCDI>local
        self.cdiobjs.add(local)
        self.thisptr.add_new_position(basecdi.thisptr, x, y, alpha)
    def initialiseEstimate(self):
        """! Initialise the reconstruction estimate
        """
        print "initialise estimate"
        self.thisptr.initialise_estimate()
    def iterate(self):
        """! perform one iteration
        """ 
        self.thisptr.iterate()
    def setIterationsPerCycle(self, iterations):
        """! Set number of iterations per cycle
        @param iterations The number of iterations
        """
        self.thisptr.set_iterations_per_cycle(iterations)
    def setFeedbackParameter(self, beta):
        """! Set the Feedback parameter
        @param beta
        """
        self.thisptr.set_feedback_parameter(beta)
    def setAmplificationFactor(self, gamma):
        """! Set the amplification factor 
        @param gamma
        """
        self.thisptr.set_amplification_factor(gamma)
    def setProbeScaling(self, n_probe, alpha):
        """!Set the probing scale
        @param n_probe The probe number
        @param alpha The scaling 
        """
        self.thisptr.set_probe_scaling(n_probe, alpha)
    def setTransmission(self, PyComplex2D trans):
        """! Set the transmission function
        @param trans A PyComplex2D
        """
        self.thisptr.set_transmission(deref(trans.thisptr))
    def adjustPositions(self, type, forwards=True, x_min=-50, x_max=50, y_min=-50, y_max=50, step_size=4.0):
        """!Adjust Positions
        @param type CROSSCORRELATION or MINIMUMERROR
        @param forwards a boolean, which direction to fft default: True
        @param x_min default: -50
        @param x_max default: 50
        @param y_min default: -50
        @param y_max default 50
        @param step_size default: 4.0 
        """
        self.thisptr.adjust_positions(type, forwards, x_min, x_max, y_min, y_max, step_size)
    def getFinalXPosition(self, n_probe):
        """! Get the final x position of the nth probe
        @param n_probe The probe number
        """
        return self.thisptr.get_final_x_position(n_probe)
    def getFinalYPosition(self, n_probe):
        """! get the final y position of the nth probe 
        @param n_probe The probe number
        """
        return self.thispr.get_final_y_position(n_probe)
    def getTransmission(self):
        """! 
        gets the transmission function
        @return A PyComplex2D
        """
        cdef Complex_2D * tmp = self.thisptr.get_transmission().clone()
        nx=tmp.get_size_x()
        ny=tmp.get_size_y()
        trans=PyComplex2D(nx,ny)
        trans.thisptr=tmp
        return trans

            
        
        
CROSSCORRELATION=phasediversecdi.CROSS_CORRELATION
MINIMUMERROR=phasediversecdi.MINIMUM_ERROR