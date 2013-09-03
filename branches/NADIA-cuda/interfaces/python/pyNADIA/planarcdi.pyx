cimport planarcdi
cimport double2d
cimport complex2d
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from planarcdi cimport PlanarCDI, PyPlanarCDI
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint

from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef class PyPlanarCDI:
     def __cinit__(self,PyComplex2D object_estimate, int n_best):
         self.thisptr=new PlanarCDI(deref(object_estimate.thisptr),n_best)
         self.nx=object_estimate.thisptr.get_size_x()
         self.ny=object_estimate.thisptr.get_size_y()
     def getIntensityAutocorrelation(self):
         autoc = PyDouble2D(self.nx,self.ny)
         self.thisptr.get_intensity_autocorrelation(deref(autoc.thisptr))
         return autoc
     # all these are inherited from the baseCDI class
     def iterate(self):
         self.thisptr.iterate()
     def setRelaxationParameter(self, double param):
         self.thisptr.set_relaxation_parameter(param)
     def getBestResult(self,double error, int index=0):
         result = PyComplex2D(self.nx,self.ny) 
         result.thisptr=self.thisptr.get_best_result(error,index)
         return result
     def initialiseEstimate(self,int seed=0):
         self.thisptr.initialise_estimate(seed)
     def setSupport(self,PyDouble2D support, bool soften=False):
         self.thisptr.set_support(deref(support.thisptr),soften)
     def setIntensity(self,PyDouble2D intensity):
         self.thisptr.set_intensity(deref(intensity.thisptr))
     def setBeamStop(self,PyDouble2D beamstop):
         self.thisptr.set_beam_stop(deref(beamstop.thisptr))
     def getSizeX(self):
         return self.nx
     def getSizeY(self):
         return self.ny
     def getExitSurfaceWave(self):
         esw = PyComplex2D(self.nx,self.ny)
         cdef Complex_2D * tmp = new Complex_2D(self.thisptr.get_exit_surface_wave())
         esw.thisptr=tmp
         return esw
     def setAlgorithm(self,string algorithm):
         self.thisptr.set_algorithm(self.thisptr.getAlgFromName(algorithm))
     def setCustomAlgorithm(self, ms):
         if ms.length()!=10:
            print "argument must be an array of 10 doubles"
            exit 
         self.thisptr.set_custom_algorithm(ms[0],ms[1],ms[2],ms[3],ms[4],\
                ms[5],ms[6],ms[7],ms[8],ms[9])
     def printAlgorithm(self):
         self.thisptr.print_algorithm()
     def getError(self):
         return self.thisptr.get_error()
     def applyShrinkwrap(self,double gauss_width=1.5,double threshold=0.1):
         self.thisptr.apply_shrinkwrap(gauss_width,threshold)
     def getSupport(self):
         support = PyDouble2D()
         cdef Double_2D * tmp = new Double_2D(self.thisptr.get_support())
         support.thisptr=tmp
         return support
     def applySupport(self,PyComplex2D c):
         self.thisptr.apply_support(deref(c.thisptr))
     def projectIntensity(self, PyComplex2D c):
         self.thisptr.project_intensity(deref(c.thisptr))
     def scaleIntensity(self,PyComplex2D c):
         self.thisptr.scale_intensity(deref(c.thisptr))
     def propagateToDetector(self,PyComplex2D c):
         self.thisptr.propagate_to_detector(deref(c.thisptr))
     def propagateFromDetector(self,PyComplex2D c):
         self.thisptr.propagate_from_detector(deref(c.thisptr))
     def setFFTWType(self,int type):    
         self.thisptr.set_fftw_type(type)
     def setComplexConstraint(self,PyTransmissionConstraint ptc):
         self.thisptr.set_complex_constraint(deref(ptc.thisptr))
    
         