cimport polycdi
cimport double2d
cimport complex2d
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from polycdi cimport PolyCDI, PyPolyCDI 
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint

from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef class PyPolyCDI:
     def __cinit__(self,PyComplex2D object_estimate):
         self.thisptr=new PolyCDI(deref(object_estimate.thisptr))
         self.nx=object_estimate.thisptr.get_size_x()
         self.ny=object_estimate.thisptr.get_size_y()
     def __dealloc__(self):
         del self.thisptr
     def initialiseEstimate(self, int seed=0):
         self.thisptr.initialise_estimate(seed)
     def scaleIntensity(self,PyComplex2D c):
         self.thisptr.scale_intensity(deref(c.thisptr))
     def expandWL(self,PyComplex2D c):
         self.thisptr.expand_wl(deref(c.thisptr))
     def sumIntensity(self, pycomp2dlist):
         cdef vector[Complex_2D] * comp2dlist = new vector[Complex_2D]()
         result = PyDouble2D()
         for c in pycomp2dlist:
             comp2dlist.push_back(deref((<PyComplex2D>c).thisptr))
         cdef Double_2D tmp=self.thisptr.sum_intensity(deref(comp2dlist))
         result.thisptr=&tmp
         return result
     def iterate(self):
         return self.thisptr.iterate()
     def updateTransmission(self):
         self.thisptr.update_transmission()
     def setIterationsPerCycle(self, int iterations):
         self.thisptr.set_iterations_per_cycle(iterations)
     def getTransmission(self):
         trans = PyComplex2D(self.nx,self.ny)
         cdef Complex_2D * tmp = new Complex_2D(self.thisptr.get_transmission())
         trans.thisptr=tmp
         return trans
     def setTransmission(self, PyComplex2D newTrans):
         self.thisptr.set_transmission(deref(newTrans.thisptr))
     def setSpectrum(self, PyDouble2D spec):
         self.thisptr.set_spectrum(deref(spec.thisptr))
     def setSpectrum(self, string file_name):
         self.thisptr.set_spectrum(file_name)
     def getIntensity(self):
         intens=PyDouble2D()
         cdef Double_2D * tmp = new Double_2D(self.thisptr.get_intensity())
         intens.thisptr=tmp
         return intens
     def propagateModesToDetector(self):
         modes=PyDouble2D()
         cdef Double_2D * tmp = new Double_2D(self.thisptr.propagate_modes_to_detector())
         modes.thisptr=tmp
         return modes
     def getMode(self, int modenum):
         mode =PyComplex2D(self.nx,self.ny)
         cdef Complex_2D * tmp = new Complex_2D(self.thisptr.get_mode(modenum))
         mode.thisptr=tmp
         return mode
     def propagateToDetector(self,PyComplex2D c):
         self.thisptr.propagate_to_detector(deref(c.thisptr))
     def propagateFromDetector(self,PyComplex2D c):
         self.thisptr.propagate_from_detector(deref(c.thisptr))
             
     def getBestResult(self,double & error, int index=0):
         cdef Complex_2D * cresult
         cresult= self.thisptr.get_best_result(error,index)
         result=PyComplex2D(cresult.get_size_x(),cresult.get_size_y())
     def setSupport(self,PyDouble2D support, bool soften=False):
         self.thisptr.set_support(deref(support.thisptr),soften)
     def setIntensity(self,PyDouble2D intensity):
         self.thisptr.set_intensity(deref(intensity.thisptr))
     def setBeamStop(self,PyDouble2D beamstop):
         self.thisptr.set_beam_stop(deref(beamstop.thisptr))
     def getError(self):
         return self.thisptr.get_error()
     def getSupport(self):
         result =PyDouble2D()
     def setAlgorithm(self,int alg):
         self.thisptr.set_algorithm(alg)
     def resetBest(self):
         self.thisptr.reset_best()
     def set_fftw_type(self,int type):
         self.thisptr.set_fftw_type(type)  

         