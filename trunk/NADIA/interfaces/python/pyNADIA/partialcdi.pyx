cimport polycdi
cimport double2d
cimport complex2d
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from partialcdi cimport PartialCDI, PyPartialCDI 
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint

from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef class PyPartialCDI:
     def __cinit__(self,PyComplex2D object_estimate, double lcx, double lcy, \
                                    double pxsize, double pysize, double energy, double zsd, int nleg, int nmode, unsigned int n_best):
         self.thisptr=new PartialCDI(deref(object_estimate.thisptr),lcx,lcy,pxsize,pysize,energy,zsd,nleg,nmode,n_best)
         self.nx=object_estimate.thisptr.get_size_x()
         self.ny=object_estimate.thisptr.get_size_y()
     def __dealloc__(self):
         del self.thisptr
     def initialiseEstimate(self, int seed=0):
         self.thisptr.initialise_estimate(seed)
     def iterate(self):
         return self.thisptr.iterate()
     def getTransmission(self):
         trans = PyComplex2D(self.nx,self.ny)
         cdef Complex_2D * tmp = new Complex_2D(self.thisptr.get_transmission())
         trans.thisptr=tmp
         return trans
     def setTransmission(self, PyComplex2D newTrans):
         self.thisptr.set_transmission(deref(newTrans.thisptr))
     def setThreshold(self, double new_threshold):
         self.thisptr.set_threshold(new_threshold)
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