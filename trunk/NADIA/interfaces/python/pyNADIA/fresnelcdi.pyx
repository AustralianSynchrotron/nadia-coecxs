from fresnelcdi cimport FresnelCDI, PyFresnelCDI
from double2d cimport Double_2D, PyDouble2D
from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp cimport bool
from complex2d cimport Complex_2D, PyComplex2D


cdef class PyFresnelCDI: 
     #cdef nx,ny
     #cdef FresnelCDI *thisptr    
   # def __cinit__(self):
   #     if self.thisptr:
  #         del self.thisptr
  #      self.thisptr = new FresnelCDI()        
     def __cinit__(self,PyComplex2D object_estimate,PyComplex2D wf,double wavelength,double ftd,double fts,double pixel_size,double norm):
        if self.thisptr:
            del self.thisptr
        cdef FresnelCDI *newptr 
        self.thisptr = newptr
        self.thisptr = new FresnelCDI(deref(object_estimate.thisptr), deref(wf.thisptr),wavelength,ftd,fts,pixel_size,norm)
        self.nx=object_estimate.thisptr.get_size_x()
        self.ny=object_estimate.thisptr.get_size_y()   
     def __dealloc__(self):
        if self.thisptr:
            del self.thisptr     
     def iterate(self):
         self.thisptr.iterate()
     def setSupport(self,PyDouble2D support):
        cdef Double_2D csupport = deref(support.thisptr)
        self.thisptr.set_support(csupport,False)
     def getTransmissionFunction(self):
        trans=PyComplex2D(self.thisptr.get_illumination_at_sample().get_size_x(),self.thisptr.get_illumination_at_sample().get_size_y())
        self.thisptr.get_transmission_function(deref(trans.thisptr))
        return trans
     def applyShrinkwrap(self,int x, int y):
        self.thisptr.apply_shrinkwrap(x,y)
     def getIlluminationAtSample(self):     
        result = PyComplex2D(self.nx,self.ny)
        result.thisptr = self.thisptr.get_illumination_at_sample().clone()
        return result 
     def getError(self):
         return self.thisptr.get_error()
     def initialiseEstimate(self,int n=0):
         self.thisptr.initialise_estimate(n)
     def setIntensity(self, PyDouble2D intensity):
         self.thisptr.set_intensity(deref(intensity.thisptr))
     def setAlgorithm(self,string algorithm):
         self.thisptr.set_algorithm(self.thisptr.getAlgFromName(algorithm))
     def getBestResult(self,double & error, int index=0):
         cdef Complex_2D * cresult
         cresult= self.thisptr.get_best_result(error,index)
         result=PyComplex2D(cresult.get_size_x(),cresult.get_size_y())
     def setSupport(self,PyDouble2D support, bool soften=False):
         self.thisptr.set_support(deref(support.thisptr),soften)
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
     
                  
         