from double2d cimport Double_2D,PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from libcpp.string cimport string
from libcpp cimport bool
from transmissionconstraint cimport TransmissionConstraint,ComplexConstraint
from cython.operator cimport dereference as deref

cdef class PyBaseCDI:
     pass
    # def __cinit__(self, PyComplex2D complex, int n_best=0):
    #     self.thisptr=new BaseCDI(deref(complex.thisptr),n_best)
    # def iterate(self):
   #      return self.thisptr.iterate()
   #  def getBestResult(self,double & error, int index=0):
   #      cdef Complex_2D * cresult
   #      cresult= self.thisptr.get_best_result(error,index)
   #      result=PyComplex2D(cresult.get_size_x(),cresult.get_size_y())
   #  def setSupport(self,PyDouble2D support, bool soften=False):
   #      self.thisptr.set_support(deref(support.thisptr),soften)
  #   def setIntensity(self,PyDouble2D intensity):
  #       self.thisptr.set_intensity(deref(intensity.thisptr))
   #  def setBeamStop(self,PyDouble2D beamstop):
   #      self.thisptr.set_beam_stop(deref(beamstop.thisptr))
  #   def getError(self):
  #       return self.thisptr.get_error()
   #  def getSupport(self):
   #      result =PyDouble2D()
   #  def setAlgorithm(self,int alg):
  #       self.thisptr.set_algorithm(alg)
  #   def resetBest(self):
  #       self.thisptr.reset_best()
  #   def set_fftw_type(self,int type):
  #       self.thisptr.set_fftw_type(type)
     
