import phasediversecdi
cimport double2d
cimport complex2d
from basecdi cimport BaseCDI,PyBaseCDI
from double2d cimport Double_2D, PyDouble2D
from complex2d cimport Complex_2D, PyComplex2D
from phasediversecdi cimport PhaseDiverseCDI, PyPhaseDiverseCDI
from transmissionconstraint cimport TransmissionConstraint, PyTransmissionConstraint

from libcpp.string cimport string
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef class PyPhaseDiverseCDI:
     def __cinit__(self,double beta=1.0,double gamma=1.0, bool parallel=False, int granularity=1):
         self.thisptr=new PhaseDiverseCDI(beta,gamma,parallel,granularity)
     def __dealloc__(self):
         del self.thisptr
     def addNewPosition(self, PyBaseCDI local, double x=0.0, double y=1.0, double alpha=1.0):
         self.thisptr.add_new_position(local.thisptr,x,y,alpha)
     def initialiseEstimate(self):
         self.thisptr.initialise_estimate()
     def iterate(self): 
         self.thisptr.iterate()
     def setIterationsPerCycle(self,int iterations):
         self.thisptr.set_iterations_per_cycle(iterations)
     def setFeedbackParameter(self, double beta):
         self.thisptr.set_feedback_parameter(beta)
     def setAmplificationFactor(self,double gamma):
         self.thisptr.set_amplification_factor(gamma)
     def setProbeScaling(self, int n_probe,double alpha):
         self.thisptr.set_probe_scaling(n_probe,alpha)
     def setTransmission(self, PyComplex2D trans):
         self.thisptr.set_transmission(deref(trans.thisptr))
     def adjustPositions(self, int type, bool forwards=True, int x_min=-50, int x_max=50, int y_min=-50, int y_max=50, double step_size=4.0):
         self.thisptr.adjust_positions(type,forwards,x_min,x_max,y_min,y_max,step_size)
     def getFinalXPosition(self,int n_probe):
         return self.thisptr.get_final_x_position(n_probe)    
     def getFinalYPosition(self, int n_probe):
         return self.thispr.get_final_y_position(n_probe)
     