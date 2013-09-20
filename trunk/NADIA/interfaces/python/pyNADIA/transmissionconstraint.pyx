cimport transmissionconstraint
from double2d cimport PyDouble2D, Double_2D
from complex2d cimport PyComplex2D,Complex_2D
from transmissionconstraint cimport TransmissionConstraint, ComplexConstraint,PyComplexConstraint
from libcpp cimport bool
from cython.operator cimport dereference as deref

cdef class PyComplexConstraint:
    
    def __cinit__(self,PyDouble2D region, double alpha1,double alpha2):
        self.thisptr=new ComplexConstraint(deref(region.thisptr),alpha1,alpha2)
        
    def __dealloc__(self):
        del self.thisptr
        
    def setFixedC(self, c_mean):
        self.thisptr.set_fixed_c(c_mean)
    def setCMean(self, c_mean):
        self.thisptr.set_c_mean(c_mean)
    def getCMean(self):
        return self.thisptr.get_c_mean()
    def setAlpha1(self, double alpha1):
        self.thisptr.set_alpha1(alpha1)
    def setAlpha2(self,double alpha2):
        self.thisptr.set_alpha2(alpha2)
    def getNewMag(self,old_mag,old_phase):
        self.thiptr.get_new_mag(old_mag,old_phase)
    def getNewPhase(self,old_mag,old_phase):
        self.thisptr.get_new_phase(old_mag,old_phase)
    def getRegion(self):
        region = PyDouble2D()
        region.thisptr=self.thisptr.get_region()
        return region
    
cdef class PyTransmissionConstraint:

    def __cinit__(self):
        self.thisptr=new TransmissionConstraint()
    def __dealloc__(self):
        del self.thisptr
    def deleteComplexConstraintRegions(self):
        self.thisptr.delete_complex_constraint_regions()
    def addComplexConstraint(self, PyComplexConstraint nc):
        self.thisptr.add_complex_constraint(deref(nc.thisptr))
    def setChargeFlipping(self, enable, flip_sign=1):
        self.thisptr.set_charge_flipping(enable,flip_sign)
    def setEnforceUnity(self,enable):
        self.thisptr.set_enforce_unity(enable)
    def applyConstraint(self, PyComplex2D transmission):
        self.thisptr.apply_constraint(deref(transmission.thisptr))

       