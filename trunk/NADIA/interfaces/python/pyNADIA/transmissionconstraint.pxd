from double2d cimport Double_2D
from complex2d cimport Complex_2D
from libcpp cimport bool

cdef extern from "TransmissionConstraint.h":
    cdef cppclass ComplexConstraint:
        ComplexConstraint(Double_2D & region, double alpha1, double alpha2)
        void set_fixed_c(double c_mean)
        void set_c_mean(double c_mean)
        double get_c_mean()
        void set_alpha1(double alpha1)
        void set_alpha2(double alpha2)
        double get_new_mag(double old_mag, double old_phase)
        double get_new_phase(double old_mag, double old_phase)
        Double_2D * get_region()     
    cdef cppclass TransmissionConstraint:
        TransmissionConstraint()
        void delete_complex_constraint_regions()
        void add_complex_constraint(ComplexConstraint & new_constraint)
        void set_charge_flipping(bool enable, int flip_sign=1)
        void set_enform_unity(bool enable)
        void set_custom_constraint(void (*custom_constraint)(Complex_2D & transmission))
        void apply_constraint(Complex_2D & transmission)

cdef class PyComplexConstraint:
     cdef ComplexConstraint *thisptr

cdef class PyTransmissionConstraint:
     cdef TransmissionConstraint *thisptr