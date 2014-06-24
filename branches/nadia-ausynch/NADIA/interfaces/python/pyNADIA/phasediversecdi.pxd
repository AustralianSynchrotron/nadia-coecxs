cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D
from basecdi cimport BaseCDI
from libcpp.string cimport string
from libcpp cimport bool
from libcpp cimport set
cimport cython

cdef extern from "PhaseDiverseCDI.h":
    cdef cppclass PhaseDiverseCDI:
        PhaseDiverseCDI(double beta, double gamma, bool parallel,int granularity)
        void add_new_position(BaseCDI * local, double x,double y,double alpha)
        void remove_positions()
        void initialise_estimate()
        int iterate()
        void set_iterations_per_cycle(int iterations)
        void set_feedback_parameter(double beta)
        void set_amplification_factor(double gamma)
        void set_probe_scaling(int n_probe, double alpha)
        Complex_2D * get_transmission()
        void set_transmission(Complex_2D & new_transmission)
        void adjust_positions(int type, bool forwards, int x_min, int x_max, int y_min, int y_max, double step_size)
        double get_final_x_position(int n_probe)
        double get_final_y_position(int n_probe)
    enum:
        CROSS_CORRELATION "PhaseDiverseCDI::CROSS_CORRELATION"
        MINIMUM_ERROR "PhaseDiverseCDI::MINIMUM_ERROR"

cdef class PyPhaseDiverseCDI:
    cdef PhaseDiverseCDI *thisptr
    cdef set cdiobjs