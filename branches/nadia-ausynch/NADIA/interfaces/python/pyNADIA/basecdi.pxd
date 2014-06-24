from double2d cimport Double_2D
from complex2d cimport Complex_2D
from libcpp.string cimport string
from libcpp cimport bool
from transmissionconstraint cimport TransmissionConstraint

cdef extern from "BaseCDI.h":
    cdef enum:
        ER, BIO, BOO, HIO, DM, SF, ASR, HPR, RAAR, CUSTOM
    cdef enum:
        PSF, PFS, PS, PF, PI 
    cdef cppclass BaseCDI:
        # BaseCDI(Complex_2D & complex) except+
        # BaseCDI(Complex_2D & complex, unsigned int n_best) except+
        int getAlgFromName(string algorithm_name)
        # int iterate()
        Complex_2D * get_best_result(double & error)
        Complex_2D * get_best_result(double & error, int index)
        void set_intensity(const Double_2D & detector_intensity)
        void set_beam_stop(const Double_2D & beam_stop_region)
        void set_relaxation_parameter(double relaxation_parameter)
        int get_size_x()
        int get_size_y()
        const Complex_2D & get_exit_surface_wave()
        void set_exit_surface_wave(Complex_2D & esw)
        void set_algorithm(int alg)
        void set_custom_algorithm(double m1, double m2, double m3, \
                double m4, double m5, double m6, \
                double m7, double m8, \
                double m9, double m10)
        void print_algorithm()
        double get_error()
        void apply_shrinkwrap()
        void apply_shrinkwrap(double gauss_width)
        void apply_shrinkwrap(double gauss_width, double threshold)
        const Double_2D & get_support()
        # void apply_support(Complex_2D & c)
        void project_intensity(Complex_2D & c)
        # void scale_intensity(Complex_2D & c)
        void set_fftw_type(int type)
        void set_complex_constraint(TransmissionConstraint & trans_constraint)
        void reset_best()

cdef class PyBaseCDI:
    cdef BaseCDI * thisptr
