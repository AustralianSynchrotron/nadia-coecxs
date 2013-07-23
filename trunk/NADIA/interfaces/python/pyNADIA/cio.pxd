from libcpp.string cimport string
from libcpp cimport bool
cimport double2d
cimport complex2d
from double2d cimport Double_2D
from complex2d cimport Complex_2D

cdef extern from "io.h":
     cdef int read(string filename, Double_2D & data)
     cdef int read_ppm(string file_name, Double_2D & data)
     cdef int read_spec(string file_name, Double_2D & data)
     cdef int read_tiff(string file_name, Double_2D & data)
     cdef int read_hdf4(string file_name, Double_2D & data)
     cdef int read_dbin(string file_name, int nx, int ny, Double_2D & data)
     cdef int read_cplx(string file_name, Complex_2D & complex)
     cdef int write_ppm(string file_name, const Double_2D & data)
     cdef int write_ppm(string file_name, const Double_2D & data, bool log_scale)
     cdef int write_ppm(string file_name, const Double_2D & data, bool log_scale, double min, double max)
     cdef int write_spec(string file_name, const Double_2D & data)
     cdef int write_spec(string file_name, const Double_2D & data, int nx, int ny)
     cdef int write_dbin(string file_name, const Double_2D & data)
     cdef int write_cplx(string file_name, const Complex_2D & complex)
     cdef int write_tiff(string file_name, const Double_2D & data)
     cdef int write_tiff(string file_name, const Double_2D & data, bool log_scale)
     cdef int write_tiff(string file_name, const Double_2D & data, bool log_scale, double min, double max)
    


