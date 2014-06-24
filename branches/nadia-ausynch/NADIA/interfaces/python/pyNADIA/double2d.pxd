cimport cython

cdef extern from "Double_2D.h":
    cdef cppclass Real_2D[T]:
        Real_2D() except+
        Real_2D(int x_size, int y_size) except+
        Real_2D(const Real_2D & object) except+
        void allocate_memory(int x_size, int y_size)
        void copy(const Real_2D[T] & other)
        T get_max()
        int get_size_x()
        int get_size_y()
        T get(int x, int y)
        T get_sum()
        T get_abs_sum()
        T get_min()
        void add(const Real_2D[T] & other)
        void add(const Real_2D[T] & other, double norm)
        void scale(double scale_factor)
        void square()
        void sq_root()
        void set(int x,int y, T value)
        T * c_array()
IF DOUBLE_PRECISION !='1':
        ctypedef Real_2D[float] Double_2D
ELSE:
        ctypedef Real_2D[double] Double_2D

cdef class PyDouble2D:
    cdef Double_2D *thisptr
    cdef object __weakref__