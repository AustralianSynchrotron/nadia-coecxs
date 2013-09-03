cimport cython
cdef extern from "Double_2D.h":
    cdef cppclass Real_2D[T]:
        Real_2D() except+
        Real_2D(int x_size, int y_size) except+
        Real_2D(const Real_2D & object) except+
        #Real_2D& operator=(const Real_2D& rhs)
        void allocate_memory(int x_size, int y_size)
        void copy(const Real_2D[T] & other)
        void set(int x, int y, T value)
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
        T * array
        void copy(const Real_2D[T] & other_array)
    ctypedef Real_2D[float] Double_2D

cdef class PyDouble2D:
     cdef Double_2D *thisptr