"""! @package pyNADIA.basecdi
Empty placeholder Python wrapper BaseCDI class. 
"""

cdef class PyBaseCDI:
    """!
    This class is only useful for use by PyPhaseDiverseCDI methods as a base class for PyPlanarCDI
    and PyFresnelCDI.
    """
    def __dealloc__(self):
        del self.thisptr

