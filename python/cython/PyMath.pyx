import numpy
#cimport numpy as np
cimport cython

from libcpp.vector cimport vector as std_vector

from Math cimport *

cdef class PyMath:
    cdef Math *thisptr

    def __cinit__(self):
        self.thisptr = new Math()

    def __dealloc__(self):
        del self.thisptr

    def E1(self, double x):
        return self.thisptr.E1(x)

    def deBoerD(self, double x):
        return self.thisptr.deBoerD(x)

    def deBoerL0(self, double mu1, double mu2, double muj, double density = 0.0, double thickness = 0.0):
        """
        The case the product density * thickness is 0.0 is for calculating the thick target limit
        """
        return self.thisptr.deBoerL0(mu1, mu2, muj, density, thickness)
