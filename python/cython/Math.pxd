#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.vector cimport vector as std_vector

cdef extern from "fisx_math.h":
    cdef cppclass Math:
        Math()

        double E1(double) except +

        double deBoerD(double) except +
