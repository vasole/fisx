#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.vector cimport vector as std_vector

cdef extern from "fisx_math.h" namespace "fisx":
    cdef cppclass Math:
        Math()

        double E1(double) except +
        
        double En(int, double) except +

        double deBoerD(double) except +

        double deBoerL0(double, double, double, double, double) except +

        double erf(double)

        double erfc(double)

        double hypermet(double, double, double, double, double, 
                        double, double, double, double) except + 
