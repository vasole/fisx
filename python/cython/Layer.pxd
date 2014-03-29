#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *

cdef extern from "fisx_layer.h":
    cdef cppclass Layer:
        Layer(std_string, double, double, double)

        std_vector[double] getTransmission(std_vector[double], Elements, double)
