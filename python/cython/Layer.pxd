#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map
from libcpp.map cimport pair as std_pair

from Elements cimport *
from Material cimport *

cdef extern from "fisx_layer.h":
    cdef cppclass Layer:
        Layer(std_string, double, double, double)

        std_vector[double] getTransmission(std_vector[double], Elements, double)

        void setMaterial(Material)

        std_vector[std_pair[std_string, double]] getPeakFamilies(double, Elements) except +
