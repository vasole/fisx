#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

cdef extern from "fisx_material.h" namespace "fisx":
    cdef cppclass Material:
        Material(std_string, double, double, std_string) except +

        void setName(std_string) except +

        void setComposition(std_map[std_string, double]) except +

        void setComposition(std_vector[std_string], std_vector[double])  except +

        std_map[std_string, double] getComposition() except +