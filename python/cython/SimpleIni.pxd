#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

cdef extern from "fisx_simpleini.h" namespace "fisx":
    cdef cppclass SimpleIni:
        SimpleIni() except +
        SimpleIni(std_string) except +
        void setFileName(std_string) except +
        std_vector[std_string] getSections() except +
        std_map[std_string, std_string] readSection(std_string) except +