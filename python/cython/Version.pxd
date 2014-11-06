#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string

cdef extern from "fisx_version.h" namespace "fisx":

    std_string fisxVersion()
