#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.pair cimport pair as std_pair

cdef extern from "fisx_simplespecfile.h":
    cdef cppclass SimpleSpecfile:
        SimpleSpecfile(std_string) except +
        void setFileName(std_string) except +
        int getNumberOfScans() except +
        std_vector[std_string] getScanHeader(int scanIndex) except +
        std_vector[std_string] getScanLabels(int scanIndex) except +
        std_vector[std_vector[double] ] getScanData(int scanIndex) except +
