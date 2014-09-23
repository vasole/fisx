#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from SimpleSpecfile cimport *
    
cdef class PySimpleSpecfile:
    cdef SimpleSpecfile *thisptr

    def __cinit__(self, name):
        name = toBytes(name)
        self.thisptr = new SimpleSpecfile(name)

    def __dealloc__(self):
        del self.thisptr

    def getNumberOfScans(self):
        return self.thisptr.getNumberOfScans()

    #def getScanHeader(self, int scanIndex):
    #    return self.thisptr.getScanHeader(scanIndex)

    def getScanLabels(self, int scanIndex):
        return self.thisptr.getScanLabels(scanIndex)

    def getScanData(self, int scanIndex):
        return self.thisptr.getScanData(scanIndex)
