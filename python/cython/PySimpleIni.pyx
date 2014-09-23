#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from SimpleIni cimport *
    
cdef class PySimpleIni:
    cdef SimpleIni *thisptr

    def __cinit__(self, name):
        name = toBytes(name)
        self.thisptr = new SimpleIni(name)

    def __dealloc__(self):
        del self.thisptr

    def getKeys(self):
        return self.thisptr.getSections()

    def readKey(self, key):
        key = toBytes(key)
        return self.thisptr.readSection(key)
