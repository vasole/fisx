#cimport numpy as np
cimport cython

from cython.operator cimport dereference as deref
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Material cimport *

cdef class PyMaterial:
    cdef Material *thisptr

    def __cinit__(self, std_string materialName, double density=1.0, double thickness=1.0, std_string comment=""):
        self.thisptr = new Material(materialName, density, thickness, comment)

    def __dealloc__(self):
        del self.thisptr

    def setName(self, std_string name):
        self.thisptr.setName(name)

    def setCompositionFromLists(self, std_vector[std_string] elementList, std_vector[double] massFractions):
        self.thisptr.setComposition(elementList, massFractions)

    def setComposition(self, std_map[std_string, double] composition):
        self.thisptr.setComposition(composition)

    def getComposition(self):
        return self.thisptr.getComposition()
