import sys
cimport cython

from cython.operator cimport dereference as deref
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Material cimport *

cdef class PyMaterial:
    cdef Material *thisptr

    def __cinit__(self, materialName, double density=1.0, double thickness=1.0, comment=""):
        materialName = toBytes(materialName)
        comment = toBytes(comment)
        self.thisptr = new Material(materialName, density, thickness, comment)

    def __dealloc__(self):
        del self.thisptr

    def setName(self, name):
        name = toBytes(name)
        self.thisptr.setName(name)

    def setCompositionFromLists(self, elementList, std_vector[double] massFractions):
        if sys.version > "3.0":
            elementList = [toBytes(x) for x in elementList]
        self.thisptr.setComposition(elementList, massFractions)

    def setComposition(self, composition):
        if sys.version > "3.0":
            composition = toBytesKeys(composition)
        self.thisptr.setComposition(composition)

    def getComposition(self):
        return self.thisptr.getComposition()
