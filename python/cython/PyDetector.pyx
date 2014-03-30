import numpy
#cimport numpy as np
cimport cython

from cython.operator cimport dereference as deref
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *
from Detector cimport *

cdef class PyDetector:
    cdef Detector *thisptr

    def __cinit__(self, std_string materialName, double density=1.0, double thickness=1.0, double funny=1.0):
        self.thisptr = new Detector(materialName, density, thickness, funny)

    def __dealloc__(self):
        del self.thisptr

    def getTransmission(self, energies, PyElements elementsLib, double angle=90.):
        if not hasattr(energies, "__len__"):
            energies = numpy.array([energies], numpy.float)
        return self.thisptr.getTransmission(energies, deref(elementsLib.thisptr), angle)

    def setActiveArea(self, double area):
        self.thisptr.setActiveArea(area)

    def setDiameter(self, double value):
        self.thisptr.setDiameter(value)

    def getActiveArea(self):
        return self.thisptr.getActiveArea()

    def getDiameter(self):
        return self.thisptr.getDiameter()

    def setDistance(self, double value):
        self.thisptr.setDistance(value)

    def getDistance(self):
        return self.thisptr.getDistance()
