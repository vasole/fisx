import numpy
import sys
cimport cython

from cython.operator cimport dereference as deref
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *
from Detector cimport *
from FisxCythonTools import toBytes, toBytesKeys, toBytesKeysAndValues, toStringKeys, toStringKeysAndValues

cdef class PyDetector:
    cdef Detector *thisptr

    def __cinit__(self, materialName, double density=1.0, double thickness=1.0, double funny=1.0):
        self.thisptr = new Detector(toBytes(materialName), density, thickness, funny)

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

    def setMaximumNumberOfEscapePeaks(self, int n):
        self.thisptr.setMaximumNumberOfEscapePeaks(n)

    def getEscape(self, double energy, PyElements elementsLib, std_string label="", int update=1):
        if sys.version < "3.0":
            if update:
                return self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label, 1)
            else:
                return self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label, 0)
        else:
            if update:
                return toStringKeysAndValues(self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label, 1))
            else:
                return toStringKeysAndValues(self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label, 0))
