#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2018 European Synchrotron Radiation Facility
#
# This file is part of the fisx X-ray developed by V.A. Sole
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
import numpy
import sys
cimport cython

from cython.operator cimport dereference as deref
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *
from Detector cimport *
from FisxCythonTools import toBytes, toBytesKeys, toBytesKeysAndValues, toString,  toStringKeys, toStringKeysAndValues, toStringList

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

    def getEscape(self, double energy, PyElements elementsLib, label="", int update=1):
        label_ = toBytes(label)
        if sys.version < "3.0":
            if update:
                return self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label_, 1)
            else:
                return self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label_, 0)
        else:
            if update:
                return toStringKeysAndValues(self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label_, 1))
            else:
                return toStringKeysAndValues(self.thisptr.getEscape(energy, deref(elementsLib.thisptr), label_, 0))

    def getEscapePeakEnergyThreshold(self):
        return self.thisptr.getEscapePeakEnergyThreshold()

    def getEscapePeakIntensityThreshold(self):
        return self.thisptr.getEscapePeakIntensityThreshold()

    def getEscapePeakNThreshold(self):
        return self.thisptr.getEscapePeakNThreshold()

    def getEscapePeakAlphaIn(self):
        return self.thisptr.getEscapePeakAlphaIn()

    def getThickness(self):
        return self.thisptr.getThickness()

    def getDensity(self):
        return self.thisptr.getDensity()

    def getComposition(self, PyElements elementsLib):
        if sys.version < "3.0":
            return self.thisptr.getComposition(deref(elementsLib.thisptr))
        else:
            return toStringKeysAndValues(self.thisptr.getComposition(deref(elementsLib.thisptr)))
