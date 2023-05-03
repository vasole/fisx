#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2023 European Synchrotron Radiation Facility
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

from Beam cimport *
from Detector cimport *
from fisx.FisxCythonTools import toBytes, toBytesKeys, toBytesKeysAndValues, toString,  toStringKeys, toStringKeysAndValues, toStringList

cdef class PyBeam:
    cdef Beam *thisptr

    def __cinit__(self):
        self.thisptr = new Beam()

    def __dealloc__(self):
        del self.thisptr

    def setBeam(self, energies, weights=None):
        cdef std_vector[int] characteristic = std_vector[int]()
        cdef std_vector[double] divergency = std_vector[double]()
        if not hasattr(energies, "__len__"):
            energies = numpy.array([energies], dtype=numpy.float64)
        else:
            energies = numpy.array(energies, dtype=numpy.float64, copy=False)
        if weights:
            if not hasattr(weights, "__len__"):
                weights = numpy.array([weights], numpy.float64)
            else:
                weights = numpy.array(weights, dtype=numpy.float64, copy=False)
        else:
            weights = numpy.ones(energies.shape, dtype=numpy.float64)
        return self.thisptr.setBeam(energies,
                                    weights,
                                    characteristic,
                                    divergency)

    def getBeamAsDoubleVectors(self):
        output = self.thisptr.getBeamAsDoubleVectors()
        return output[0], output[1]
