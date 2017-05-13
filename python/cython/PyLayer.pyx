#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2017 European Synchrotron Radiation Facility
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
#cimport numpy as np
cimport cython

from cython.operator cimport dereference as deref
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map
from libcpp.map cimport pair as std_pair
from operator import itemgetter

from Elements cimport *
from Material cimport *
from Layer cimport *

cdef class PyLayer:
    cdef Layer *thisptr

    def __cinit__(self, materialName, double density=1.0, double thickness=1.0, double funny=1.0):
        self.thisptr = new Layer(toBytes(materialName), density, thickness, funny)

    def __dealloc__(self):
        del self.thisptr

    def getComposition(self, PyElements elementsLib):
        """
        getComposition(elementsLib)

        Given a reference to an elements library, it gives back a dictionary where the keys are the
        elements and the values the mass fractions.
        """
        return self.thisptr.getComposition(deref(elementsLib.thisptr))

    def getTransmission(self, energies, PyElements elementsLib, double angle=90.):
        """
        getTransmission(energies, ElementsLibraryInstance, angle=90.)

        Given a list of energies and a reference to an elements library returns
        the layer transmission according to the incident angle (default 90.)
        """
        if not hasattr(energies, "__len__"):
            energies = numpy.array([energies], numpy.float)
        return self.thisptr.getTransmission(energies, deref(elementsLib.thisptr), angle)

    def setMaterial(self, PyMaterial material):
        """
        setMaterial(MaterialInstance)

        Set the material of the layer. It has to be an instance!
        """
        self.thisptr.setMaterial(deref(material.thisptr))

    def getPeakFamilies(self, double energy, PyElements elementsLib):
        """
        getPeakFamilies(energy, ElementsLibraryInstance)

        Given an energy and a reference to an elements library return dictionarys.
        The key is the peak family ("Si K", "Pb L1", ...) and the value the binding energy.
        """
        tmpResult = self.thisptr.getPeakFamilies(energy, deref(elementsLib.thisptr))
        return sorted(tmpResult, key=itemgetter(1))

