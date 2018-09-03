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

    def getName(self):
        return toString(self.thisptr.getName())

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
        return toStringKeys(self.thisptr.getComposition())
