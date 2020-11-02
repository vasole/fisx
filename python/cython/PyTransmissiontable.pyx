#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2020 European Synchrotron Radiation Facility
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

from TransmissionTable cimport *

cdef class PyTransmissionTable:
    cdef TransmissionTable *thisptr

    def __cinit__(self):
        self.thisptr = new TransmissionTable()

    def __dealloc__(self):
        del self.thisptr

    def getName(self):
        return toString(self.thisptr.getName())

    def setName(self, name):
        name = toBytes(name)
        self.thisptr.setName(name)

    def getComment(self):
        return toString(self.thisptr.getName())

    def setComment(self, comment):
        comment = toBytes(comment)
        self.thisptr.setComment(comment)

    def setTransmissionTableFromLists(self,
                             std_vector[double] energy, \
                             std_vector[double] transmission, \
                             name="",
                             comment=""):
        self.thisptr.setTransmissionTable(energy,
                                          transmission, \
                                          toBytes(name),
                                          toBytes(comment))

    def setTransmissionTableFromLists(self,
                             std_map[double, double] transmissionTable, \
                             name="",
                             comment=""):
        self.thisptr.setTransmissionTable(transmissionTable,
                                          toBytes(name),
                                          toBytes(comment))

    def getTransmissionTable(self):
        return self.thisptr.getTransmissionTable
