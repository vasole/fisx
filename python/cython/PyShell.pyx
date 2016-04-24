#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2016 European Synchrotron Radiation Facility
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

#from libcpp.string cimport string as std_string
#from libcpp.vector cimport vector as std_vector
#from libcpp.map cimport map as std_map

from Shell cimport *

cdef class PyShell:
    cdef Shell *thisptr

    def __cinit__(self, name):
        name = toBytes(name)
        self.thisptr = new Shell(name)

    def __dealloc__(self):
        del self.thisptr

    def setRadiativeTransitions(self, transitions, std_vector[double] values):
        if sys.version > "3.0":
            transitions = [toBytes(x) for x in transitions]
        self.thisptr.setRadiativeTransitions(transitions, values)

    def setNonradiativeTransitions(self, transitions, std_vector[double] values):
        if sys.version > "3.0":
            transitions = [toBytes(x) for x in transitions]
        self.thisptr.setNonradiativeTransitions(transitions, values)

    def getAugerRatios(self):
        return self.thisptr.getAugerRatios()

    def getCosterKronigRatios(self):
        return self.thisptr.getCosterKronigRatios()

    def getFluorescenceRatios(self):
        return self.thisptr.getFluorescenceRatios()

    def getRadiativeTransitions(self):
        return self.thisptr.getRadiativeTransitions()

    def getNonradiativeTransitions(self):
        return self.thisptr.getNonradiativeTransitions()

    def getDirectVacancyTransferRatios(self, subshell):
        return self.thisptr.getDirectVacancyTransferRatios(toBytes(subshell))

    def setShellConstants(self, shellConstants):
        if sys.version > "3.0":
            shellConstants = toBytesKeys(shellConstants)
        self.thisptr.setShellConstants(shellConstants)

    def getShellConstants(self):
        return self.thisptr.getShellConstants()
