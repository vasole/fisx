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
import numpy
#cimport numpy as np
cimport cython

from libcpp.vector cimport vector as std_vector

from Math cimport *

cdef class PyMath:
    cdef Math *thisptr

    def __cinit__(self):
        self.thisptr = new Math()

    def __dealloc__(self):
        del self.thisptr

    def E1(self, double x):
        return self.thisptr.E1(x)

    def En(self, int n, double x):
        return self.thisptr.En(n, x)

    def deBoerD(self, double x):
        return self.thisptr.deBoerD(x)

    def deBoerL0(self, double mu1, double mu2, double muj, double density = 0.0, double thickness = 0.0):
        """
        The case the product density * thickness is 0.0 is for calculating the thick target limit
        """
        return self.thisptr.deBoerL0(mu1, mu2, muj, density, thickness)

    def erf(self, double x):
        """
        Calculate the error function erf(x)
        """
        return self.thisptr.erf(x)

    def erfc(self, double x):
        """
        Calculate the complementary error function erfc(x)
        """
        return self.thisptr.erfc(x)
