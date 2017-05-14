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

    def deBoerX(self, double p, double q, double d1, double d2, double mu_1_j, double mu_2_j, double mu_b_d_t = 0.0):
        """
        static double deBoerX(const double & p, const double & q, \
                              const double & d1, const double & d2, \
                              const double & mu_1_j, const double & mu_2_j, \
                              const double & mu_b_j_d_t = 0.0);
        For multilayers
        p and q following article
        d1 is the product density * thickness of fluorescing layer
        d2 is the product density * thickness of layer j originating the secondary excitation
        mu_1_j is the mass attenuation coefficient of fluorescing layer at j excitation energy
        mu_2_j is the mass attenuation coefficient of layer j at j excitation energy
        mu_b_d_t is the sum of the products mu * density * thickness of layers between layer i and j
        """
        return self.thisptr.deBoerX(p, q, d1, d2, mu_1_j, mu_2_j, mu_b_d_t)

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
