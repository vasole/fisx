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
#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Detector cimport *
from Elements cimport *
from Layer cimport *

cdef extern from "fisx_xrf.h" namespace "fisx":
    cdef cppclass XRF:
        XRF() except +
        XRF(std_string) except +
        void readConfigurationFromFile(std_string) except +
        void setBeam(std_vector[double], std_vector[double], std_vector[int], std_vector[double]) except +
        void setBeam(double, double) except +
        void setBeamFilters(std_vector[Layer]) except +
        void setSample(std_vector[Layer], int) except +
        void setAttenuators(std_vector[Layer]) except +
        void setGeometry(double, double, double) except +
        void setDetector(Detector) except +
        double getGeometricEfficiency(int) except +

        std_map[std_string, std_map[std_string, double]] getFluorescence(std_string, \
                Elements, int, std_string, int, int, double) except +

        std_map[std_string, std_map[int, std_map[std_string, std_map[std_string, double]]]] \
                getMultilayerFluorescence(std_vector[std_string], Elements, int, int, int, double) except +

        std_map[std_string, std_map[int, std_map[std_string, std_map[std_string, double]]]] \
                getMultilayerFluorescence(std_string, \
                                          Elements, int, std_string, int, int, int, double) except +
