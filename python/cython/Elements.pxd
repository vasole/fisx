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
#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map
from libcpp.pair cimport pair as std_pair

from Material cimport *

cdef extern from "fisx_elements.h" namespace "fisx":
    cdef cppclass Elements:
        Elements(std_string) except +
        Elements(std_string, short) except +
        Elements(std_string, std_string, std_string) except +

        std_vector[std_string] getElementNames()

        std_vector[std_string] getMaterialNames()

        std_map[std_string, double] getComposition(std_string) except +

        void addMaterial(Material, int) except +

        void setShellConstantsFile(std_string, std_string) except +

        void setShellRadiativeTransitionsFile(std_string, std_string) except +

        void setShellNonradiativeTransitionsFile(std_string, std_string) except +

        std_string getShellConstantsFile(std_string) except +

        std_string getShellRadiativeTransitionsFile(std_string) except +

        std_string getShellNonradiativeTransitionsFile(std_string) except +

        void setMassAttenuationCoefficients(std_string,\
                                            std_vector[double], \
                                            std_vector[double], \
                                            std_vector[double], \
                                            std_vector[double], \
                                            std_vector[double]) except +

        void setMassAttenuationCoefficientsFile(std_string) except +

        std_map[std_string, double] getMassAttenuationCoefficients(std_string, double) except +

        std_map[std_string, std_vector[double]] getMassAttenuationCoefficients(std_string) except +

        std_map[std_string, std_vector[double]]\
                            getMassAttenuationCoefficients(std_string, std_vector[double]) except +
        
        std_map[std_string, std_vector[double]] getMassAttenuationCoefficients(std_map[std_string, double],\
                                                                               double) except +

        std_map[std_string, std_vector[double]]\
                            getMassAttenuationCoefficients(std_map[std_string, double],\
                                                           std_vector[double]) except +

        std_vector[std_map[std_string, std_map[std_string, double]]] getExcitationFactors( \
                            std_string element,
                            std_vector[double] energy,
                            std_vector[double] weights) except +

        std_vector[std_pair[std_string, double]] getPeakFamilies(std_string, double) except +

        std_vector[std_pair[std_string, double]] getPeakFamilies(std_vector[std_string], double) except +

        std_map[std_string, double] getBindingEnergies(std_string) except +

        std_map[std_string, std_map [std_string, double]] getEscape(std_map[std_string, double],\
                                              double, double, double, int, double, double) except +

        void updateEscapeCache(std_map[std_string, double],\
                                       std_vector[double], double, double, int, double, double) except +
        
        std_map[std_string, double] getEmittedXRayLines(std_string, double) except +

        std_map[std_string, double] getRadiativeTransitions(std_string, std_string) except +

        std_map[std_string, double] getNonradiativeTransitions(std_string, std_string) except +

        std_map[std_string, double] getShellConstants(std_string, std_string) except +

        void setElementCascadeCacheEnabled(std_string, int) except +

        void emptyElementCascadeCache(std_string) except +

        void fillCache(std_string, std_vector[double]) except +

        void updateCache(std_string, std_vector[double]) except +

        void setCacheEnabled(std_string, int) except +

        void clearCache(std_string) except +

        int isCacheEnabled(std_string) except +

        int isEscapeCacheEnabled() except +

        void setEscapeCacheEnabled(int) except+

        int getCacheSize(std_string) except+

        void removeMaterials()
