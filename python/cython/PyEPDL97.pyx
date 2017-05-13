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

from EPDL97 cimport *
    
cdef class PyEPDL97:
    cdef EPDL97 *thisptr

    def __cinit__(self, name=None):
        if name is None:
            from fisx import DataDir
            name = DataDir.FISX_DATA_DIR
        self.thisptr = new EPDL97(toBytes(name))

    def __dealloc__(self):
        del self.thisptr

    def setDataDirectory(self, name):
        self.thisptr.setDataDirectory(toBytes(name))

    def setBindingEnergies(self, int z, std_map[std_string, double] energies):
        self.thisptr.setBindingEnergies(z, energies)

    def getBindingEnergies(self, int z):
        return toStringKeys(self.thisptr.getBindingEnergies(z))
    
    def getMassAttenuationCoefficients(self, z, energy=None):
        if energy is None:
            return toStringKeys(self._getDefaultMassAttenuationCoefficients(z))
        elif hasattr(energy, "__len__"):
            return toStringKeys(self._getMultipleMassAttenuationCoefficients(z, energy))
        else:
            return toStringKeys(self._getMultipleMassAttenuationCoefficients(z, [energy]))

    def _getDefaultMassAttenuationCoefficients(self, int z):
        return self.thisptr.getMassAttenuationCoefficients(z)

    def _getSingleMassAttenuationCoefficients(self, int z, double energy):
        return self.thisptr.getMassAttenuationCoefficients(z, energy)

    def _getMultipleMassAttenuationCoefficients(self, int z, std_vector[double] energy):
        return self.thisptr.getMassAttenuationCoefficients(z, energy)
                                       
    def getPhotoelectricWeights(self, z, energy):
        if hasattr(energy, "__len__"):
            return toStringKeys(self._getMultiplePhotoelectricWeights(z, energy))
        else:
            return toStringKeys(self._getMultiplePhotoelectricWeights(z, [energy]))

    def _getSinglePhotoelectricWeights(self, int z, double energy):
        return self.thisptr.getPhotoelectricWeights(z, energy)

    def _getMultiplePhotoelectricWeights(self, int z, std_vector[double] energy):
        return self.thisptr.getPhotoelectricWeights(z, energy)
