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

from cython.operator cimport dereference as deref
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from XRF cimport *
from Layer cimport *
    
cdef class PyXRF:
    cdef XRF *thisptr

    def __cinit__(self, std_string configurationFile=""):
        if len(configurationFile):
            self.thisptr = new XRF(configurationFile)
        else:
            self.thisptr = new XRF()

    def __dealloc__(self):
        del self.thisptr

    def readConfigurationFromFile(self, std_string fileName):
        self.thisptr.readConfigurationFromFile(fileName)

    def setBeam(self, energies, weights=None, characteristic=None, divergency=None):
        if not hasattr(energies, "__len__"):
            if divergency is None:
                divergency = 0.0
            self._setSingleEnergyBeam(energies, divergency)
        else:
            if weights is None:
                weights = [1.0] * len(energies)
            elif not hasattr(weights, "__len__"):
                weights = [weights]
            if characteristic is None:
                characteristic = [1] * len(energies)
            if divergency is None:
                divergency = [0.0] * len(energies)

            self._setBeam(energies, weights, characteristic, divergency)

    def _setSingleEnergyBeam(self, double energy, double divergency):
        self.thisptr.setBeam(energy, divergency)

    def _setBeam(self, std_vector[double] energies, std_vector[double] weights, \
                       std_vector[int] characteristic, std_vector[double] divergency):
        self.thisptr.setBeam(energies, weights, characteristic, divergency)

    def setBeamFilters(self, layerList):
        """
        Due to wrapping constraints, the filter list must have the form:
        [[Material name or formula0, density0, thickness0, funny factor0],
         [Material name or formula1, density1, thickness1, funny factor1],
         ...
         [Material name or formulan, densityn, thicknessn, funny factorn]]

        Unless you know what you are doing, the funny factors must be 1.0
        """
        cdef std_vector[Layer] container
        if len(layerList):
            if len(layerList[0]) == 4:
                for name, density, thickness, funny in layerList:
                    container.push_back(Layer(toBytes(name), density, thickness, funny))
            else:
                for name, density, thickness in layerList:
                    container.push_back(Layer(toBytes(name), density, thickness, 1.0))
        self.thisptr.setBeamFilters(container)

    def setSample(self, layerList, referenceLayer=0):
        """
        Due to wrapping constraints, the list must have the form:
        [[Material name or formula0, density0, thickness0, funny factor0],
         [Material name or formula1, density1, thickness1, funny factor1],
         ...
         [Material name or formulan, densityn, thicknessn, funny factorn]]

        Unless you know what you are doing, the funny factors must be 1.0
        """
        cdef std_vector[Layer] container
        if len(layerList[0]) == 4:
            for name, density, thickness, funny in layerList:
                container.push_back(Layer(toBytes(name), density, thickness, funny))
        else:
            for name, density, thickness in layerList:
                container.push_back(Layer(toBytes(name), density, thickness, 1.0))
        self.thisptr.setSample(container, referenceLayer)


    def setAttenuators(self, layerList):
        """
        Due to wrapping constraints, the filter list must have the form:
        [[Material name or formula0, density0, thickness0, funny factor0],
         [Material name or formula1, density1, thickness1, funny factor1],
         ...
         [Material name or formulan, densityn, thicknessn, funny factorn]]

        Unless you know what you are doing, the funny factors must be 1.0
        """
        cdef std_vector[Layer] container
        if len(layerList[0]) == 4:
            for name, density, thickness, funny in layerList:
                container.push_back(Layer(toBytes(name), density, thickness, funny))
        else:
            for name, density, thickness in layerList:
                container.push_back(Layer(toBytes(name), density, thickness, 1.0))
        self.thisptr.setAttenuators(container)

    def setDetector(self, PyDetector detector):
        self.thisptr.setDetector(deref(detector.thisptr))

    def setGeometry(self, double alphaIn, double alphaOut, double scatteringAngle = -90.0):
        if scatteringAngle < 0.0:
            self.thisptr.setGeometry(alphaIn, alphaOut, alphaIn + alphaOut)
        else:
            self.thisptr.setGeometry(alphaIn, alphaOut, scatteringAngle)

    def getMultilayerFluorescence(self, elementFamilyLayer, PyElements elementsLibrary, \
                            int secondary = 0, int useGeometricEfficiency = 1, int useMassFractions = 0, \
                            secondaryCalculationLimit = 0.0):
        """
        Input
        elementFamilyLayer - Vector of strings. Each string represents the information we are interested on.
        "Cr"     - We want the information for Cr, for all line families and sample layers
        "Cr K"   - We want the information for Cr, for the family of K-shell emission lines, in all layers.
        "Cr K 0" - We want the information for Cr, for the family of K-shell emission lines, in layer 0.
        elementsLibrary - Instance of library to be used for all the Physical constants
        secondary - Flag to indicate different levels of secondary excitation to be considered.
                    0 Means not considered
                    1 Consider only intralayer secondary excitation
                    2 Consider intralayer and interlayer secondary excitation
        useGeometricEfficiency - Take into account solid angle or not. Default is 1 (yes)

        useMassFractions - If 0 (default) the output corresponds to the requested information if the mass
        fraction of the element would be one on each calculated sample layer. To get the actual signal, one
        has to multiply bthe rates by the actual mass fraction of the element on each sample layer.
                           If set to 1 the rate will be already corrected by the actual mass fraction.

        Return a complete output of the form
        [Element Family][Layer][line]["energy"] - Energy in keV of the emission line
        [Element Family][Layer][line]["primary"] - Primary rate prior to correct for detection efficiency
        [Element Family][Layer][line]["secondary"] - Secondary rate prior to correct for detection efficiency
        [Element Family][Layer][line]["rate"] - Overall rate
        [Element Family][Layer][line]["efficiency"] - Detection efficiency
        [Element Family][Layer][line][element line layer] - Secondary rate (prior to correct for detection efficiency)
        due to the fluorescence from the given element, line and layer index composing the map key.
        """
        if sys.version > "3.0":
            elementFamilyLayer = [toBytes(x) for x in elementFamilyLayer]
            return toStringKeysAndValues(self.thisptr.getMultilayerFluorescence(elementFamilyLayer, \
                            deref(elementsLibrary.thisptr), \
                            secondary, useGeometricEfficiency, \
                            useMassFractions, secondaryCalculationLimit))
        else:
            return self.thisptr.getMultilayerFluorescence(elementFamilyLayer, \
                            deref(elementsLibrary.thisptr), \
                            secondary, useGeometricEfficiency, \
                            useMassFractions, secondaryCalculationLimit)

    def getFluorescence(self, elementName, PyElements elementsLibrary, \
                            int sampleLayer = 0, lineFamily="K", int secondary = 0, \
                            int useGeometricEfficiency = 1, int useMassFractions = 0, \
                            double secondaryCalculationLimit = 0.0):
        if sys.version > "3.0":
            elementName = toBytes(elementName)
            lineFamily = toBytes(lineFamily)
            return toStringKeysAndValues(self.thisptr.getMultilayerFluorescence(elementName, deref(elementsLibrary.thisptr), \
                            sampleLayer, lineFamily, secondary, useGeometricEfficiency, useMassFractions, \
                            secondaryCalculationLimit))
        else:
            return self.thisptr.getMultilayerFluorescence(elementName, deref(elementsLibrary.thisptr), \
                            sampleLayer, lineFamily, secondary, useGeometricEfficiency, useMassFractions, \
                            secondaryCalculationLimit)

    def getGeometricEfficiency(self, int layerIndex = 0):
        return self.thisptr.getGeometricEfficiency(layerIndex)
