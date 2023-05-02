#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2023 European Synchrotron Radiation Facility
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
from TransmissionTable cimport *
from Beam cimport *

cdef class PyXRF:
    cdef XRF *thisptr

    def __cinit__(self, std_string configurationFile=""):
        if len(configurationFile):
            self.thisptr = new XRF(configurationFile)
        else:
            self.thisptr = new XRF()

    def __dealloc__(self):
        del self.thisptr

    def readConfigurationFromFile(self, fileName):
        """
        Read the configuration from a PyMca .cfg ot .fit file
        """
        self.thisptr.readConfigurationFromFile(toBytes(fileName))

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
            if isinstance(layerList[0], PyLayer):
                for layer in layerList:
                    self._addLayerToLayerVector(layer, container)
            elif len(layerList[0]) == 4:
                for name, density, thickness, funny in layerList:
                    container.push_back(Layer(toBytes(name), density, thickness, funny))
            else:
                for name, density, thickness in layerList:
                    container.push_back(Layer(toBytes(name), density, thickness, 1.0))
        self.thisptr.setBeamFilters(container)

    def setUserBeamFilters(self, pyTransmissionTableList):
        """
        Provide a list of already instantiated transmision tables to be used
        as filters between beam and sample
        """
        self._fillTransmissionTable(pyTransmissionTableList, "filter")

    def _fillTransmissionTable(self, pyTransmissionTableList, function):
        if function not in ["filter", "attenuator"]:
            raise ValueError("Please specify usage as filter or as attenuator")
        if len(pyTransmissionTableList):
            if hasattr(pyTransmissionTableList[0], "getTransmissionTable"):
                instantiated = True
            else:
                instantiated = False
        else:
            instantiated = True
        cdef std_vector[TransmissionTable] container
        cdef TransmissionTable t
        if instantiated:
            for item in pyTransmissionTableList:
                name = item.getName()
                comment = item.getComment()
                table = item.getTransmissionTable()
                t = TransmissionTable()
                t.setTransmissionTable(table, name, comment)
                container.push_back(t)
        else:
            for item in pyTransmissionTableList:
                t = TransmissionTable()
                if len(item) == 4:
                    t.setTransmissionTable(item[0],
                                           item[1],
                                           toBytes(item[2]),
                                           toBytes(item[3]))
                elif hasattr(item[0], "keys"):
                    t.setTransmissionTable(item[0],
                                           toBytes(item[1]),
                                           toBytes(item[2]))
                else:
                    raise ValueError("Not appropriate input type or length")
                container.push_back(t)
        if function == "filter":
            self.thisptr.setUserBeamFilters(container)
        else:
            self.thisptr.setUserAttenuators(container)

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
        if len(layerList):
            if isinstance(layerList[0], PyLayer):
                for layer in layerList:
                    self._addLayerToLayerVector(layer, container)
            elif len(layerList[0]) == 4:
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
        if len(layerList):
            if isinstance(layerList[0], PyLayer):
                for layer in layerList:
                    self._addLayerToLayerVector(layer, container)
            elif len(layerList[0]) == 4:
                for name, density, thickness, funny in layerList:
                    container.push_back(Layer(toBytes(name), density, thickness, funny))
            else:
                for name, density, thickness in layerList:
                    container.push_back(Layer(toBytes(name), density, thickness, 1.0))
        return self.thisptr.setAttenuators(container)

    cdef void _addLayerToLayerVector(self, PyLayer layer, std_vector[Layer]& container):
        container.push_back(deref(layer.thisptr))

    def setUserAttenuators(self, pyTransmissionTableList):
        """
        Provide a list of in which each item is either an already instantiated
        PyTransmissionTable or each item is a list of the arguments of the
        PyTransmissionTable method setTransmissionTable.
        This transmission tables will be used as filters between sample and
        detector
        """
        self._fillTransmissionTable(pyTransmissionTableList, "attenuator")

    def setDetector(self, PyDetector detector):
        self.thisptr.setDetector(deref(detector.thisptr))

    def setGeometry(self, double alphaIn, double alphaOut, double scatteringAngle = -90.0):
        if scatteringAngle < 0.0:
            self.thisptr.setGeometry(alphaIn, alphaOut, alphaIn + alphaOut)
        else:
            self.thisptr.setGeometry(alphaIn, alphaOut, scatteringAngle)

    def getLayerComposition(self, PyLayer layerInstance, PyElements elementsLibrary):
        return toStringKeys(self.thisptr.getLayerComposition(deref(layerInstance.thisptr),
                                                deref(elementsLibrary.thisptr)))

    def getLayerMassAttenuationCoefficients(self, PyLayer layerInstance, energies, \
                                            PyElements elementsLibrary):
        if not hasattr(energies, "__len__"):
            return toStringKeys(self._getLayerMassAttenuationCoefficientsSingle( \
                                    layerInstance, \
                                    energies, elementsLibrary))
        else:
            return toStringKeys(self._getLayerMassAttenuationCoefficientsMultiple( \
                                    layerInstance,
                                    energies,
                                    elementsLibrary))

    def _getLayerMassAttenuationCoefficientsSingle(self, PyLayer layerInstance, \
                                            double energy, \
                                            PyElements elementsLibrary):
        cdef std_map[std_string, double] composition
        composition.clear()
        return self.thisptr.getLayerMassAttenuationCoefficients( \
                                    deref(layerInstance.thisptr), \
                                    energy, deref(elementsLibrary.thisptr), \
                                    composition)

    def _getLayerMassAttenuationCoefficientsMultiple(self, PyLayer layerInstance, \
                                            std_vector[double] energies, \
                                            PyElements elementsLibrary):
        cdef std_map[std_string, double] composition
        composition.clear()
        return self.thisptr.getLayerMassAttenuationCoefficients( \
                                    deref(layerInstance.thisptr), \
                                    energies, deref(elementsLibrary.thisptr), \
                                    composition)

    def getLayerTransmission(self, PyLayer layerInstance, energies, \
                             PyElements elementsLibrary, angle=90.):
        if not hasattr(energies, "__len__"):
            return self._getLayerTransmissionSingle( \
                                    layerInstance, \
                                    energies, \
                                    elementsLibrary, \
                                    angle)
        else:
            return self._getLayerTransmissionMultiple( \
                                    layerInstance, \
                                    energies, \
                                    elementsLibrary, \
                                    angle)

    def _getLayerTransmissionSingle(self, PyLayer layerInstance, double energy, \
                             PyElements elementsLibrary, double angle):
        cdef std_map[std_string, double] composition
        composition.clear()
        return self.thisptr.getLayerTransmission( \
                                    deref(layerInstance.thisptr), \
                                    energy, \
                                    deref(elementsLibrary.thisptr), \
                                    angle, \
                                    composition)

    def _getLayerTransmissionMultiple(self, PyLayer layerInstance, \
                             std_vector[double] energies, \
                             PyElements elementsLibrary, double angle):
        cdef std_map[std_string, double] composition
        composition.clear()
        return self.thisptr.getLayerTransmission( \
                                    deref(layerInstance.thisptr), \
                                    energies, \
                                    deref(elementsLibrary.thisptr), \
                                    angle, \
                                    composition)

    def getLayerPeakFamilies(self, PyLayer layerInstance, double energy, \
                             PyElements elementsLibrary):
        cdef std_map[std_string, double] composition
        composition.clear()
        return toStringKeys(self.thisptr.getLayerPeakFamilies( \
                                    deref(layerInstance.thisptr), \
                                    energy, deref(elementsLibrary.thisptr), \
                                    composition))

    def getMultilayerFluorescence(self, elementFamilyLayer, PyElements elementsLibrary, \
                            int secondary = 0, int useGeometricEfficiency = 1, int useMassFractions = 0, \
                            double secondaryCalculationLimit = 0.0):
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
        cdef std_vector[std_string] elementFamilyLayerVector
        cdef std_map[std_string, std_map[int, std_map[std_string, std_map[std_string, double]]]] result
        if sys.version > "3.0":
            for x in elementFamilyLayer:
                elementFamilyLayerVector.push_back(toBytes(x))
            with nogil:
                result = self.thisptr.getMultilayerFluorescence( \
                            elementFamilyLayerVector, \
                            deref(elementsLibrary.thisptr), \
                            secondary, useGeometricEfficiency, \
                            useMassFractions, secondaryCalculationLimit)
            return toStringKeysAndValues(result)
        else:
            return self.thisptr.getMultilayerFluorescence(elementFamilyLayer, \
                            deref(elementsLibrary.thisptr), \
                            secondary, useGeometricEfficiency, \
                            useMassFractions, secondaryCalculationLimit)

    def getFluorescence(self, elementNames, PyElements elementsLibrary, \
                            sampleLayer = 0, lineFamily="K", int secondary = 0, \
                            int useGeometricEfficiency = 1, int useMassFractions = 0, \
                            double secondaryCalculationLimit = 0.0,
                            beam = None):

        """
        Input
        elementNames - Single string or Vector of strings. Each string represents the information we are interested on.
        "Cr"         - We want the information for Cr
        ["Cr", "Fe"] - We want the information for Cr and Fe.
        elementsLibrary - Instance of library to be used for all the Physical constants
        sampleLayer - Single integer or vector of integers representing the layers where the calculation has to take place.
                      A negative value implies the calculation will take places in all layers
                      0 Means calculation on top layer, 1 the second, and so on
                      The program expects either a single index or as many indices as element names provided.
        lineFamily - Single string or Vector of strings representing the peak families to calculate.
                     The program expects as many peak families as element names provided.
        secondary  - Flag to indicate different levels of secondary excitation to be considered.
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
        cdef std_vector[std_string] elementNamesVector
        cdef std_vector[int] sampleLayerIndicesVector
        cdef std_vector[std_string] lineFamiliesVector
        cdef std_map[std_string, std_map[int, std_map[std_string, std_map[std_string, double]]]] result
        cdef Beam beamInstance = Beam();
        cdef std_vector[double] beamEnergies
        cdef std_vector[double] beamWeights
        cdef std_vector[int] dummyIntVec
        cdef std_vector[double] dummyDoubleVec

        if hasattr(elementNames[0], "__len__"):
            # we have received a list of elements
            pass
        else:
            # we have a single element, convert to list
            elementNames = [elementNames]

        if hasattr(sampleLayer, "__len__"):
            # we have received a list of layer indices
            pass
        else:
            # we should have received an integer, convert to list
            sampleLayer = [sampleLayer]

        if len(lineFamily) and hasattr(lineFamily, "__len__"):
            # we have received a list of peak families
            pass
        else:
            # we should have a peak family, convert to list
            lineFamily = [lineFamily]

        if beam:
            # for the time being only a list expected
            if not hasattr(beam, "__len__"):
                # assume a single energy
                beamEnergies.push_back(beam)
            elif len(beam):
                if hasattr(beam[0], "__len__"):
                    # at the very least have energies and weights
                    for item in beam:
                        beamEnergies.push_back(item[0])
                        beamWeights.push_back(item[1])
                else:
                    # we have a list of energies
                    for item in beam:
                        beamEnergies.push_back(item)
            beamInstance.setBeam(beamEnergies, beamWeights, dummyIntVec, dummyDoubleVec)

        # check the sizes match
        if len(lineFamily) != len(elementNames):
            raise IndexError("Number of elements should match the number of requested peak families")

        if len(sampleLayer) > 1:
            if len(sampleLayer) != len(elementNames):
                raise IndexError("Provide a single layer index or as many indices as elements")

        # fill the vectors
        for x in elementNames:
            elementNamesVector.push_back(toBytes(x))

        for x in sampleLayer:
            sampleLayerIndicesVector.push_back(x)

        for x in lineFamily:
            lineFamiliesVector.push_back(toBytes(x))

        with nogil:
            result = self.thisptr.getMultilayerFluorescence(elementNamesVector, deref(elementsLibrary.thisptr), \
                            sampleLayerIndicesVector, lineFamiliesVector, secondary, useGeometricEfficiency, useMassFractions, \
                            secondaryCalculationLimit, beamInstance)

        return toStringKeysAndValues(result)

    def getGeometricEfficiency(self, int layerIndex = 0):
        return self.thisptr.getGeometricEfficiency(layerIndex)
