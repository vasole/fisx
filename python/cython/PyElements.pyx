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

from operator import itemgetter
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *
from Material cimport *

__doc__ = """

Initialization with XCOM mass attenuation cross sections:

import os
from fisx import DataDir
dataDir = DataDir.FISX_DATA_DIR
bindingEnergies = os.path.join(dataDir, "BindingEnergies.dat")
xcomFile = os.path.join(dataDir, "XCOM_CrossSections.dat")
xcom = Elements(dataDir, bindingEnergies, xcomFile)

"""
cdef class PyElements:
    cdef Elements *thisptr

    def __cinit__(self, directoryName="",
                        bindingEnergiesFile="",
                        crossSectionsFile="",
                        pymca=0):
        if len(directoryName) == 0:
            from fisx import DataDir
            directoryName = DataDir.FISX_DATA_DIR
        directoryName = toBytes(directoryName)
        if pymca:
            pymca = 1
            self.thisptr = new Elements(directoryName, pymca)
        else:
            bindingEnergiesFile = toBytes(bindingEnergiesFile)
            crossSectionsFile = toBytes(crossSectionsFile)
            if len(bindingEnergiesFile):
                self.thisptr = new Elements(directoryName, bindingEnergiesFile, crossSectionsFile)
            else:
                self.thisptr = new Elements(directoryName)
                if len(crossSectionsFile):
                    self.thisptr.setMassAttenuationCoefficientsFile(crossSectionsFile)

    def initializeAsPyMca(self):
        """
        Configure the instance to use the same set of data as PyMca.
        """
        import os
        try:
            from PyMca5 import getDataFile
        except ImportError:
            # old fashion way with duplicated data in PyMca and in fisx
            return self.__initializeAsPyMcaOld()

        from fisx import DataDir
        directoryName = DataDir.FISX_DATA_DIR
        bindingEnergies = getDataFile("BindingEnergies.dat")
        xcomFile = getDataFile("XCOM_CrossSections.dat")
        del self.thisptr
        self.thisptr = new Elements(toBytes(directoryName), toBytes(bindingEnergies), toBytes(xcomFile))
        for shell in ["K", "L", "M"]:
            shellConstantsFile = getDataFile(shell+"ShellConstants.dat")
            self.thisptr.setShellConstantsFile(toBytes(shell),
                                               toBytes(shellConstantsFile))
        for shell in ["K", "L", "M"]:
            radiativeRatesFile = getDataFile(shell+"ShellRates.dat")
            self.thisptr.setShellRadiativeTransitionsFile(toBytes(shell), toBytes(radiativeRatesFile))

    def __initializeAsPyMcaOld(self):
        """
        Configure the instance to use the same set of data as PyMca.
        """
        import os
        try:
            from fisx import DataDir
            directoryName = DataDir.FISX_DATA_DIR
            from PyMca5 import PyMcaDataDir
            dataDir = PyMcaDataDir.PYMCA_DATA_DIR
        except ImportError:
            from fisx import DataDir
            directoryName = DataDir.FISX_DATA_DIR
            dataDir = directoryName
        bindingEnergies = os.path.join(dataDir, "BindingEnergies.dat")
        xcomFile = os.path.join(dataDir, "XCOM_CrossSections.dat")
        del self.thisptr
        self.thisptr = new Elements(toBytes(directoryName), toBytes(bindingEnergies), toBytes(xcomFile))
        for shell in ["K", "L", "M"]:
            shellConstantsFile = os.path.join(dataDir, shell+"ShellConstants.dat")
            self.thisptr.setShellConstantsFile(toBytes(shell), toBytes(shellConstantsFile))

        for shell in ["K", "L", "M"]:
            radiativeRatesFile = os.path.join(dataDir, shell+"ShellRates.dat")
            self.thisptr.setShellRadiativeTransitionsFile(toBytes(shell), toBytes(radiativeRatesFile))

    def getElementNames(self):
        return toStringList(self.thisptr.getElementNames())

    def getMaterialNames(self):
        return toStringList(self.thisptr.getMaterialNames())

    def getComposition(self, materialOrFormula):
        if sys.version < "3.0":
            return self.thisptr.getComposition(toBytes(materialOrFormula))
        else:
            return toStringKeys(self.thisptr.getComposition(toBytes(materialOrFormula)))

    def __dealloc__(self):
        del self.thisptr

    def addMaterial(self, PyMaterial material, int errorOnReplace=1):
        self.thisptr.addMaterial(deref(material.thisptr), errorOnReplace)

    def setShellConstantsFile(self, mainShellName, fileName):
        """
        Load main shell (K, L or M) constants from file (fluorescence and Coster-Kronig yields)
        """
        self.thisptr.setShellConstantsFile(toBytes(mainShellName), toBytes(fileName))

    def getShellConstantsFile(self, mainShellName):
        if sys.version < "3.0":
            return self.thisptr.getShellConstantsFile(mainShellName)
        else:
            return toString(self.thisptr.getShellConstantsFile(toBytes(mainShellName)))

    def setShellRadiativeTransitionsFile(self, mainShellName, fileName):
        """
        Load main shell (K, L or M) X-ray emission rates from file.
        The library normalizes internally.
        """
        self.thisptr.setShellRadiativeTransitionsFile(toBytes(mainShellName), toBytes(fileName))

    def getShellRadiativeTransitionsFile(self, mainShellName):
        if sys.version < "3.0":
            return self.thisptr.getShellRadiativeTransitionsFile(mainShellName)
        else:
            return toString(self.thisptr.getShellRadiativeTransitionsFile(toBytes(mainShellName)))

    def getShellNonradiativeTransitionsFile(self, mainShellName):
        if sys.version < "3.0":
            return self.thisptr.getShellNonradiativeTransitionsFile(mainShellName)
        else:
            return toString(self.thisptr.getShellNonradiativeTransitionsFile(toBytes(mainShellName)))

    def setMassAttenuationCoefficients(self,
                                       std_string element,
                                       std_vector[double] energies,
                                       std_vector[double] photo,
                                       std_vector[double] coherent,
                                       std_vector[double] compton,
                                       std_vector[double] pair):
        self.thisptr.setMassAttenuationCoefficients(element,
                                                    energies,
                                                    photo,
                                                    coherent,
                                                    compton,
                                                    pair)
    def setMassAttenuationCoefficientsFile(self, crossSectionsFile):
        self.thisptr.setMassAttenuationCoefficientsFile(toBytes(crossSectionsFile))
    
    def _getSingleMassAttenuationCoefficients(self, std_string element,
                                                     double energy):
        if sys.version < "3.0":
            return self.thisptr.getMassAttenuationCoefficients(element, energy)
        else:
            return toStringKeys(self.thisptr.getMassAttenuationCoefficients(element, energy))

    def _getElementDefaultMassAttenuationCoefficients(self, std_string element):
        if sys.version < "3.0":
            return self.thisptr.getMassAttenuationCoefficients(element)
        else:
            return toStringKeys(self.thisptr.getMassAttenuationCoefficients(element))

    def getElementMassAttenuationCoefficients(self, element, energy=None):
        if energy is None:
            return self._getElementDefaultMassAttenuationCoefficients(toBytes(element))
        elif hasattr(energy, "__len__"):
            return self._getMultipleMassAttenuationCoefficients(toBytes(element),
                                                                       energy)
        else:
            return self._getMultipleMassAttenuationCoefficients(toBytes(element),
                                                                       [energy])

    def _getMultipleMassAttenuationCoefficients(self, std_string element,
                                                       std_vector[double] energy):
        if sys.version < "3.0":
            return self.thisptr.getMassAttenuationCoefficients(element, energy)
        else:
            return toStringKeys(self.thisptr.getMassAttenuationCoefficients(element, energy))

    def getMassAttenuationCoefficients(self, name, energy=None):
        """
        name can be an element, a formula or a material composition given as a dictionary:
            key is the element name
            fraction is the mass fraction of the element.

        WARNING: The library renormalizes in order to make sure the sum of mass
                 fractions is 1.

        It gives back the mass attenuation coefficients at the given energies as a map where
        the keys are the different physical processes and the values are lists of the 
        calculated values via log-log interpolation in the internal table.
        """
        if hasattr(name, "keys"):
            return self._getMaterialMassAttenuationCoefficients(toBytes(name), energy)
        elif energy is None:
            return self._getElementDefaultMassAttenuationCoefficients(toBytes(name))
        elif hasattr(energy, "__len__"):
            return self._getMultipleMassAttenuationCoefficients(toBytes(name), energy)
        else:
            # do not use the "single" version to have always the same signature
            return self._getMultipleMassAttenuationCoefficients(toBytes(name), [energy])

    def getExcitationFactors(self, name, energy, weight=None):
        """
        getExcitationFactors(name, energy, weight=None)	
        Given energy(s) and (optional) weight(s), for the specfified element, this method returns
        the emitted X-ray already corrected for cascade and fluorescence yield.
        It is the equivalent of the excitation factor in D.K.G. de Boer's paper.
        """
        if hasattr(energy, "__len__"):
            if weight is None:
                weight = [1.0] * len(energy)
            return self._getExcitationFactors(toBytes(name), energy, weight)[0]
        else:
            energy = [energy]
            if weight is None:
                weight = [1.0]
            else:
                weight = [weight]
            return self._getExcitationFactors(toBytes(name), energy, weight)

    def _getMaterialMassAttenuationCoefficients(self, elementDict, energy):
        """
        elementDict is a dictionary of the form:
        elmentDict[key] = fraction where:
            key is the element name
            fraction is the mass fraction of the element.

        WARNING: The library renormalizes in order to make sure the sum of mass
                 fractions is 1.
        """

        if hasattr(energy, "__len__"):
            return self._getMassAttenuationCoefficients(elementDict, energy)
        else:
            return self._getMassAttenuationCoefficients(elementDict, [energy])

    def _getMassAttenuationCoefficients(self, std_map[std_string, double] elementDict,
                                              std_vector[double] energy):
        return self.thisptr.getMassAttenuationCoefficients(elementDict, energy)

    def _getExcitationFactors(self, std_string element,
                                   std_vector[double] energies,
                                   std_vector[double] weights):
        if sys.version < "3.0":
            return self.thisptr.getExcitationFactors(element, energies, weights)
        else:
            return [toStringKeysAndValues(x) for x in self.thisptr.getExcitationFactors(element, energies, weights)]

    def getPeakFamilies(self, nameOrVector, energy):
        """
        getPeakFamilies(nameOrVector, energy)

        Given an energy and a reference to an elements library return dictionarys.
        The key is the peak family ("Si K", "Pb L1", ...) and the value the binding energy.
        """
        if type(nameOrVector) in [type([]), type(())]:
            if sys.version < "3.0":
                return sorted(self._getPeakFamiliesFromVectorOfElements(nameOrVector, energy), key=itemgetter(1))
            else:
                nameOrVector = [toBytes(x) for x in nameOrVector]
                return [(toString(x[0]), x[1]) for x in \
                        sorted(self._getPeakFamiliesFromVectorOfElements(nameOrVector, energy), key=itemgetter(1))]
        else:
            if sys.version < "3.0":
                return sorted(self._getPeakFamilies(toBytes(nameOrVector), energy), key=itemgetter(1))
            else:
                return [(toString(x[0]), x[1]) for x in \
                        sorted(self._getPeakFamilies(toBytes(nameOrVector), energy), key=itemgetter(1))]

    def _getPeakFamilies(self, std_string name, double energy):
        return self.thisptr.getPeakFamilies(name, energy)

    def _getPeakFamiliesFromVectorOfElements(self, std_vector[std_string] elementList, double energy):
        return self.thisptr.getPeakFamilies(elementList, energy)

    def getBindingEnergies(self, elementName):
        if sys.version < "3.0":
            return self.thisptr.getBindingEnergies(elementName)
        else:
            return toStringKeys(self.thisptr.getBindingEnergies(toBytes(elementName)))

    def getEscape(self, composition, double energy, double energyThreshold=0.010,
                                        double intensityThreshold=1.0e-7,
                                        int nThreshold=4 ,
                                        double alphaIn=90.,
                                        double thickness=0.0):
        if sys.version_info < (3, ):
            return self.thisptr.getEscape(toBytesKeys(composition), energy, energyThreshold, intensityThreshold, nThreshold,
                                         alphaIn, thickness)
        else:
            result = toStringKeys(self.thisptr.getEscape(toBytesKeys(composition), energy, energyThreshold,
                                intensityThreshold, nThreshold, alphaIn, thickness))
            keyList =list(result.keys())
            for key in keyList:
                result[key] = toStringKeys(result[key])
            return result

    def updateEscapeCache(self, composition, std_vector[double] energyList, double energyThreshold=0.010,
                                        double intensityThreshold=1.0e-7,
                                        int nThreshold=4 ,
                                        double alphaIn=90.,
                                        double thickness=0.0):
        self.thisptr.updateEscapeCache(toBytesKeys(composition), energyList, energyThreshold, intensityThreshold, nThreshold,
                                      alphaIn, thickness)

    def getShellConstants(self, elementName, subshell):
        if sys.version < "3.0":
            return self.thisptr.getShellConstants(elementName, subshell)
        else:
            return toStringKeys(self.thisptr.getShellConstants(toBytes(elementName), toBytes(subshell)))

    def getEmittedXRayLines(self, elementName, double energy=1000.):
        if sys.version < "3.0":
            return self.thisptr.getEmittedXRayLines(elementName, energy)
        else:
            return toStringKeys(self.thisptr.getEmittedXRayLines(toBytes(elementName), energy))

    def getRadiativeTransitions(self, elementName, subshell):
        if sys.version < "3.0":
            return self.thisptr.getRadiativeTransitions(elementName, subshell)
        else:
            return toStringKeys(self.thisptr.getRadiativeTransitions(toBytes(elementName), toBytes(subshell)))

    def getNonradiativeTransitions(self, elementName, subshell):
        if sys.version < "3.0":
            return self.thisptr.getNonradiativeTransitions(elementName, subshell)
        else:
            return toStringKeys(self.thisptr.getNonradiativeTransitions(toBytes(elementName), toBytes(subshell)))

    def setElementCascadeCacheEnabled(self, elementName, int flag = 1):
        self.thisptr.setElementCascadeCacheEnabled(toBytes(elementName), flag)

    def emptyElementCascadeCache(self, elementName):
        self.thisptr.emptyElementCascadeCache(toBytes(elementName))

    def fillCache(self, elementName, std_vector[double] energy):
        """
        Optimization methods to keep the calculations at a set of energies
        in cache.
        Clear the calculation cache of given element and fill it at the
        selected energies
        """
        self.thisptr.fillCache(toBytes(elementName), energy)

    def updateCache(self, elementName, std_vector[double] energy):
        """
        Update the element cache with those energy values not already present.
        The existing values will be kept.
        """
        self.thisptr.updateCache(toBytes(elementName), energy)

    def setCacheEnabled(self, elementName, int flag = 1):
        """
        Enable or disable the use of the stored calculations (if any).
        It does not clear the cache when disabling.
        """
        self.thisptr.setCacheEnabled(toBytes(elementName), flag)

    def setEscapeCacheEnabled(self, int flag = 1):
        """
        Enable or disable the use of the stored calculations (if any).
        It does not clear the cache when disabling.
        """
        self.thisptr.setEscapeCacheEnabled(flag)

    def clearCache(self, elementName):
        """
        Clear the calculation cache
        """
        self.thisptr.clearCache(toBytes(elementName))

    def isCacheEnabled(self, elementName):
        """
        Return 1 or 0 if the calculation cache is enabled or not
        """
        return self.thisptr.isCacheEnabled(toBytes(elementName))

    def isEscapeCacheEnabled(self):
        """
        Return 1 or 0 if the calculation cache is enabled or not
        """
        return self.thisptr.isEscapeCacheEnabled()

    def getCacheSize(self, elementName):
        """
        Return the number of energies for which the calculations are stored
        """
        return self.thisptr.getCacheSize(toBytes(elementName))

    def removeMaterials(self):
        self.thisptr.removeMaterials()

