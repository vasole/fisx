#import numpy as np
#cimport numpy as np
cimport cython

from operator import itemgetter
from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *
from Material cimport *

__doc__ = """

Initialization with XCOM parameters from fisx:

import os
from fisx import DataDir
dataDir = DataDir.DATA_DIR
bindingEnergies = os.path.join(dataDir, "BindingEnergies.dat")
xcomFile = os.path.join(dataDir, "XCOM_CrossSections.dat")
xcom = Elements(dataDir, bindingEnergies, xcomFile)

Initialization with XCOM parameters from PyMca:

import os
from PyMca5 import PyMcaDataDir
dataDir = PyMcaDataDir.PYMCA_DATA_DIR
bindingEnergies = os.path.join(dataDir, "BindingEnergies.dat")
xcomFile = os.path.join(dataDir, "XCOM_CrossSections.dat")
# This is needed because PyMca does not have the non-radiative rates
from fisx import DataDir
dataDir = DataDir.DATA_DIR
xcom = Elements(dataDir, bindingEnergies, xcomFile)

"""
cdef class PyElements:
    cdef Elements *thisptr

    def __cinit__(self, std_string directoryName="",
                        std_string bindingEnergiesFile="",
                        std_string crossSectionsFile=""):
        if len(directoryName) == 0:
            from fisx import DataDir
            directoryName = DataDir.DATA_DIR
        if len(bindingEnergiesFile):
            self.thisptr = new Elements(directoryName, bindingEnergiesFile)
        else:
            self.thisptr = new Elements(directoryName)
        if len(crossSectionsFile):
            self.thisptr.setMassAttenuationCoefficientsFile(crossSectionsFile)

    def initializeAsPyMca(self):
        import os
        try:
            from fisx import DataDir
            directoryName = DataDir.DATA_DIR
            from PyMca5 import PyMcaDataDir
            dataDir = PyMcaDataDir.PYMCA_DATA_DIR
        except ImportError:
            from fisx import DataDir
            directoryName = DataDir.DATA_DIR
            dataDir = directoryName
        bindingEnergies = os.path.join(dataDir, "BindingEnergies.dat")
        xcomFile = os.path.join(dataDir, "XCOM_CrossSections.dat")
        del self.thisptr
        self.thisptr = new Elements(directoryName, bindingEnergies, xcomFile)
        for shell in ["K", "L", "M"]:
            shellConstantsFile = os.path.join(dataDir, shell+"ShellConstants.dat")
            self.thisptr.setShellConstantsFile(shell, shellConstantsFile)

        for shell in ["K", "L", "M"]:
            radiativeRatesFile = os.path.join(dataDir, shell+"ShellRates.dat")
            self.thisptr.setShellRadiativeTransitionsFile(shell, radiativeRatesFile)

    def getElementNames(self):
        return self.thisptr.getElementNames()

    def getComposition(self, std_string materialOrFormula):
        return self.thisptr.getComposition(materialOrFormula)

    def __dealloc__(self):
        del self.thisptr

    def addMaterial(self, PyMaterial material, int errorOnReplace=1):
        self.thisptr.addMaterial(deref(material.thisptr), errorOnReplace)

    def setShellConstantsFile(self, std_string mainShellName, std_string fileName):
        """
        Load main shell (K, L or M) constants from file (fluorescence and Coster-Kronig yields)
        """
        self.thisptr.setShellConstantsFile(mainShellName, fileName)

    def getShellConstantsFile(self, std_string mainShellName):
        return self.thisptr.getShellConstantsFile(mainShellName)

    def setShellRadiativeTransitionsFile(self, std_string mainShellName, std_string fileName):
        """
        Load main shell (K, L or M) X-ray emission rates from file.
        The library normalizes internally.
        """
        self.thisptr.setShellRadiativeTransitionsFile(mainShellName, fileName)

    def getShellRadiativeTransitionsFile(self, std_string mainShellName):
        return self.thisptr.getShellRadiativeTransitionsFile(mainShellName)

    def getShellNonradiativeTransitionsFile(self, std_string mainShellName):
        return self.thisptr.getShellNonradiativeTransitionsFile(mainShellName)

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
    def setMassAttenuationCoefficientsFile(self, std_string crossSectionsFile):
        self.thisptr.setMassAttenuationCoefficientsFile(crossSectionsFile)
    
    def _getSingleMassAttenuationCoefficients(self, std_string element,
                                                     double energy):
        return self.thisptr.getMassAttenuationCoefficients(element, energy)

    def _getElementDefaultMassAttenuationCoefficients(self, std_string element):
        return self.thisptr.getMassAttenuationCoefficients(element)

    def getElementMassAttenuationCoefficients(self, element, energy=None):
        if energy is None:
            return self._getElementDefaultMassAttenuationCoefficients(element)            
        elif hasattr(energy, "__len__"):
            return self._getMultipleMassAttenuationCoefficients(element,
                                                                       energy)
        else:
            return self._getMultipleMassAttenuationCoefficients(element,
                                                                       [energy])

    def _getMultipleMassAttenuationCoefficients(self, std_string element,
                                                       std_vector[double] energy):
        return self.thisptr.getMassAttenuationCoefficients(element, energy)
                                       

    def getMassAttenuationCoefficients(self, name, energy=None):
        if hasattr(name, "keys"):
            return self._getMaterialMassAttenuationCoefficients(name, energy)
        elif energy is None:
            return self._getElementDefaultMassAttenuationCoefficients(name)
        elif hasattr(energy, "__len__"):
            return self._getMultipleMassAttenuationCoefficients(name, energy)
        else:
            # do not use the "single" version to have always the same signature
            return self._getMultipleMassAttenuationCoefficients(name, [energy])

    def getExcitationFactors(self, name, energy, weight=None):
        if hasattr(energy, "__len__"):
            if weight is None:
                weight = [1.0] * len(energy)
            return self._getExcitationFactors(name, energy, weight)[0]
        else:
            energy = [energy]
            if weight is None:
                weight = [1.0]
            else:
                weight = [weight]
            return self._getExcitationFactors(name, energy, weight)

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
        return self.thisptr.getExcitationFactors(element, energies, weights)

    def getPeakFamilies(self, nameOrVector, energy):
        if type(nameOrVector) in [type([]), type(())]:
            return sorted(self._getPeakFamiliesFromVectorOfElements(nameOrVector, energy), key=itemgetter(1))
        else:
            return sorted(self._getPeakFamilies(nameOrVector, energy), key=itemgetter(1))

    def _getPeakFamilies(self, std_string name, double energy):
        return self.thisptr.getPeakFamilies(name, energy)

    def _getPeakFamiliesFromVectorOfElements(self, std_vector[std_string] elementList, double energy):
        return self.thisptr.getPeakFamilies(elementList, energy)

    def getBindingEnergies(self, std_string elementName):
        return self.thisptr.getBindingEnergies(elementName)

    def getEscape(self, std_map[std_string, double] composition, double energy, double energyThreshold=0.010,
                                        double intensityThreshold=1.0e-7,
                                        int nThreshold=4 ,
                                        double alphaIn=90.,
                                        double thickness=0.0):
        return self.thisptr.getEscape(composition, energy, energyThreshold, intensityThreshold, nThreshold,
                                      alphaIn, thickness)

    def getShellConstants(self, std_string elementName, std_string subshell):
        return self.thisptr.getShellConstants(elementName, subshell)

    def getRadiativeTransitions(self, std_string elementName, std_string subshell):
        return self.thisptr.getRadiativeTransitions(elementName, subshell)

    def getNonradiativeTrnasitions(self, std_string elementName, std_string subshell):
        return self.thisptr.getNonradiativeTransitions(elementName, subshell)

    def setElementCascadeCacheEnabled(self, std_string elementName, int flag = 1):
        self.thisptr.setElementCascadeCacheEnabled(elementName, flag)

    def emptyElementCascadeCache(self, std_string elementName):
        self.thisptr.emptyElementCascadeCache(elementName)
