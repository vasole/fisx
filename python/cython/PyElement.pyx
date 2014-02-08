#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Element cimport *
    
cdef class PyElement:
    cdef Element *thisptr

    def __cinit__(self, std_string name, z = 0):
        self.thisptr = new Element(name, z)

    def __dealloc__(self):
        del self.thisptr

    def setName(self, std_string name):
        self.thisptr.setName(name)

    def setAtomicNumber(self, int z):
        self.thisptr.setAtomicNumber(z)

    def getAtomicNumber(self):
        return self.thisptr.getAtomicNumber()

    def setBindingEnergies(self, std_map[std_string, double] energies):
        self.thisptr.setBindingEnergies(energies)

    def getBindingEnergies(self):
        return self.thisptr.getBindingEnergies()
    
    def setMassAttenuationCoefficients(self,
                                       std_vector[double] energies,
                                       std_vector[double] photo,
                                       std_vector[double] coherent,
                                       std_vector[double] compton,
                                       std_vector[double] pair):
        self.thisptr.setMassAttenuationCoefficients(energies,
                                                    photo,
                                                    coherent,
                                                    compton,
                                                    pair)
    
    def _getDefaultMassAttenuationCoefficients(self):
        return self.thisptr.getMassAttenuationCoefficients()

    def _getSingleMassAttenuationCoefficients(self, double energy):
        return self.thisptr.getMassAttenuationCoefficients(energy)

    def getMassAttenuationCoefficients(self, energy=None):
        if energy is None:
            return self._getDefaultMassAttenuationCoefficients()
        elif hasattr(energy, "__len__"):
            return self._getMultipleMassAttenuationCoefficients(energy)
        else:
            return self._getMultipleMassAttenuationCoefficients([energy])

    def _getMultipleMassAttenuationCoefficients(self, std_vector[double] energy):
        return self.thisptr.getMassAttenuationCoefficients(energy)
                                       
    def setRadiativeTransitions(self, std_string shell,
                                std_vector[std_string] labels,
                                std_vector[double] values):
        self.thisptr.setRadiativeTransitions(shell, labels, values)

    def getRadiativeTransitions(self, std_string shell):
        return self.thisptr.getRadiativeTransitions(shell)

    def setNonradiativeTransitions(self, std_string shell,
                                   std_vector[std_string] labels,
                                   std_vector[double] values):
        self.thisptr.setNonradiativeTransitions(shell, labels, values)

    def getNonradiativeTransitions(self, std_string shell):
        return self.thisptr.getNonradiativeTransitions(shell)

    def setShellConstants(self, std_string shell,
                          std_map[std_string, double] valuesDict):
        self.thisptr.setShellConstants(shell, valuesDict)

    def getShellConstants(self, std_string shell):
        return self.thisptr.getShellConstants(shell)

    #def getXRayLines(self, std_string shell):
    #    return self.thisptr.getXRayLines(shell)
                          
    def getXRayLinesFromVacancyDistribution(self,
                            std_map[std_string, double] vacancyDict):
        return self.thisptr.getXRayLinesFromVacancyDistribution(\
                                vacancyDict)
    
