#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *
    
cdef class PyElements:
    cdef Elements *thisptr

    def __cinit__(self, std_string directoryName):
        self.thisptr = new Elements(directoryName)

    def __dealloc__(self):
        del self.thisptr

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
    
    def _getElementSingleMassAttenuationCoefficients(self, std_string element,
                                                     double energy):
        return self.thisptr.getMassAttenuationCoefficients(element, energy)

    def _getElementDefaultMassAttenuationCoefficients(self, std_string element):
        return self.thisptr.getMassAttenuationCoefficients(element)

    def getElementMassAttenuationCoefficients(self, element, energy=None):
        if energy is None:
            return self._getElementDefaultMassAttenuationCoefficients(element)            
        elif hasattr(energy, "__len__"):
            return self._getElementMultipleMassAttenuationCoefficients(element,
                                                                       energy)
        else:
            return self._getElementMultipleMassAttenuationCoefficients(element,
                                                                       [energy])

    def _getElementMultipleMassAttenuationCoefficients(self, std_string element,
                                                       std_vector[double] energy):
        return self.thisptr.getMassAttenuationCoefficients(element, energy)
                                       

    def getMassAttenuationCoefficients(self, elementDict, energy):
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
