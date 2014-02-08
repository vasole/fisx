#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from EPDL97 cimport *
    
cdef class PyEPDL97:
    cdef EPDL97 *thisptr

    def __cinit__(self, std_string name):
        self.thisptr = new EPDL97(name)

    def __dealloc__(self):
        del self.thisptr

    def setDataDirectory(self, std_string name):
        self.thisptr.setDataDirectory(name)

    def setBindingEnergies(self, int z, std_map[std_string, double] energies):
        self.thisptr.setBindingEnergies(z, energies)

    def getBindingEnergies(self, int z):
        return self.thisptr.getBindingEnergies(z)
    
    def getMassAttenuationCoefficients(self, z, energy=None):
        if energy is None:
            return self._getDefaultMassAttenuationCoefficients(z)            
        elif hasattr(energy, "__len__"):
            return self._getMultipleMassAttenuationCoefficients(z, energy)
        else:
            return self._getMultipleMassAttenuationCoefficients(z, [energy])

    def _getDefaultMassAttenuationCoefficients(self, int z):
        return self.thisptr.getMassAttenuationCoefficients(z)

    def _getSingleMassAttenuationCoefficients(self, int z, double energy):
        return self.thisptr.getMassAttenuationCoefficients(z, energy)

    def _getMultipleMassAttenuationCoefficients(self, int z, std_vector[double] energy):
        return self.thisptr.getMassAttenuationCoefficients(z, energy)
                                       
    def getPhotoelectricWeights(self, z, energy):
        if hasattr(energy, "__len__"):
            return self._getMultiplePhotoelectricWeights(z, energy)
        else:
            return self._getMultiplePhotoelectricWeights(z, [energy])

    def _getSinglePhotoelectricWeights(self, int z, double energy):
        return self.thisptr.getPhotoelectricWeights(z, energy)

    def _getMultiplePhotoelectricWeights(self, int z, std_vector[double] energy):
        return self.thisptr.getPhotoelectricWeights(z, energy)
