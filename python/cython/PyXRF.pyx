#import numpy as np
#cimport numpy as np
cimport cython

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

    def setBeam(self, energies, weights, characteristic=None, divergency=None):
        if not hasattr(energies, "__len__"):
            energies = [energies]
        if not hasattr(weights, "__len__"):
            weights = [weights]
        if characteristic is None:
            characteristic = [1] * len(energies)
        if divergency is None:
            divergency = [0.0] * len(energies)

        self._setBeam(energies, weights, characteristic, divergency)


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
        for name, density, thickness, funny in layerList:
            container.push_back(Layer(name, density, thickness, funny))
        self.thisptr.setBeamFilters(container)

    def setSample(self, layerList, referenceLayer=0):
        """
        Due to wrapping constraints, the filter list must have the form:
        [[Material name or formula0, density0, thickness0, funny factor0],
         [Material name or formula1, density1, thickness1, funny factor1],
         ...
         [Material name or formulan, densityn, thicknessn, funny factorn]]

        Unless you know what you are doing, the funny factors must be 1.0
        """
        cdef std_vector[Layer] container
        for name, density, thickness, funny in layerList:
            container.push_back(Layer(name, density, thickness, funny))
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
        for name, density, thickness, funny in layerList:
            container.push_back(Layer(name, density, thickness, funny))
        self.thisptr.setAttenuators(container)
