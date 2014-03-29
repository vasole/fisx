#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from XRF cimport *
    
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
