#import numpy as np
#cimport numpy as np
cimport cython

#from libcpp.string cimport string as std_string
#from libcpp.vector cimport vector as std_vector
#from libcpp.map cimport map as std_map

from Shell cimport *

cdef class PyShell:
    cdef Shell *thisptr

    def __cinit__(self, name):
        self.thisptr = new Shell(name)

    def __dealloc__(self):
        del self.thisptr

    def setRadiativeTransitions(self, std_vector[std_string] transitions, std_vector[double] values):
        self.thisptr.setRadiativeTransitions(transitions, values)

    def setNonradiativeTransitions(self, std_vector[std_string] transitions, std_vector[double] values):
        self.thisptr.setNonradiativeTransitions(transitions, values)

    def getAugerRatios(self):
        return self.thisptr.getAugerRatios()

    def getCosterKronigRatios(self):
        return self.thisptr.getCosterKronigRatios()

    def getFluorescenceRatios(self):
        return self.thisptr.getFluorescenceRatios()

    def getRadiativeTransitions(self):
        return self.thisptr.getRadiativeTransitions()

    def getNonradiativeTransitions(self):
        return self.thisptr.getNonradiativeTransitions()

    def getDirectVacancyTransferRatios(self, std_string subshell):
        return self.thisptr.getDirectVacancyTransferRatios(subshell)

    def setShellConstants(self, std_map[std_string, double] shellConstants):
        self.thisptr.setShellConstants(shellConstants)

    def getShellConstants(self):
        return self.thisptr.getShellConstants()
