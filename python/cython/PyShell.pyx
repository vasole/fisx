import sys
cimport cython

#from libcpp.string cimport string as std_string
#from libcpp.vector cimport vector as std_vector
#from libcpp.map cimport map as std_map

from Shell cimport *

cdef class PyShell:
    cdef Shell *thisptr

    def __cinit__(self, name):
        name = toBytes(name)
        self.thisptr = new Shell(name)

    def __dealloc__(self):
        del self.thisptr

    def setRadiativeTransitions(self, transitions, std_vector[double] values):
        if sys.version > "3.0":
            transitions = [toBytes(x) for x in transitions]
        self.thisptr.setRadiativeTransitions(transitions, values)

    def setNonradiativeTransitions(self, transitions, std_vector[double] values):
        if sys.version > "3.0":
            transitions = [toBytes(x) for x in transitions]
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

    def getDirectVacancyTransferRatios(self, subshell):
        return self.thisptr.getDirectVacancyTransferRatios(toBytes(subshell))

    def setShellConstants(self, shellConstants):
        if sys.version > "3.0":
            shellConstants = toBytesKeys(shellConstants)
        self.thisptr.setShellConstants(shellConstants)

    def getShellConstants(self):
        return self.thisptr.getShellConstants()
