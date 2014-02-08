#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

cdef extern from "fisx_shell.h":
    cdef cppclass Shell:
        Shell(std_string) except +
        void setRadiativeTransitions(std_vector[std_string], std_vector[double]) except +
        std_map[std_string, double] & getShellConstants()
        std_map[std_string, double] & getAugerRatios()
        std_map[std_string, std_map[std_string, double]] & getCosterKronigRatios()
        std_map[std_string, double] & getFluorescenceRatios()
        std_map[std_string, double] & getRadiativeTransitions()
        std_map[std_string, double] & getNonradiativeTransitions()
        std_map[std_string, double] getDirectVacancyTransferRatios(std_string) except +
        void setNonradiativeTransitions(std_vector[std_string], std_vector[double]) except +
        void setShellConstants(std_map[std_string, double]) except +

