#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Shell cimport *

cdef extern from "fisx_element.h" namespace "fisx":
    cdef cppclass Element:
        Element(std_string, int) except +

        void setName(std_string) except +
        void setAtomicNumber(int)  except +
        int getAtomicNumber()
        
        void setBindingEnergies(std_map[std_string, double])  except +
        # void setBindingEnergies(std_vector[std_string], std_vector[double])
        std_map[std_string, double] & getBindingEnergies()

        void setMassAttenuationCoefficients(std_vector[double],\
                                            std_vector[double],\
                                            std_vector[double],\
                                            std_vector[double],\
                                            std_vector[double]) except +

        std_map[std_string, double] \
            extractEdgeEnergiesFromMassAttenuationCoefficients(std_vector[double],\
                                                                std_vector[double])
        
        void setTotalMassAttenuationCoefficient(std_vector[double],\
                                                std_vector[double]) except +
                                                
        std_map[std_string, std_vector[double]] getMassAttenuationCoefficients() except +
        std_map[std_string, double] getMassAttenuationCoefficients(double) except +
        std_map[std_string, std_vector[double]]\
                            getMassAttenuationCoefficients(std_vector[double]) except +


        void setRadiativeTransitions(std_string ,\
                                     std_vector[std_string],\
                                     std_vector[double])  except +
        
        std_map[std_string, double] getRadiativeTransitions(std_string)  except +

        void setNonradiativeTransitions(std_string subshell,
                                            std_vector[std_string],
                                        std_vector[double])  except +
        std_map[std_string, double] getNonradiativeTransitions(std_string)  except +

        void setShellConstants(std_string, std_map[std_string, double] )  except +
        std_map[std_string, double] getShellConstants(std_string )  except +

        #std_map[std_string, std_map[std_string, double]]\
        #                    getXRayLines(std_string)  except +
        std_map[std_string, std_map[std_string, double]]\
            getXRayLinesFromVacancyDistribution(std_map[std_string, double])  except +

        # HOW TO DO IT??????
        Shell & getShellInstance(std_string)  except +
