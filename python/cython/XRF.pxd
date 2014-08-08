#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Detector cimport *
from Elements cimport *
from Layer cimport *

cdef extern from "fisx_xrf.h":
    cdef cppclass XRF:
        XRF() except +
        XRF(std_string) except +
        void readConfigurationFromFile(std_string) except +
        void setBeam(std_vector[double], std_vector[double], std_vector[int], std_vector[double]) except +
        void setBeam(double, double) except +
        void setBeamFilters(std_vector[Layer]) except +
        void setSample(std_vector[Layer], int) except +
        void setAttenuators(std_vector[Layer]) except +
        void setGeometry(double, double, double) except +
        void setDetector(Detector) except +
        double getGeometricEfficiency(int) except +

        std_map[std_string, std_map[std_string, double]] getFluorescence(std_string, \
                Elements, int, std_string, int, int) except +

        std_map[std_string, std_map[int, std_map[std_string, std_map[std_string, double]]]] \
                getMultilayerFluorescence(std_vector[std_string], Elements, int, int, int) except +

        std_map[std_string, std_map[int, std_map[std_string, std_map[std_string, double]]]] \
                getMultilayerFluorescence(std_string, \
                                          Elements, int, std_string, int, int, int) except +
