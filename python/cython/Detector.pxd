#import numpy as np
#cimport numpy as np
cimport cython

from libcpp.string cimport string as std_string
from libcpp.vector cimport vector as std_vector
from libcpp.map cimport map as std_map

from Elements cimport *

cdef extern from "fisx_detector.h":
    cdef cppclass Detector:
        Detector(std_string, double, double, double) except +

        std_vector[double] getTransmission(std_vector[double], Elements, double) except +

        void setActiveArea(double) except +

        void setDiameter(double) except +

        double getActiveArea() except +

        double getDiameter() except +

        void setDistance(double) except +

        double getDistance() except +

        std_map[std_string, std_map[std_string, double]] getEscape(double, Elements, std_string, int) except +

        void setMaximumNumberOfEscapePeaks(int) except +