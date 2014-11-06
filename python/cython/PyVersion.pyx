cimport cython

#from libcpp.string cimport string as std_string

from Version cimport fisxVersion as _fisxVersion

def fisxVersion():
    return _fisxVersion()