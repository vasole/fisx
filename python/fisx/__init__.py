from ._fisx import PySimpleIni as SimpleIni
from ._fisx import PySimpleSpecfile as SimpleSpecfile
from ._fisx import PyEPDL97 as EPDL97
from ._fisx import PyShell as Shell
from ._fisx import PyElement as Element
from ._fisx import PyElements as Elements
from ._fisx import PyLayer as Layer
from ._fisx import PyDetector as Detector
from ._fisx import PyXRF as XRF
from ._fisx import PyMath as Math
from ._fisx import PyMaterial as Material
from ._fisx import fisxVersion

__version__ = fisxVersion()

def version():
    return __version__
