#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2016 European Synchrotron Radiation Facility
#
# This file is part of the fisx X-ray developed by V.A. Sole
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
__author__ = "V.A. Sole - ESRF Data Analysis"
import unittest
import sys
import os
import tempfile

import numpy
import time
import os
import sys

try:
    from PyMca5.PyMca import PyMcaEPDL97
    PYMCA = True
except ImportError:
    # do not compare with PyMca because it is not installed
    PYMCA = False

ElementList= ['H', 'He', 
            'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe',
            'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
            'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
            'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
            'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce',
            'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 
            'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 
            'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 
            'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 
            'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 
            'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 
            'Bh', 'Hs', 'Mt']

def getSymbol(z):
    return ElementList[z-1]

def getZ(ele):
    return ElementList.index(ele) + 1

class testEPDL97(unittest.TestCase):
    def setUp(self):
        """
        import the module
        """
        try:
            from fisx import EPDL97
            self.epdl97 = EPDL97
        except:
            self.edpl97 = None

    def tearDown(self):
        self.epdl97 = None

    def testEPDL97Import(self):
        self.assertTrue(self.epdl97 is not None,
                        'Unsuccessful fisx.EPDL97 import')
    if PYMCA:
        def testEPDL97BindingVersusPyMcaEPDL97Binding(self):
            from fisx import DataDir
            dirname = DataDir.FISX_DATA_DIR
            epdl = self.epdl97(dirname)
            for i in range(1, 99):
                pymca = PyMcaEPDL97.EPDL97_DICT[getSymbol(i)]['binding']
                cpp = epdl.getBindingEnergies(i)
                for key in cpp:
                    difference = abs(cpp[key] - pymca[key])
                    self.assertTrue(difference < 1.0e-7,
                        "Element %s, shell %s , difference %f" %\
                            (getSymbol(i), key, difference))
            epdl = None

        def testEPDL97MuVersusPyMcaEPDL97Mu(self):
            from fisx import DataDir
            dirname = DataDir.FISX_DATA_DIR
            epdl = self.epdl97(dirname)
            x = numpy.linspace(1.0, 80., 157)
            for j in range(1, 100):
                pymca = PyMcaEPDL97.getElementCrossSections(getSymbol(j), x)
                cpp = epdl.getMassAttenuationCoefficients(j, x)
                for key in cpp:
                    for i in range(len(x)):
                        if key == "photoelectric":
                            delta = cpp[key][i] - pymca["photo"][i]
                        else:
                            delta = cpp[key][i] - pymca[key][i]
                        if cpp[key][i] > 0:
                            delta = 100. * (abs(delta)/cpp[key][i])
                            tol = 1.0e-5
                        else:
                            delta = abs(delta)
                            tol = 1.0e-7
                        self.assertTrue(delta < tol,
                            "z = %d, effect = %s, energy = %f, delta = %f" %\
                                (j, key, x[i], delta))
            epdl = None

        def testEPDL97PartialVersusPyMcaPartial(self):
            from PyMca5.PyMca import PyMcaEPDL97
            from fisx import DataDir
            dirname = DataDir.FISX_DATA_DIR
            epdl = self.epdl97(dirname)
            pymca = PyMcaEPDL97.getElementCrossSections('Pb')
            idx = numpy.nonzero((pymca['energy'] > 1.0) & \
                                (pymca['energy'] < 100))[0]
            x = pymca['energy'][idx]
            # test partial photoelectric cross sections
            shellList = ["K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"]
            for j in range(1, 100):
                cpp = epdl.getPhotoelectricWeights(j, x)
                for i in range(len(x)):
                    pymca = PyMcaEPDL97.getPhotoelectricWeights(getSymbol(j),
                                                                    shellList,
                                                                    x[i])
                    for ikey in range(len(shellList)):
                        key = shellList[ikey]
                        delta = abs(pymca[ikey]- cpp[key][i])
                        tol = 1.0e-7
                        self.assertTrue(delta < tol,
                            "Default E, z = %d, shell = %s, energy = %f, delta = %f" %\
                                (j, key, x[i], delta))
            print("\nPhotoelectric weights OK at default energies")
            x += 0.123
            for j in range(1, 100):
                cpp = epdl.getPhotoelectricWeights(j, x)
                for i in range(len(x)):
                    if x[i] < 1.0:
                        continue
                    pymca = PyMcaEPDL97.getPhotoelectricWeights(getSymbol(j),
                                                                shellList,
                                                                x[i])
                    for ikey in range(len(shellList)):
                        key = shellList[ikey]
                        delta = abs(pymca[ikey]- cpp[key][i])
                        tol = 1.0e-7
                        self.assertTrue(delta < tol,
                            "z = %d, shell = %s, energy = %f, delta = %f" %\
                                (j, key, x[i], delta))
            print("Photoelectric weights OK at any energy")

def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(\
            unittest.TestLoader().loadTestsFromTestCase(testEPDL97))
    else:
        # use a predefined order
        testSuite.addTest(testEPDL97("testEPDL97Import"))
        if PYMCA:
            testSuite.addTest(testEPDL97("testEPDL97BindingVersusPyMcaEPDL97Binding"))
            testSuite.addTest(\
                testEPDL97("testEPDL97MuVersusPyMcaEPDL97Mu"))
            testSuite.addTest(\
                testEPDL97("testEPDL97PartialVersusPyMcaPartial"))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
