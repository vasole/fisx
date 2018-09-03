#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2018 European Synchrotron Radiation Facility
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

class testElements(unittest.TestCase):
    def setUp(self):
        """
        import the module
        """
        try:
            from fisx import Elements
            self.elements = Elements
        except:
            self.elements = None

    def tearDown(self):
        self.elements = None

    def testElementsImport(self):
        self.assertTrue(self.elements is not None,
                        'Unsuccessful fisx.Elements import')

    def testElementsInstantiation(self):
        try:
            elementsInstance = self.elements()
        except:
            elementsInstance = None
            print("Instantiation error: ",
                    sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2])
        self.assertTrue(elementsInstance is not None,
                        'Unsuccesful Elements() instantiation')

    def testElementsResults(self):
        elementsInstance = self.elements()
        elementsInstance.initializeAsPyMca()
        
        elementNames = elementsInstance.getElementNames()
        if sys.version_info >= (3,):
            for item in elementNames:
                self.assertTrue(isinstance(item, str),
                                "getElementNames expected str received %s" % \
                                type(item))
            
        for item in ElementList:
            self.assertTrue(item in elementNames,
                            "Element %s not found in library" % item)

        # Check the composition is returned as string and not bytes under
        # Python 3
        composition = elementsInstance.getComposition("NaI")
        if sys.version_info >= (3,):
            for key in composition:
                self.assertTrue(isinstance(key, str),
                                "Expected string, received %s" % type(key))

        # check the returned keys are correct
        self.assertTrue(len(list(composition.keys())) == 2,
                        "Incorrect number of keys returned")
        for key in ["Na", "I"]:
            self.assertTrue(key in composition,
                            "key %s not found" % key)

        for shell in ["K", "L", "M"]:
            fileName = elementsInstance.getShellRadiativeTransitionsFile(shell)
            self.assertTrue(isinstance(fileName, str),
                            "Expected string, received %s" % type(fileName))

        for shell in ["K", "L", "M"]:
            fileName = elementsInstance.getShellNonradiativeTransitionsFile(shell)
            self.assertTrue(isinstance(fileName, str),
                            "Expected string, received %s" % type(fileName))

        escapeDict = elementsInstance.getEscape(composition, 100.)
        if sys.version_info >= (3,):
            for key in escapeDict:
                self.assertTrue(isinstance(key, str),
                        "getEscape expected string key, received %s" % type(key))
                if isinstance(escapeDict[key], dict):
                    for key2 in escapeDict[key]:
                        self.assertTrue(isinstance(key2, str),
                            "Expected string subkey, received %s" % type(key2))

def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(\
            unittest.TestLoader().loadTestsFromTestCase(testElements))
    else:
        # use a predefined order
        testSuite.addTest(testElements("testElementsImport"))
        testSuite.addTest(testElements("testElementsInstantiation"))
        testSuite.addTest(testElements("testElementsResults"))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
