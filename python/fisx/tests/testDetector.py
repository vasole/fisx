#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2018 European Synchrotron Radiation Facility
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

class testDetector(unittest.TestCase):
    def setUp(self):
        """
        import the module
        """
        try:
            from fisx import Detector
            self._module = Detector
        except:
            self._module = None

    def tearDown(self):
        self._module = None

    def testDetectorImport(self):
        self.assertTrue(self._module is not None,
                        'Unsuccessful fisx.Detector import')

    def testDetectorInstantiation(self):
        try:
            instance = self._module("Dummy")
            instance = self._module(b"Dummy")
        except:
            instance = None
            print("Instantiation error: ",
                    sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2])
        self.assertTrue(instance is not None,
                        'Unsuccesful Detector() instantiation')

    def testDetectorResults(self):
        from fisx import Elements
        from fisx import Detector

        elementsInstance = Elements()
        elementsInstance.initializeAsPyMca()
        
        # Make a detector from a formula
        detectorName = "NaI"
        detectorDensity = 4.0
        detectorThickness = 0.00350
        detectorInstance = Detector(detectorName,
                                    detectorDensity,
                                    detectorThickness)

        # Check the composition is returned as string and not bytes under
        # Python 3
        composition = detectorInstance.getComposition(elementsInstance)
        if sys.version_info >= (3, ):
            for key in composition:
                self.assertTrue(isinstance(key, str),
                                "Expected string, received %s" % type(key))

        # check the returned keys are correct
        self.assertTrue(len(list(composition.keys())) == 2,
                        "Incorrect number of keys returned")
        for key in ["Na", "I"]:
            self.assertTrue(key in composition,
                            "key %s not found" % key)

        thickness = detectorInstance.getThickness()
        self.assertTrue(abs( thickness - detectorThickness) < 1.0e-7,
                       "Wrong detector thickness!")

        density = detectorInstance.getDensity()
        self.assertTrue(abs(density - detectorDensity) < 1.0e-7,
                       "Wrong detector density!")

        # check the return of the getEscape method
        escape = detectorInstance.getEscape(100.0, elementsInstance)
        if sys.version_info >= (3, ):
            for key in escape:
                self.assertTrue(isinstance(key, str),
                                "Expected string, received %s" % type(key))

def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(\
            unittest.TestLoader().loadTestsFromTestCase(testDetector))
    else:
        # use a predefined order
        testSuite.addTest(testDetector("testDetectorImport"))
        testSuite.addTest(testDetector("testDetectorInstantiation"))
        testSuite.addTest(testDetector("testDetectorResults"))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
