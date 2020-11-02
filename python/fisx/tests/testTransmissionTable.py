#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2020 European Synchrotron Radiation Facility
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
import numpy

class testTransmissionTable(unittest.TestCase):
    def setUp(self):
        """
        import the module
        """
        try:
            from fisx import TransmissionTable
            self.transmissionTable = TransmissionTable
        except:
            self.transmissionTable = None

    def tearDown(self):
        self.transmissionTable = None

    def testTransmissionTableImport(self):
        self.assertTrue(self.transmissionTable is not None,
                        'Unsuccessful fisx.TransmissionTable import')

    def testTransmissionTableInstantiation(self):
        try:
            instance = self.transmissionTable()
        except:
            instance = None
            print("Instantiation error: ",
                    sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2])
        self.assertTrue(instance is not None,
                        'Unsuccesful TransmissionTable() instantiation')

    def testTransmissionTableDefaults(self):
        instance = self.transmissionTable()
        self.assertTrue(instance.getName() == "",
                        "Incorrect default name")
        self.assertTrue(instance.getComment() == "",
                        "Incorrect default comment")
        transmissionTable = instance.getTransmissionTable() 
        self.assertTrue( transmissionTable == {0.0:1.0},
                        "Incorrect default transmission table")
        t = instance.getTransmission(-1.0)
        self.assertTrue( t == 0.0,
                        "Incorrect transmission %f. Expected 0.0" % t)
        t = instance.getTransmission(1.0)
        self.assertTrue( t == 1.0,
                        "Incorrect transmission %f. Expected 1.0" % t)
        

    def testTransmissionTableResults(self):
        instance = self.transmissionTable()
        self.assertTrue(instance.getName() == "",
                        "Incorrect default name")
        self.assertTrue(instance.getComment() == "",
                        "Incorrect default comment")
        energyList = [0.1, 0.2, 0.3, 0.4, 0.5] 
        tranmissionList = [0.1, 0.2, 0.3, 0.4, 0.5]
        instance.setName("My table")
        self.assertTrue(instance.getName() == "My table",
                        "Incorrect name")
        comment = instance.getComment()
        self.assertTrue(comment == "",
                        "Incorrect default comment got <%s>" % comment)
        instance.setTransmissionTableFromLists(energyList,
                                               transmissionList,
                                               name="",
                                               comment="My comment")
        self.assertTrue(instance.getName() == "My table",
                        "Incorrect name")
        self.assertTrue(instance.getComment() == "My comment",
                        "Incorrect comment")
        v = instance.getTransmission(0.9)
        self.assertTrue(v == 0.0,
                        "Incorrect transmission %f. Expected 0.0" % v)
        v = instance.getTransmission(1.0)
        self.assertTrue(v == 1.0,
                        "Incorrect transmission %f. Expected 1s.0" % v)
        e = [0.15, 0.22, 0.41]
        v = instance.getTransmission(e)
        for i in range(len(e)):
            self.assertTrue(v[i] == e[i],
                            "Incorrect transmission %f expected % f" % \
                                (v[i], e[i]))
        
def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(\
            unittest.TestLoader().loadTestsFromTestCase(testTransmissionTable))
    else:
        # use a predefined order
        for methodName in ["testTransmissionTableImport",
                           "testTransmissionTableInstantiation",
                           "testTransmissionTableDefaults",
                           "testTransmissionTableResults"]:
            testSuite.addTest(testTransmissionTable(methodName))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
