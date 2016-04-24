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

class testXRF(unittest.TestCase):
    def setUp(self):
        """
        import the module
        """
        try:
            from fisx import XRF
            self._module = XRF
        except:
            self._module = None

    def tearDown(self):
        self._module = None

    def testXRFImport(self):
        self.assertTrue(self._module is not None,
                        'Unsuccessful fisx.XRF import')

    def testXRFInstantiation(self):
        try:
            instance = self._module()
        except:
            instance = None
            print("Instantiation error: ",
                    sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2])
        self.assertTrue(instance is not None,
                        'Unsuccesful XRF() instantiation')

    def testXRFResults(self):
        from fisx import Elements
        from fisx import Material
        from fisx import Detector
        from fisx import XRF

        elementsInstance = Elements()
        elementsInstance.initializeAsPyMca()
        # After the slow initialization (to be made once), the rest is fairly fast.
        xrf = XRF()
        xrf.setBeam(16.0) # set incident beam as a single photon energy of 16 keV
        xrf.setBeamFilters([["Al1", 2.72, 0.11, 1.0]]) # Incident beam filters
        # Steel composition of Schoonjans et al, 2012 used to generate table I
        steel = {"C":  0.0445,
                 "N":  0.04,
                 "Si": 0.5093,
                 "P":  0.02,
                 "S":  0.0175,
                 "V":  0.05,
                 "Cr":18.37,
                 "Mn": 1.619,
                 "Fe":64.314, # calculated by subtracting the sum of all other elements
                 "Co": 0.109,
                 "Ni":12.35,
                 "Cu": 0.175,
                 "As": 0.010670,
                 "Mo": 2.26,
                 "W":  0.11,
                 "Pb": 0.001}
        SRM_1155 = Material("SRM_1155", 1.0, 1.0)
        SRM_1155.setComposition(steel)
        elementsInstance.addMaterial(SRM_1155)
        xrf.setSample([["SRM_1155", 1.0, 1.0]]) # Sample, density and thickness
        xrf.setGeometry(45., 45.)               # Incident and fluorescent beam angles
        detector = Detector("Si1", 2.33, 0.035) # Detector Material, density, thickness
        detector.setActiveArea(0.50)            # Area and distance in consistent units
        detector.setDistance(2.1)               # expected cm2 and cm.
        xrf.setDetector(detector)
        Air = Material("Air", 0.0012048, 1.0)
        Air.setCompositionFromLists(["C1", "N1", "O1", "Ar1", "Kr1"],
                                    [0.0012048, 0.75527, 0.23178, 0.012827, 3.2e-06])
        elementsInstance.addMaterial(Air)
        xrf.setAttenuators([["Air", 0.0012048, 5.0, 1.0],
                            ["Be1", 1.848, 0.002, 1.0]]) # Attenuators
        fluo = xrf.getMultilayerFluorescence(["Cr K", "Fe K", "Ni K"],
                                             elementsInstance,
                                             secondary=2,
                                             useMassFractions=1)
        print("\nElement   Peak          Energy       Rate      Secondary  Tertiary")
        for key in fluo:
            for layer in fluo[key]:
                peakList = list(fluo[key][layer].keys())
                peakList.sort()
                for peak in peakList:
                    # energy of the peak
                    energy = fluo[key][layer][peak]["energy"]
                    # expected measured rate
                    rate = fluo[key][layer][peak]["rate"]
                    # primary photons (no attenuation and no detector considered)
                    primary = fluo[key][layer][peak]["primary"]
                    # secondary photons (no attenuation and no detector considered)
                    secondary = fluo[key][layer][peak]["secondary"]
                    # tertiary photons (no attenuation and no detector considered)
                    tertiary = fluo[key][layer][peak].get("tertiary", 0.0)
                    # correction due to secondary excitation
                    enhancement2 = (primary + secondary) / primary
                    enhancement3 = (primary + secondary + tertiary) / primary
                    print("%s   %s    %.4f     %.3g     %.5g    %.5g" % \
                                       (key, peak + (13 - len(peak)) * " ", energy,
                                       rate, enhancement2, enhancement3))
                    # compare against expected values from Schoonjans et al.
                    testXMI = True
                    if (key == "Cr K") and peak.startswith("KL3"):
                        second = 1.626
                        third = 1.671
                    elif (key == "Cr K") and peak.startswith("KM3"):
                        second = 1.646
                        third = 1.694
                    elif (key == "Fe K") and peak.startswith("KL3"):
                        second = 1.063
                        third = 1.064
                    elif (key == "Fe K") and peak.startswith("KL3"):
                        second = 1.065
                        third = 1.066
                    else:
                        testXMI = False
                    if testXMI:
                        discrepancy = 100 * (abs(second-enhancement2)/second)
                        self.assertTrue(discrepancy < 1.5,
                            "%s %s secondary discrepancy = %.1f %%" % \
                            (key, peak, discrepancy))
                        discrepancy = 100 * (abs(third-enhancement3)/third)
                        self.assertTrue(discrepancy < 1.5,
                            "%s %s tertiary discrepancy = %.1f %%" % \
                            (key, peak, discrepancy))

def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(\
            unittest.TestLoader().loadTestsFromTestCase(testXRF))
    else:
        # use a predefined order
        testSuite.addTest(testXRF("testXRFImport"))
        testSuite.addTest(testXRF("testXRFInstantiation"))
        testSuite.addTest(testXRF("testXRFResults"))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
