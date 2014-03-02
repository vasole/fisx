__author__ = "V.A. Sole - ESRF Data Analysis"
import unittest
import os
import sys

class testData(unittest.TestCase):
    ELEMENTS = ['H', 'He', 
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

    def setUp(self):
        """
        Get the data directory
        """
        self._importSuccess = False
        try:
            from fisx import DataDir
            self._importSuccess = True
            self.dataDir = DataDir.DATA_DIR
        except:
            self.dataDir = None

        try:
            self.docDir = DataDir.DOC_DIR
        except:
            self.docDir = None

    def testDataDirectoryPresence(self):
        self.assertTrue(self._importSuccess,
                        'Unsuccessful fisx.DataDir import')
        self.assertTrue(self.dataDir is not None,
                        'Unassigned fisx.DataDir.DATA_DIR')
        self.assertTrue(os.path.exists(self.dataDir),
                        'Directory "%s" does not exist' % self.dataDir)
        self.assertTrue(os.path.isdir(self.dataDir),
                        '"%s" expected to be a directory' % self.dataDir)

    def testDataFilePresence(self):
        # Testing file presence
        self.testDataDirectoryPresence()
        for fname in ['BindingEnergies.dat',
                      'EADL97_BindingEnergies.dat',
                      'EADL97_KShellConstants.dat',
                      'EADL97_KShellNonradiativeRates.dat',
                      'EADL97_KShellRadiativeRates.dat',
                      'EADL97_LShellConstants.dat',
                      'EADL97_LShellNonradiativeRates.dat',
                      'EADL97_LShellRadiativeRates.dat',
                      'EADL97_MShellConstants.dat',
                      'EADL97_MShellNonradiativeRates.dat',
                      'EADL97_MShellRadiativeRates.dat',
                      'EPDL97_CrossSections.dat',
                      'KShellRates.dat',
                      'LShellRates.dat',
                      'MShellRates.dat',
                      'XCOM_CrossSections.dat']:
            actualName = os.path.join(self.dataDir, fname)
            self.assertTrue(os.path.exists(actualName),
                            'File "%s" does not exist.' % actualName)
            self.assertTrue(os.path.isfile(actualName),
                            'File "%s" is not an actual file.' % actualName)
    
    def testDocDirectoryPresence(self):
        self.assertTrue(self._importSuccess,
                        'Unsuccessful fisx.DataDir import')
        self.assertTrue(self.docDir is not None,
                        'Unassigned fisx.DataDir.DOC_DIR')
        #self.assertTrue(os.path.exists(self.docDir),
        #                'Directory "%s" does not exist' % self.docDir)
        #self.assertTrue(os.path.isdir(self.docDir),
        #                '"%s" expected to be a directory' % self.docDir)
        
def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(unittest.TestLoader().loadTestsFromTestCase(testData))
    else:
        # use a predefined order
        testSuite.addTest(testData("testDataDirectoryPresence"))
        testSuite.addTest(testData("testDataFilePresence"))
        #testSuite.addTest(testData("testDocDirectoryPresence"))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
