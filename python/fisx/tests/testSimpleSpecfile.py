__author__ = "V.A. Sole - ESRF Data Analysis"
import unittest
import sys
import os
import gc
import tempfile

class testSimpleSpecfile(unittest.TestCase):
    def setUp(self):
        """
        import the module
        """
        try:
            from fisx import SimpleSpecfile
            self.specfileClass = SimpleSpecfile
        except:
            self.specfileClass = None
        if self.specfileClass is not None:
            text  = "#F \n"
            text += "\n"
            text += "#S 10  Undefined command 0\n"
            text += "#N 3\n"
            text += "#L First label  Second label  Third label\n"
            text += "10  100  1000\n"
            text += "20  400  8000\n"
            text += "30  900  270000\n"
            text += "\n"
            text += "#S 20  Undefined command 1\n"
            text += "#N 3\n"
            text += "#L First  Second  Third\n"
            text += "1.3  1  1\n"
            text += "2.5  4  8\n"
            text += "3.7  9  27\n"
            text += "\n"
            tmpFile = tempfile.mkstemp(text=False)
            if sys.version < '3.0':
                os.write(tmpFile[0], text)
            else:
                os.write(tmpFile[0], bytes(text, 'utf-8'))
            os.close(tmpFile[0])
            self.fname = tmpFile[1]

    def tearDown(self):
        """clean up any possible files"""
        # make sure the file handle is free
        self._sf = None
        self._scan = None
        # this should free the handle
        gc.collect()
        if self.specfileClass is not None:
            if os.path.exists(self.fname):
                os.remove(self.fname)

    def testSimpleSpecfileImport(self):
        #"""Test successful import"""
        self.assertTrue(self.specfileClass is not None,
                        'Unsuccessful fisx.SimpleSpecfile import')

    def testSimpleSpecfileReading(self):
        #"""Test specfile readout"""
        self.testSimpleSpecfileImport()
        self._sf = self.specfileClass(self.fname)

        # test the number of found scans
        nScans = self._sf.getNumberOfScans()
        self.assertEqual(nScans, 2,
                         'Expected to read 2 scans, read %s' %\
                         nScans)

        # test scan iteration selection method
        labels = self._sf.getScanLabels(1)
        expectedLabels = ['First', 'Second', 'Third']
        self.assertEqual(len(labels), 3,
                         'Expected to read 3 scans, got %s' % len(labels))
        for i in range(3):
            self.assertEqual(labels[i], expectedLabels[i],
                    'Read "%s" instead of "%s"' %\
                     (labels[i], expectedLabels[i]))

        gc.collect()

    def testSimpleSpecfileReadingCompatibleWithUserLocale(self):
        #"""Test specfile compatible with C locale"""
        self.testSimpleSpecfileImport()
        self._sf = self.specfileClass(self.fname)
        data = self._sf.getScanData(1)
        self._sf = None
        self.assertEqual(data[0] [0], 1.3,
                    'Read %f instead of %f' %\
                    (data[0][0], 1.3))
        self.assertEqual(data[1] [0], 2.5,
                    'Read %f instead of %f' %\
                    (data[1] [0], 2.5))
        self.assertEqual(data[2] [0], 3.7,
                    'Read %f instead of %f' %\
                    (data[2] [0], 3.7))
        gc.collect()

    def testSimpleSpecfileVersusPyMca(self):
        import glob
        try:
            from PyMca5 import PyMcaDataDir
            from PyMca5.PyMcaIO import specfilewrapper as Specfile
        except ImportError:
            print("\n ****** Skipping: PyMca not installed ******* \n")
            return
        fileList = glob.glob(os.path.join(PyMcaDataDir.PYMCA_DATA_DIR, '*.dat'))
        self.testSimpleSpecfileImport()
        for filename in fileList:
            self._sf = self.specfileClass(filename)
            a = self._sf
            nScans = a.getNumberOfScans()
            sf = Specfile.Specfile(filename)
            for i in range(nScans):
                scan = sf[i]
                labels0 = a.getScanLabels(i)
                labels1 = scan.alllabels()
                for j in range(len(labels0)):
                    self.assertEqual(labels0[j], labels1[j],
                                "Read <%s> instead of <%s>" %\
                                     (labels0[j], labels1[j]))
                data0 = a.getScanData(i)
                data1 = scan.data().T
                M = abs(data0 - data1).max()
                if M > 1.0e-8:
                    raise ValueError("Error reading data")
        

def getSuite(auto=True):
    testSuite = unittest.TestSuite()
    if auto:
        testSuite.addTest(\
            unittest.TestLoader().loadTestsFromTestCase(testSimpleSpecfile))
    else:
        # use a predefined order
        testSuite.addTest(testSimpleSpecfile("testSimpleSpecfileImport"))
        testSuite.addTest(testSimpleSpecfile("testSimpleSpecfileReading"))
        testSuite.addTest(\
            testSimpleSpecfile("testSimpleSpecfileReadingCompatibleWithUserLocale"))
        testSuite.addTest(\
            testSimpleSpecfile("testSimpleSpecfileVersusPyMca"))
    return testSuite

def test(auto=False):
    unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    test()
