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
import unittest
import os
import sys
import glob
import unittest

def getSuite(auto=True):
    pythonFiles = glob.glob(os.path.join(os.path.dirname(__file__), "*.py"))
    sys.path.insert(0, os.path.dirname(__file__))
    testSuite = unittest.TestSuite()
    for fname in pythonFiles:
        if os.path.basename(fname) in ["__init__.py", "testAll.py"]:
            continue
        modName = os.path.splitext(os.path.basename(fname))[0]
        try:
            module = __import__(modName)
        except ImportError:
            print("Failed to import %s" % fname)
            continue
        if hasattr(module, "getSuite"):
            testSuite.addTest(module.getSuite(auto))
    return testSuite

def main(auto=True):
    return unittest.TextTestRunner(verbosity=2).run(getSuite(auto=auto))

if __name__ == '__main__':
    if len(sys.argv) > 1:
        auto = False
    else:
        auto = True
    if main(auto).wasSuccessful():
        print("Test suite was successful")
        sys.exit(0)
    else:
        print("Test suite failed")
        sys.exit(1)
