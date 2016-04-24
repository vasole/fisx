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
if 0:
    import os
    if os.path.exists('PyMca'):
        if os.path.exists('setup.py'):
            if os.path.exists('py2app_setup.py'):
                txt ='Tests cannnot be imported from top source directory'
                raise ImportError(txt)
    from PyMca.tests.TestAll import main as testAll
    from PyMca.tests.ConfigDictTest import test as testConfigDict
    from PyMca.tests.EdfFileTest import test as testEdfFile
    from PyMca.tests.ElementsTest import test as testElements
    from PyMca.tests.GefitTest import test as testGefit
    from PyMca.tests.PCAToolsTest import test as testPCATools
    from PyMca.tests.SpecfileTest import test as testSpecfile
    from PyMca.tests.specfilewrapperTest import test as testSpecfilewrapper

from fisx.tests.testAll import main as testAll
from fisx.tests.testXRF import test as testXRF
from fisx.tests.testSimpleSpecfile import test as testSimpleSpecfile
from fisx.tests.testEPDL97 import test as testEPDL97
from fisx.tests.testDataDir import test as testDataDir
from fisx.tests.testElements import test as testElements

