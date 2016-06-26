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
import os
# this will be filled by the setup
FISX_DATA_DIR = 'DATA_DIR_FROM_SETUP'
# this is to be filled by the setup
FISX_DOC_DIR = 'DOC_DIR_FROM_SETUP'

# this is used in build directory
if not os.path.exists(FISX_DATA_DIR):
    tmp_dir = os.path.dirname(os.path.abspath(__file__))
    old_tmp_dir = tmp_dir + "dummy"
    basename = "fisx_data"
    FISX_DATA_DIR = os.path.join(tmp_dir, "fisx", basename)
    while (len(FISX_DATA_DIR) > 14) and (tmp_dir != old_tmp_dir):
        if os.path.exists(FISX_DATA_DIR):
            break
        old_tmp_dir = tmp_dir
        tmp_dir = os.path.dirname(tmp_dir)
        FISX_DATA_DIR = os.path.join(tmp_dir, "fisx", basename)

if not os.path.exists(FISX_DATA_DIR):
    FISX_DATA_DIR = os.getenv("FISX_DATA_DIR")
    if FISX_DATA_DIR is not None:
        if not os.path.exists(FISX_DATA_DIR):
            raise IOError('%s directory set from environent not found' % FISX_DATA_DIR)
        else:
            txt = "WARNING: Taking FISX_DATA_DIR from environement.\n"
            txt += "Use it at your own risk."
            print(txt)
    else:
        raise IOError('%s directory not found' % basename)

if not os.path.exists(FISX_DOC_DIR):
    tmp_dir = os.path.dirname(os.path.abspath(__file__))
    old_tmp_dir = tmp_dir + "dummy"
    basename = "fisx_data"
    FISX_DOC_DIR = os.path.join(tmp_dir,basename)
    while (len(FISX_DOC_DIR) > 14) and (tmp_dir != old_tmp_dir):
        if os.path.exists(FISX_DOC_DIR):
            break
        old_tmp_dir = tmp_dir
        tmp_dir = os.path.dirname(tmp_dir)
        FISX_DOC_DIR = os.path.join(tmp_dir, "fisx", basename)

if not os.path.exists(FISX_DOC_DIR):
    FISX_DOC_DIR = os.getenv("FISX_DOC_DIR")
    if FISX_DOC_DIR is not None:
        if not os.path.exists(FISX_DOC_DIR):
            raise IOError('%s directory not found' % basename)
    else:
        # use the data dir as doc dir
        FISX_DOC_DIR = FISX_DATA_DIR
