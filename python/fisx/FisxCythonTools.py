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
import sys

if sys.version < "3.0":
    # Cython handles str and bytes properly under Python 2.x
    def toBytes(inputArgument, encoding="utf-8"):
        return inputArgument

    def toString(inputArgument, encoding="utf-8"):
        return inputArgument

    def toBytesKeys(inputDict, encoding="utf-8"):
        return inputDict

    def toBytesKeysAndValues(inputDict, encoding="utf-8"):
        return inputDict

    def toStringKeys(inputDict, encoding="utf-8"):
        return inputDict

    def toStringKeysAndValues(inputDict, encoding="utf-8"):
        return inputDict

    def toStringList(inputList, encoding="utf-8"):
        return inputList

else:
    def toBytes(inputArgument, encoding="utf-8"):
        if hasattr(inputArgument, "encode"):
            return inputArgument.encode(encoding)
        else:
            # I do not check for being already a bytes instance
            return inputArgument

    def toString(inputArgument, encoding="utf-8"):
        if hasattr(inputArgument, "decode"):
            return inputArgument.decode(encoding)
        else:
            # I do not check for being already a string instance
            return inputArgument

    def toBytesKeys(inputDict, encoding="utf-8"):
        return dict((key.encode(encoding), value) if hasattr(key, "encode") \
                    else (key, value) for key, value in inputDict.items())

    def toBytesKeysAndValues(inputDict, encoding="utf-8"):
        if not isinstance(inputDict, dict):
            return inputDict
        return dict((key.encode(encoding), toByteKeysAndValues(value)) if hasattr(key, "encode") \
                    else (key, toByteKeysAndValues(value)) for key, value in inputDict.items())

    def toStringKeysAndValues(inputDict, encoding="utf-8"):
        if not isinstance(inputDict, dict):
            return inputDict
        return dict((key.decode(encoding), toStringKeysAndValues(value)) if hasattr(key, "decode") \
                    else (key, toStringKeysAndValues(value)) for key, value in inputDict.items())

    def toStringKeys(inputDict, encoding="utf-8"):
        return dict((key.decode(encoding), value) if hasattr(key, "decode") \
                    else (key, value) for key, value in inputDict.items())

    def toStringList(inputList, encoding="utf-8"):
        return list([key.decode(encoding) if hasattr(key, "decode") \
                      else key for key in inputList])
