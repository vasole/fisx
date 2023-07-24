#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2023 European Synchrotron Radiation Facility
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
#ifndef FISX_XRFFUNCTION_H
#define FISX_XRFFUNCTION_H
#include "fisx_xrf.h"
#include <map>
#include <vector>
#include <iostream>

namespace fisx
{

class XRFFunction
{

typedef std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
        expectedLayerEmissionType;
typedef std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > >::const_iterator \
        iteratorExpectedLayerEmissionType;

public:
    /*!
    Default constructor
    */
    XRFFunction();

    /*!
    Constructor with configuration file
    */
    XRFFunction(const std::string & configurationFile);

    /*!
    Constructor with configuration instance
    */
    XRFFunction(const XRFConfig & configurationInstance);

    /*!
    Constructor with configuration file
    */
    XRFFunction(const XRF xrfInstance);

    /*!
    Evaluate function given parameters and configuration
    */
    std::vector<double> calculate(std::vector< double> x, std::map<std::string, double> parameters, XRFConfig configurationInstance);

    /*!
    Evaluate function given parameters using default configuration
    */
    std::vector<double> calculate(std::vector< double> x, std::map<std::string, double> parameters);

private:
    /*!
    The internal XRF instance
    */
    XRF xrf;

    /*!
    Some optimization
    */
    expectedLayerEmissionType lastMultilayerFluorescence;
};

} // namespace fisx

#endif // FISX_XRFFUNCTION_H
