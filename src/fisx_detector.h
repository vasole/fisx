#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2018 European Synchrotron Radiation Facility
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
#ifndef FISX_DETECTOR_H
#define FISX_DETECTOR_H
#include "fisx_layer.h"

namespace fisx
{

/*!
  \class Detector
  \brief Class describing the detector.
*/

class Detector:public Layer
{
public:

    /*!
    layer like initialization.
    The detector is assumed to be cylindrical and the diameter is calculated.
    */
    Detector(const std::string & name="",  const double & density = 0.0,
                                           const double & thickness = 0.0,
                                           const double & funnyFactor = 1.0);

    void setMaterial(const std::string & materialName);

    void setMaterial(const Material & material);



    /*!
    Active area in cm2.
    The detector is assumed to be cylindrical and the diameter is calculated.
    */
    void setActiveArea(const double & area);

    /*!
    Diameter in cm2.
    For the time being the detector is assumed to be cylindrical.
    */
    void setDiameter(const double & diameter);

    /*!
    Returns the active area in cm2
    */
    double getActiveArea() const;

    /*!
    Returns the diameter in cm2
    */
    const double & getDiameter() const;

    /*!
    Sets the distance to reference layer in cm2
    */
    void setDistance(const double & distance);

    /*!
    Returns the distance to reference layer in cm2
    */
    const double & getDistance() const;

    /*!
    Make sure the internal escape peak cache is empty
    */
    void clearEscapePeakCache();

    /*!
    Getter to internal variable
    */
    const double getEscapePeakEnergyThreshold() const {return this->escapePeakEnergyThreshold;};
    const double getEscapePeakIntensityThreshold() const {return this->escapePeakIntensityThreshold;};
    const int getEscapePeakNThreshold() const {return this->escapePeakNThreshold;};
    const double getEscapePeakAlphaIn() const {return this->escapePeakAlphaIn;};

    /*!
    Returns escape peak energy and rate per detected photon of given energy.

    The optional arguments label and update serve for caching purposes.
    To make use of the cache, the label must not be empty (it should be used to identify the elementsLibrary)
    and update must be set to 0 because if it is set to 1 it will systematically empty the cascade cache.
    */
    std::map<std::string, std::map<std::string, double> > getEscape(const double & energy, \
                                                            const Elements & elementsLibrary, \
                                                            const std::string & label = "", \
                                                            const int & update = 1);
    void setMinimumEscapePeakEnergy(const double & energy);
    void setMinimumEscapePeakIntensity(const double & intensity);
    void setMaximumNumberOfEscapePeaks(const int & nPeaks);

private:
    double diameter ;
    double distance ;
    // Escape peak related parameters
    double escapePeakEnergyThreshold;
    double escapePeakIntensityThreshold;
    int escapePeakNThreshold;
    double escapePeakAlphaIn;
    std::map< std::string, std::map< double, std::map<std::string, std::map<std::string, double> > > >escapePeakCache;
    // TODO: Calibration, fano, noise, and so on.
};

} // namespace fisx

#endif //FISX_DETECTOR_H
