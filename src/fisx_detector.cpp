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
#include "fisx_detector.h"
#include <math.h>
#include <stdexcept>

namespace fisx
{

Detector::Detector(const std::string & name, const double & density, const double & thickness, \
                  const double & funnyFactor): Layer(name, density, thickness, funnyFactor)
{
    this->diameter =  0.0;
    this->distance = 10.0;
    this->escapePeakEnergyThreshold = 0.010,
    this->escapePeakIntensityThreshold = 1.0e-7;
    this->escapePeakNThreshold = 4;
    this->escapePeakAlphaIn = 90.;
    this->escapePeakCache.clear();
}

void Detector::setMaterial(const std::string & materialName)
{
    this->escapePeakCache.clear();
    this->Layer::setMaterial(materialName);
}

void Detector::setMaterial(const Material & material)
{
    this->escapePeakCache.clear();
    this->Layer::setMaterial(material);
}

void Detector::setMinimumEscapePeakEnergy(const double & energy)
{
    this->escapePeakEnergyThreshold = energy;
    this->escapePeakCache.clear();
}

void Detector::setMinimumEscapePeakIntensity(const double & intensity)
{
    this->escapePeakIntensityThreshold = intensity;
    this->escapePeakCache.clear();
}

void Detector::setMaximumNumberOfEscapePeaks(const int & nPeaks)
{
    this->escapePeakNThreshold = nPeaks;
    this->escapePeakCache.clear();
}

void Detector::clearEscapePeakCache()
{
    this->escapePeakCache.clear();
}

double Detector::getActiveArea() const
{
    double pi;
    pi = acos(-1.0);
    return (0.25 * pi) * (this->diameter * this->diameter);
}

void Detector::setActiveArea(const double & area)
{
    double pi;
    pi = acos(-1.0);
    if (area < 0)
    {
        throw std::invalid_argument("Negative detector area");
    }
    this->diameter = 2.0 * sqrt(area/pi);
}

void Detector::setDiameter(const double & diameter)
{
    if (diameter < 0)
    {
        throw std::invalid_argument("Negative detector diameter");
    }
    this->diameter = diameter;
}

void Detector::setDistance(const double & distance)
{
    if (distance <= 0)
    {
        throw std::invalid_argument("Negative detector distance");
    }
    this->distance = distance;
}

const double & Detector::getDiameter() const
{
    return this->diameter;
}

const double & Detector::getDistance() const
{
    return this->distance;
}

std::map<std::string, std::map<std::string, double> > Detector::getEscape(const double & energy,
                                                            const Elements & elementsLibrary,
                                                            const std::string & label,
                                                            const int & update)
{
    // To implement a cache at the detector level does not work because the same
    // detector can use different elementsLibrary instances and to use different
    // labels to identify them does not seem a good idea.
    // It is preferable to fill the cache at the elements library level
    if (label.size())
    {
        if (update != 0)
            this->escapePeakCache.clear();
        if (this->escapePeakCache.find(label) != this->escapePeakCache.end())
        {
            if (this->escapePeakCache[label].find(energy) != this->escapePeakCache[label].end())
            {
                return this->escapePeakCache[label][energy];
            }
        }
        this->escapePeakCache[label][energy] = elementsLibrary.getEscape(this->getComposition(elementsLibrary), \
                                                                 energy, \
                                                                 this->escapePeakEnergyThreshold, \
                                                                 this->escapePeakIntensityThreshold, \
                                                                 this->escapePeakNThreshold, \
                                                                 this->escapePeakAlphaIn,
                                                                 0);
        return this->escapePeakCache[label][energy];
    }
    else
    {
        return elementsLibrary.getEscape(this->getComposition(elementsLibrary), \
                                         energy, \
                                         this->escapePeakEnergyThreshold, \
                                         this->escapePeakIntensityThreshold, \
                                         this->escapePeakNThreshold, \
                                         this->escapePeakAlphaIn,
                                         0);
    }
}

} // namespace fisx
