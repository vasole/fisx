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
#include "fisx_xrf.h"
#include "fisx_math.h"
#include "fisx_simpleini.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace fisx
{

XRF::XRF()
{
    // initialize geometry with default parameters
    this->configuration = XRFConfig();
    this->setGeometry(45., 45.);
    //this->elements = NULL;
};

XRF::XRF(const std::string & fileName)
{
    this->readConfigurationFromFile(fileName);
    //this->elements = NULL;
}

void XRF::readConfigurationFromFile(const std::string & fileName)
{
    this->recentBeam = true;
    this->configuration.readConfigurationFromFile(fileName);
}

void XRF::setGeometry(const double & alphaIn, const double & alphaOut, const double & scatteringAngle)
{
    this->recentBeam = true;
    if (scatteringAngle < 0.0)
    {
        this->configuration.setGeometry(alphaIn, alphaOut, alphaIn + alphaOut);
    }
    else
    {
        this->configuration.setGeometry(alphaIn, alphaOut, scatteringAngle);
    }
}

void XRF::setBeam(const Beam & beam)
{
    this->recentBeam = true;
    this->configuration.setBeam(beam);
}

void XRF::setBeam(const double & energy, const double & divergency)
{
    this->recentBeam = true;
    this->configuration.setBeam(energy, divergency);
}

void XRF::setBeam(const std::vector<double> & energies, \
                 const std::vector<double> & weight, \
                 const std::vector<int> & characteristic, \
                 const std::vector<double> & divergency)
{
    this->configuration.setBeam(energies, weight, characteristic, divergency);
}

void XRF::setBeamFilters(const std::vector<Layer> &layers)
{
    this->recentBeam = true;
    this->configuration.setBeamFilters(layers);
}


void XRF::setSample(const std::vector<Layer> & layers, const int & referenceLayer)
{
    this->configuration.setSample(layers, referenceLayer);
}

void XRF::setSample(const std::string & name, \
                   const double & density, \
                   const double & thickness)
{
    std::vector<Layer> vLayer;
    vLayer.push_back(Layer(name, density, thickness, 1.0));
    this->configuration.setSample(vLayer, 0);
}

void XRF::setSample(const Layer & layer)
{
    std::vector<Layer> vLayer;
    vLayer.push_back(layer);
    this->configuration.setSample(vLayer, 0);
}

void XRF::setAttenuators(const std::vector<Layer> & attenuators)
{
    this->configuration.setAttenuators(attenuators);
}

void XRF::setDetector(const Detector & detector)
{
    this->configuration.setDetector(detector);
}

double XRF::getGeometricEfficiency(const int & sampleLayerIndex) const
{
    const Detector & detector = this->configuration.getDetector();
    const double PI = acos(-1.0);
    const double & sinAlphaOut = sin(this->configuration.getAlphaOut()*(PI/180.));
    const double & detectorDistance = detector.getDistance();
    const double & detectorDiameter = detector.getDiameter();
    double distance;
    const std::vector<Layer> & sample = this->configuration.getSample();
    std::vector<Layer>::size_type iLayer;
    const int & referenceLayerIndex = this->configuration.getReferenceLayer();
    const Layer* layerPtr;

    // if the detector diameter is zero, return 1
    if (detectorDiameter == 0.0)
    {
        return 1.0;
    }
    distance = detectorDistance;
    if ((distance == 0.0) && (sampleLayerIndex == 0))
    {
        return 0.5;
    }
    if (sampleLayerIndex < 0)
    {
        std::cout << "Negative sample layer index in getGeometricEfficiency " << sampleLayerIndex << std::endl;
        throw std::invalid_argument("Negative sample layer index in getGeometricEfficiency");
    }
    if (sampleLayerIndex != referenceLayerIndex)
    {
        if (sampleLayerIndex > referenceLayerIndex)
        {
            for (iLayer = (std::vector<Layer>::size_type) referenceLayerIndex; iLayer < (std::vector<Layer>::size_type) sampleLayerIndex; iLayer++)
            {
                layerPtr = &sample[iLayer];
                distance += (*layerPtr).getThickness() / sinAlphaOut;
            }
        }
        else
        {
            for (iLayer = (std::vector<Layer>::size_type) sampleLayerIndex; iLayer < (std::vector<Layer>::size_type) referenceLayerIndex; iLayer++)
            {
                layerPtr = &sample[iLayer];
                distance -= (*layerPtr).getThickness() / sinAlphaOut;
            }
        }
    }

    // we can calculate the geometric efficiency for the given layer
    // calculate geometric efficiency 0.5 * (1 - cos theta)
    return (0.5 * (1.0 - (distance / sqrt(pow(distance, 2) + pow(0.5 * detectorDiameter, 2)))));
}

std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                XRF::getMultilayerFluorescence( \
                const std::string & elementName, \
                const Elements & elementsLibrary, const int & sampleLayerIndex, \
                const std::string & lineFamily, const int & secondary, \
                const int & useGeometricEfficiency, const int & useMassFractions, \
                const double & optimizationFactor )
{
    std::vector<std::string> elementList;
    std::vector<std::string> familyList;
    std::vector<int> layerList;
    std::string tmpString;
    std::vector<std::string> tmpStringVector;

    elementList.push_back(elementName);

    if (lineFamily == "")
    {
        throw std::invalid_argument("Please specify K, L or M as peak family");
    }
    familyList.push_back(lineFamily);
    if (sampleLayerIndex < 0)
    {
        throw std::invalid_argument("Layer index cannot be negative");
    }
    layerList.push_back(sampleLayerIndex);
    return this->getMultilayerFluorescence(elementList, elementsLibrary, layerList, familyList, \
                                           secondary, useGeometricEfficiency, useMassFractions, \
                                           optimizationFactor);
}

double XRF::getEnergyThreshold(const std::string & elementName, const std::string & family, \
                                const Elements & elementsLibrary) const
{
    std::map<std::string, double> binding;
    binding = elementsLibrary.getBindingEnergies(elementName);
    if ((family == "K") || (family.size() == 2))
        return binding[family];

    if (family == "L")
    {
        if (binding["L3"] > 0)
            return binding["L3"];
        if (binding["L2"] > 0)
            return binding["L2"];
        return binding["L1"]; // It can be 0.0
    }

    if (family == "M")
    {
        if (binding["M5"] > 0)
            return binding["M5"];
        if (binding["M4"] > 0)
            return binding["M4"];
        if (binding["M3"] > 0)
            return binding["M3"];
        if (binding["M2"] > 0)
            return binding["M2"];
        return binding["M1"]; // It can be 0.0
    }
    return 0.0;
}

const XRFConfig & XRF::getConfiguration() const
{
    return this->configuration;

}

std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                XRF::getMultilayerFluorescence(const std::vector<std::string> & elementFamilyLayer, \
                const Elements & elementsLibrary, const int & secondary, \
                const int & useGeometricEfficiency, const int & useMassFractions, \
                const double & secondaryCalculationLimit)
{
    std::vector<std::string> elementList;
    std::vector<std::string> familyList;
    std::vector<int> layerList;
    std::vector<std::string>::size_type i;
    int layerIndex;
    std::string tmpString;
    std::vector<std::string> tmpStringVector;

    elementList.resize(elementFamilyLayer.size());
    familyList.resize(elementFamilyLayer.size());
    layerList.resize(elementFamilyLayer.size());

    for(i = 0; i < elementFamilyLayer.size(); i++)
    {
        tmpString = "";
        SimpleIni::parseStringAsMultipleValues(elementFamilyLayer[i], tmpStringVector, tmpString, ' ');
        // We should have a key of the form "Cr", "Cr K", or "Cr K 0"
        if(tmpStringVector.size() == 3)
        {
            elementList[i] = tmpStringVector[0];
            familyList[i] = tmpStringVector[1];
            if (!SimpleIni::stringConverter(tmpStringVector[2], layerIndex))
            {
                tmpString = "Unsuccessul conversion to layer integer: " + tmpStringVector[2];
                std::cout << tmpString << std::endl;
                throw std::invalid_argument(tmpString);
            }
            layerList[i] = layerIndex;
        }
        if(tmpStringVector.size() == 2)
        {
            elementList[i] = tmpStringVector[0];
            familyList[i] = tmpStringVector[1];
            layerList[i] = -1;
        }
        if(tmpStringVector.size() == 1)
        {
            elementList[i] = tmpStringVector[0];
            familyList[i] = "";
            layerList[i] = -1;
        }
    }
    return this->getMultilayerFluorescence(elementList, elementsLibrary, \
                                           layerList, familyList, secondary, useGeometricEfficiency, \
                                           useMassFractions, secondaryCalculationLimit);
}

} // namespace fisx
