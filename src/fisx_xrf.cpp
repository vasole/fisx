#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2023 European Synchrotron Radiation Facility
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
#include <algorithm>

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

void XRF::setSingleEnergyBeam(const double & energy, const double & divergency)
{
    this->recentBeam = true;
    this->configuration.setSingleEnergyBeam(energy, divergency);
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

void XRF::setUserBeamFilters(const std::vector<TransmissionTable> & userFilters)
{
    this->recentBeam = true;
    this->configuration.setUserBeamFilters(userFilters);
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

void XRF::setUserAttenuators(const std::vector<TransmissionTable> & userAttenuators)
{
    this->configuration.setUserAttenuators(userAttenuators);
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
                const double & optimizationFactor, \
                const Beam & overwritingBeam )
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
                                           optimizationFactor, overwritingBeam);
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

std::map<std::string, double> XRF::getLayerComposition(const Layer & layer, const Elements & elements) const
{
    std::map <std::string, double> composition;
    Material layerMaterial;

    if (layer.hasMaterialComposition())
    {
        // It has an actual material with the composition. However, the name of the material can be unset
        // and one cannot simply use its name via the line commented below.
        // return elements.getComposition(layer.getMaterialName(), this->getMaterials());
        return elements.getComposition(layer.getMaterial().getComposition(), this->getMaterials());
    }
    else
    {
        // It has a material name, element or formula as composition
        return elements.getComposition(layer.getMaterialName(), this->getMaterials());
    }
}

std::map<std::string, double> XRF::getLayerMassAttenuationCoefficients( \
                                                const Layer & layer,
                                                const double & energy,
                                                const Elements & elements,
                                                const std::map<std::string, double> & layerComposition) const
{
    if (layerComposition.size() > 0)
        return elements.getMassAttenuationCoefficients(layerComposition, energy, 1);
    else
        return elements.getMassAttenuationCoefficients(this->getLayerComposition(layer, elements), energy);
}

std::map<std::string, std::vector<double> > XRF::getLayerMassAttenuationCoefficients(const Layer & layer,
                                                const std::vector<double> & energy,
                                                const Elements & elements,
                                                const std::map<std::string, double> & layerComposition) const
{
    if (layerComposition.size() > 0)
    {
        // this is an optimization path where the layer composition is already in terms of elements
        // and not in terms of materials or formulas to be parsed
        return elements.getMassAttenuationCoefficients(layerComposition, energy, 1);
    }
    else
    {
        return elements.getMassAttenuationCoefficients(this->getLayerComposition(layer, elements), energy);
    }
}

double XRF::getLayerTransmission(const Layer & layer,
                                 const double & energy,
                                 const Elements & elements,
                                 const double & angle,
                                 const std::map<std::string, double> & layerComposition) const
{
    std::vector<double> energies;

    energies.resize(0);
    energies.push_back(energy);
    return this->getLayerTransmission(layer, energies, elements, angle, layerComposition)[0];
}

std::vector<double> XRF::getLayerTransmission(const Layer & layer,
                                             const std::vector<double> & energy,
                                             const Elements & elements,
                                             const double & angle,
                                             const std::map<std::string, double> & layerComposition) const
{
    const double PI = std::acos(-1.0);
    std::vector<double>::size_type i;
    std::vector<double> tmpDoubleVector;
    double tmpDouble;

    if (angle == 90.0)
    {
        tmpDouble = layer.getDensity() * layer.getThickness();
    }
    else
    {
        if (angle < 0)
            tmpDouble = std::sin(-angle * PI / 180.);
        else
            tmpDouble = std::sin(angle * PI / 180.);
        tmpDouble = layer.getDensity() * layer.getThickness() / tmpDouble;
    }

    if(tmpDouble <= 0.0)
    {
        std::string msg;
        msg = "Layer " + layer.getName() + " thickness is " + elements.toString(tmpDouble) + " g/cm2";
        throw std::runtime_error( msg );
    }

    tmpDoubleVector = this->getLayerMassAttenuationCoefficients(layer, energy, elements, layerComposition)["total"];

    for (i = 0; i < tmpDoubleVector.size(); i++)
    {
        tmpDoubleVector[i] = (1.0 - layer.getFunnyFactor()) + \
                              (layer.getFunnyFactor() * exp(-(tmpDouble * tmpDoubleVector[i])));
    }
    return tmpDoubleVector;
}

std::vector<std::pair<std::string, double> > XRF::getLayerPeakFamilies( \
                                                        const Layer & layer,
                                                        const double & energy,
                                                        const Elements & elements,
                                                        const std::map<std::string, double> & layerComposition) const
{
    std::map<std::string, double> composition;
    std::vector<std::string> elementsList;
    std::map<std::string, double>::const_iterator c_it;

    if (layerComposition.size() > 0)
        composition = layerComposition;
    else
        composition = this->getLayerComposition(layer, elements);

    for(c_it = composition.begin(); c_it != composition.end(); ++c_it)
    {
        if (std::find(elementsList.begin(), elementsList.end(), c_it->first) == elementsList.end())
        {
            elementsList.push_back(c_it->first);
        }
    }

    return elements.getPeakFamilies(elementsList, energy);

}


std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                XRF::getMultilayerFluorescence(const std::vector<std::string> & elementFamilyLayer, \
                const Elements & elementsLibrary, const int & secondary, \
                const int & useGeometricEfficiency, const int & useMassFractions, \
                const double & secondaryCalculationLimit, const Beam & overwritingBeam)
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
                                           useMassFractions, secondaryCalculationLimit, overwritingBeam);
}

void XRF::printConfiguration() const
{
    std::cout << this->getConfiguration() << std::endl;
}

} // namespace fisx
