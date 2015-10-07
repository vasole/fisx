#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014 V. Armando Sole
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

std::map<std::string, std::map<std::string, double> > XRF::getFluorescence(const std::string & elementName, \
                const Elements & elementsLibrary, const int & sampleLayerIndex, \
                const std::string & lineFamily, const int & secondary, const int & useGeometricEfficiency)
{
    // get all the needed configuration
    const Beam & beam = this->configuration.getBeam();
    std::vector<std::vector<double> >actualRays = beam.getBeamAsDoubleVectors();
    std::vector<double>::size_type iRay;
    const std::vector<Layer> & filters = this->configuration.getBeamFilters();;
    const std::vector<Layer> & sample = this->configuration.getSample();
    const std::vector<Layer> & attenuators = this->configuration.getAttenuators();
    const Layer* layerPtr;
    std::vector<Layer>::size_type iLayer;
    const Detector & detector = this->configuration.getDetector();
    const Element & element = elementsLibrary.getElement(elementName);
    std::string msg;
    std::map<std::string, std::map<std::string, double> > result;
    std::map<std::string, std::map<std::string, double> > actualResult;
    const double PI = acos(-1.0);
    const double & alphaIn = this->configuration.getAlphaIn();
    const double & alphaOut = this->configuration.getAlphaOut();
    const double & detectorDiameter = detector.getDiameter();
    double geometricEfficiency;
    double sinAlphaIn = sin(alphaIn*(PI/180.));
    double sinAlphaOut = sin(alphaOut*(PI/180.));
    double tmpDouble;
    std::string tmpString;

    std::cout << "WARNING: This method is obsolete. Use any of the getMultilayerFluorescnce ones" << std::endl;
    if (actualRays.size() == 0)
    {
        // no excitation beam
        if (sample.size() > 0)
        {
            msg = "Sample is defined but beam it is not!";
            throw std::invalid_argument( msg );
        }
        if (lineFamily.size() == 0)
        {
            msg = "No sample and no beam. Please specify family of lines to get theoretical ratios.";
            throw std::invalid_argument( msg );
        }
        // we just have to get the theoretical ratios and only deal with attenuators and detector
        if (lineFamily == "K")
        {
            // get the K lines
            result[lineFamily] =  element.getXRayLines(lineFamily);
        }
        else if ( (lineFamily == "L1") || (lineFamily == "L2") || (lineFamily == "L3"))
        {
            // get the relevant L subshell lines
            result[lineFamily] =  element.getXRayLines(lineFamily);
        }
        else if ( (lineFamily == "M1") || (lineFamily == "M2") || (lineFamily == "M3") || (lineFamily == "M4") || (lineFamily == "M5"))
        {
            result[lineFamily] =  element.getXRayLines(lineFamily);
        }
        else if ((lineFamily == "L") || (lineFamily == "M"))
        {
            std::cout << "I should assume an initial vacancy distribution given by the jumps. " << std::endl;
            msg = "Excitation energy needed in order to properly calculate intensity ratios.";
            throw std::invalid_argument( msg );
        }
        else
        {
            msg = "Excitation energy needed in order to properly calculate intensity ratios.";
            throw std::invalid_argument( msg );
        }
    }
    std::vector<double> & energies = actualRays[0];
    std::vector<double> doubleVector;

    // beam is ordered
    //maxEnergy = energies[energies.size() - 1];

    // get the beam after the beam filters
    std::vector<double> muTotal;
    muTotal.resize(energies.size());
    std::fill(muTotal.begin(), muTotal.end(), 0.0);
    doubleVector.resize(energies.size());
    std::fill(doubleVector.begin(), doubleVector.end(), 1.0);
    for (iLayer = 0; iLayer < filters.size(); iLayer++)
    {
        layerPtr = &filters[iLayer];
        doubleVector = (*layerPtr).getTransmission(energies, elementsLibrary);
        for (iRay = 0; iRay < energies.size(); iRay++)
        {
            actualRays[1][iRay] *= doubleVector[iRay];
        }
    }

    // this has sense if we put all the previous stuff cached

    std::vector<double> weights;
    weights = actualRays[1];
    std::fill(muTotal.begin(), muTotal.end(), 0.0);
    for (iLayer = 0; iLayer < sampleLayerIndex; iLayer++)
    {
        layerPtr = &sample[iLayer];
        doubleVector = (*layerPtr).getTransmission(energies, \
                                            elementsLibrary, alphaIn);
        for (iRay = 0; iRay < energies.size(); iRay++)
        {
            weights[iRay] *= doubleVector[iRay];
        }
    }

    // we can already calculate the geometric efficiency
    if ((useGeometricEfficiency != 0) && (detectorDiameter > 0.0))
    {
        // calculate geometric efficiency 0.5 * (1 - cos theta)
        geometricEfficiency = this->getGeometricEfficiency(sampleLayerIndex);
    }
    else
    {
        geometricEfficiency = 1.0 ;
    }
    // we have reached the layer we are interesed on
    // calculate its total mass attenuation coefficient at each incident energy
    std::fill(muTotal.begin(), muTotal.end(), 0.0);
    std::map<std::string, double> sampleLayerComposition;
    layerPtr = &sample[sampleLayerIndex];
    if ((*layerPtr).hasMaterialComposition())
    {
        const Material & material = (*layerPtr).getMaterial();
        doubleVector = elementsLibrary.getMassAttenuationCoefficients( \
                            material.getComposition(), energies)["total"];
        if (secondary > 0)
        {
            sampleLayerComposition = material.getComposition();
        }
    }
    else
    {
       doubleVector = elementsLibrary.getMassAttenuationCoefficients(\
                            (*layerPtr).getMaterialName(), energies)["total"];
        if (secondary > 0)
        {
            sampleLayerComposition = elementsLibrary.getComposition((*layerPtr).getMaterialName());
        }
    }
    for (iRay = 0; iRay < energies.size(); iRay++)
    {
        muTotal[iRay] = doubleVector[iRay] /sinAlphaIn;
    }

    //std::cout << this->configuration << std::endl;

    std::map<std::string, std::map<std::string, double> > tmpResult;
    std::map<std::string, std::map<std::string, double> >::const_iterator c_it;
    std::map<std::string, double>::const_iterator mapIt;
    std::map<std::string, double>::const_iterator mapIt2;
    std::map<std::string, double> muTotalFluo;
    std::map<std::string, double> detectionEfficiency;
    std::vector<double> sampleLayerEnergies;
    std::vector<std::string> sampleLayerEnergyNames;
    std::vector<double> sampleLayerRates;
    std::vector<double> sampleLayerMuTotal;
    std::vector<double>::size_type iLambda;

    iRay = energies.size();
    while (iRay > 0)
    {
        --iRay;

        if (secondary > 0)
        {
            sampleLayerEnergies.clear();
            sampleLayerEnergyNames.clear();
            sampleLayerRates.clear();
            sampleLayerMuTotal.clear();
            layerPtr = &sample[sampleLayerIndex];
            for (mapIt = sampleLayerComposition.begin(); \
                 mapIt != sampleLayerComposition.end(); ++mapIt)
            {
                // get excitation factors for each element
                tmpResult = elementsLibrary.getExcitationFactors(mapIt->first,
                                                    energies[iRay], weights[iRay]);
                //and add the energies and rates to the sampleLayerLines
                for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
                {
                    mapIt2 = c_it->second.find("energy");
                    sampleLayerEnergies.push_back(mapIt2->second);
                    mapIt2 = c_it->second.find("rate");
                    sampleLayerRates.push_back(mapIt2->second * mapIt->second);
                    tmpString = mapIt->first + " " + c_it->first;
                    sampleLayerEnergyNames.push_back(tmpString);
                }
                sampleLayerMuTotal = (*layerPtr).getMassAttenuationCoefficients(sampleLayerEnergies, \
                                                                    elementsLibrary)["total"];
            }
        }
        // energy = energies[iRay];
        // we should check the energies that have to be considered
        // now for *each* line, we have to calculate how the "rate" key is to be modified
        tmpResult = elementsLibrary.getExcitationFactors(elementName, energies[iRay], weights[iRay]);
        if (muTotalFluo.size() == 0)
        {
            layerPtr = &sample[sampleLayerIndex];
            // we have to calculate the sample total mass attenuation coefficients at the fluorescent energies
            for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
            {
                mapIt = c_it->second.find("energy");
                muTotalFluo[c_it->first] = (*layerPtr).getMassAttenuationCoefficients( \
                                                                mapIt->second, \
                                                elementsLibrary)["total"] / sinAlphaOut;
            }
            // calculate the transmission of the fluorescence photon in the way back.
            // it will be the same for each incident energy.
            // in the sample upper layers
            // in the attenuators
            // the geometric factor
            // the detector efficiency
            for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
            {
                mapIt = c_it->second.find("energy");
                tmpDouble = mapIt->second;
                detectionEfficiency[c_it->first] = 1.0;

                // transmission through upper layers
                iLayer = sampleLayerIndex;
                while (iLayer > 0)
                {
                    --iLayer;
                    layerPtr = &sample[iLayer];
                    detectionEfficiency[c_it->first] *= (*layerPtr).getTransmission( mapIt->second, \
                                                            elementsLibrary, alphaOut);
                }
                // transmission through attenuators
                for (iLayer = 0; iLayer < attenuators.size(); iLayer++)
                {
                    layerPtr = &attenuators[iLayer];
                    detectionEfficiency[c_it->first] *= (*layerPtr).getTransmission( mapIt->second, \
                                                            elementsLibrary, 90.);
                }
                //std::cout << mapIt->second << " " << c_it->first << "" << detectionEfficiency[c_it->first] << std::endl;
                // detection efficienty decomposed in geometric and intrinsic
                if (detectorDiameter > 0.0)
                {
                    // apply geometric efficiency 0.5 * (1 - cos theta)
                    detectionEfficiency[c_it->first] *= geometricEfficiency;
                }
                if (detector.hasMaterialComposition() || (detector.getMaterialName().size() > 0))
                {
                    // calculate intrinsic efficiency
                    detectionEfficiency[c_it->first] *= (1.0 - detector.getTransmission( mapIt->second, \
                                                                            elementsLibrary, 90.0));
                }
            }

            actualResult = tmpResult;
            for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
            {
                actualResult[c_it->first]["rate"] = 0.0;
            }
        }
        for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
        {
            //The self attenuation term
            tmpDouble = (muTotal[iRay] + muTotalFluo[c_it->first]);
            //std::cout << "sum of mass att coef " << tmpDouble << std::endl;
            tmpDouble = (1.0 - exp(- tmpDouble * \
                                   sample[sampleLayerIndex].getDensity() * \
                                   sample[sampleLayerIndex].getThickness())) / (tmpDouble * sinAlphaIn);
            mapIt = c_it->second.find("rate");
            //std::cout << "RATE = " << mapIt->second << std::endl;
            //std::cout << "ATT TERM = " << tmpDouble << std::endl;
            //std::cout << "EFFICIENCY = " << detectionEfficiency[c_it->first] << std::endl;
            actualResult[c_it->first]["rate"] += mapIt->second * tmpDouble  * \
                                                detectionEfficiency[c_it->first];
        }

        // probably I sould calculate this first to prevent adding small numbers to a bigger one
        if (secondary > 0)
        {
            // std::cout << "sample energies = " << sampleLayerEnergies.size() << std::endl;
            for(iLambda = 0; iLambda < sampleLayerEnergies.size(); iLambda++)
            {
                // analogous to incident beam
                //if (sampleLayerEnergies[iLambda] < this->getEnergyThreshold(elementName, lineFamily, elementsLibrary))
                //    continue;
                tmpResult = elementsLibrary.getExcitationFactors(elementName, \
                            sampleLayerEnergies[iLambda], sampleLayerRates[iLambda]);
                for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
                {
                    tmpDouble = Math::deBoerL0(muTotal[iRay],
                                               muTotalFluo[c_it->first],
                                               sampleLayerMuTotal[iLambda],
                                               sample[sampleLayerIndex].getDensity(),
                                               sample[sampleLayerIndex].getThickness());
                    /*
                    std::cout << "energy0" << energies[iRay] << "L0" << tmpDouble << std::endl;
                    std::cout << "muTotal[iRay] " << muTotal[iRay] << std::endl;
                    std::cout << "muTotalFluo[c_it->first] " << muTotalFluo[c_it->first] << std::endl;
                    std::cout << "sampleLayerMuTotal[iLambda] " << sampleLayerMuTotal[iLambda] << std::endl;
                    */
                    tmpDouble += Math::deBoerL0(muTotalFluo[c_it->first],
                                                muTotal[iRay],
                                                sampleLayerMuTotal[iLambda],
                                                sample[sampleLayerIndex].getDensity(),
                                                sample[sampleLayerIndex].getThickness());
                    tmpDouble *= (0.5/sinAlphaIn);
                    mapIt = c_it->second.find("rate");
                    actualResult[c_it->first]["rate"] += mapIt->second * tmpDouble * \
                                                         detectionEfficiency[c_it->first];
                }
            }
        }
    }
    return actualResult;
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

    if (sampleLayerIndex != referenceLayerIndex)
    {
        if (sampleLayerIndex > referenceLayerIndex)
        {
            for (iLayer = referenceLayerIndex; iLayer < sampleLayerIndex; iLayer++)
            {
                layerPtr = &sample[iLayer];
                distance += (*layerPtr).getThickness() / sinAlphaOut;
            }
        }
        else
        {
            for (iLayer = sampleLayerIndex; iLayer < referenceLayerIndex; iLayer++)
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
    std::vector<std::string>::size_type i;
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

std::map<std::string, std::vector<double> > XRF::getSpectrum(const std::vector<double> & channel, \
                const std::map<std::string, double> & detectorParameters, \
                const std::map<std::string, double> & shapeParameters, \
                const std::map<std::string, double> & peakFamilyArea, \
                const expectedLayerEmissionType & emissionRatios) const
{
    std::map<std::string, double>::const_iterator c_it;
    std::map<std::string, std::vector<double> > result;
    for (c_it = shapeParameters.begin(); c_it != shapeParameters.end(); ++c_it)
    {
        std::cout << "Key = " << c_it->first << " Value " << c_it->second << std::endl;
    }
    return result;
}

void XRF::getSpectrum(double * channel, double * energy, double *spectrum, int nChannels, \
                const std::map<std::string, double> & detectorParameters, \
                const std::map<std::string, double> & shapeParameters, \
                const std::map<std::string, double> & peakFamilyArea, \
                const expectedLayerEmissionType & emissionRatios) const
{
    int i;
    std::string detectorKeys[5] = {"Zero", "Gain", "Noise", "Fano", "QuantumEnergy"};
    std::string shapeKeys[6] = {"ShortTailArea", "ShortTailSlope", "LongTailArea", "LongTailSlope", "StepHeight", \
                                "Eta"};
    std::map<std::string, double>::const_iterator c_it;

    std::string tmpString;
    std::vector<std::string> tmpStringVector;
    double zero, gain, noise, fano, quantum;
    int layerIndex;

    double area;
    double position;
    double fwhm;
    double shortTailArea = 0.0, shortTailSlope = -1.0;
    double longTailArea = 0.0, longTailSlope = -1.0;
    double stepHeight = 0.0;
    double eta = 0.0;

    for (c_it = detectorParameters.begin(); c_it != detectorParameters.end(); ++c_it)
    {
        tmpString = c_it->first;
        SimpleIni::toUpper(tmpString);
        if ( tmpString == "ZERO")
        {
            zero = c_it->second;
            continue;
        }
        if (tmpString == "GAIN")
        {
            gain = c_it->second;
            continue;
        }
        if (tmpString == "NOISE")
        {
            noise = c_it->second;
            continue;
        }
        if (tmpString == "FANO")
        {
            fano = c_it->second;
            continue;
        }
        if (tmpString == "QUANTUMENERGY")
        {
            quantum = c_it->second;
            continue;
        }
        std::cout << "WARNING: Unused detector parameter "<< c_it->first << " with value " << c_it->second << std::endl;
    }

    for (i = 0; i < nChannels; i++)
    {
        energy[i] = zero + gain * channel[i];
    }
    for (i = 0; i < nChannels; i++)
    {
        spectrum[i] = 0.0;
    }

    for (c_it = peakFamilyArea.begin(); c_it != peakFamilyArea.end(); ++c_it)
    {
        std::map<int, std::map<std::string, std::map<std::string, double> > > ::const_iterator layerIterator;
        std::map<std::string, std::map<std::string, double> >::const_iterator lineIterator;
        std::map<std::string, double>::const_iterator ratePointer;
        std::vector<double> layerTotalSignal;
        double totalSignal;
        iteratorExpectedLayerEmissionType emissionRatiosPointer;
        emissionRatiosPointer = emissionRatios.find(c_it->first);
        // check if the description of that peak multiplet is available
        if (emissionRatiosPointer != emissionRatios.end())
        {
            // In this case emission ratios has the form "Cr K".
            // Remmeber that peakFamily can have the form "Cr K 0"
            // We have to sum all the signals, to normalize to unit area, and multiply by the supplied
            // area. This could have been already done ...
            // loop for each layer
            layerTotalSignal.clear();
            totalSignal = 0.0;
            for (layerIterator = emissionRatiosPointer->second.begin();
                 layerIterator != emissionRatiosPointer->second.end(); ++layerIterator)
            {
                layerTotalSignal.push_back(0.0);
                for (lineIterator = layerIterator->second.begin(); \
                     lineIterator != layerIterator->second.end(); ++lineIterator)
                {
                    ratePointer = lineIterator->second.find("rate");
                    if (ratePointer == lineIterator->second.end())
                    {
                        tmpString = "Keyword <rate> not found!!!";
                        std::cout << tmpString << std::endl;
                        throw std::invalid_argument(tmpString);
                    }
                    layerTotalSignal[layerTotalSignal.size() - 1] += ratePointer->second;
                }
                totalSignal += layerTotalSignal[layerTotalSignal.size() - 1];
            }
           // Now we already have area (provided) and ratio (dividing by totalSignal).
           // We can therefore calculate the signal keeping the proper ratios.
            for (layerIterator = emissionRatiosPointer->second.begin();
                 layerIterator != emissionRatiosPointer->second.end(); ++layerIterator)
            {
                for (lineIterator = layerIterator->second.begin(); \
                     lineIterator != layerIterator->second.end(); ++lineIterator)
                {
                    ratePointer = lineIterator->second.find("rate");
                    area = c_it->second * (ratePointer->second / totalSignal);
                    ratePointer = lineIterator->second.find("energy");
                    if (ratePointer == lineIterator->second.end())
                    {
                        tmpString = "Keyword <energy> not found!!!";
                        std::cout << tmpString << std::endl;
                        throw std::invalid_argument(tmpString);
                    }
                    position = ratePointer->second;
                    fwhm = Math::getFWHM(position, noise, fano, quantum);
                    for (i = 0; i < nChannels; i++)
                    {
                        spectrum[i] += Math::hypermet(energy[i], \
                                                      area, position, fwhm, \
                                                      shortTailArea, shortTailSlope, \
                                                      longTailArea, longTailSlope, stepHeight);
                    }
                }
            }
        }
        else
        {
            tmpString = "";
            SimpleIni::parseStringAsMultipleValues(c_it->first, tmpStringVector, tmpString, ' ');
            if(tmpStringVector.size() != 3)
            {
                tmpString = "Unsuccessul conversion to Element, Family, layer index: " + c_it->first;
            }

            // We should have a key of the form "Cr K 0"
            if (!SimpleIni::stringConverter(tmpStringVector[2], layerIndex))
            {
                tmpString = "Unsuccessul conversion to layer integer: " + tmpStringVector[2];
                std::cout << tmpString << std::endl;
                throw std::invalid_argument(tmpString);
            }
            // TODO: Deal with Ka, Kb, L, L1, L2, L3, ...
            tmpString = tmpStringVector[0] + " " + tmpStringVector[1];
            emissionRatiosPointer = emissionRatios.find(tmpString);
            if (emissionRatiosPointer == emissionRatios.end())
            {
                tmpString = "Undefined emission ratios for element " + tmpStringVector[0] +\
                            " family " + tmpStringVector[1];
                std::cout << tmpString << std::endl;
                throw std::invalid_argument(tmpString);
            }
            // Emission ratios has the form "Cr K" but we have received peakFamily can have the form "Cr K index"
            // We have to to normalize the signal from that element, family and layer to unit area,
            // and multiply by the supplied area. This could have been already done ...
            layerIterator = emissionRatiosPointer->second.find(layerIndex);
            if (layerIterator == emissionRatiosPointer->second.end())
            {
                tmpString = "I do not have information for layer number " + tmpStringVector[2];
                std::cout << tmpString << std::endl;
                throw std::invalid_argument(tmpString);
            }
            layerTotalSignal.clear();
            totalSignal = 0.0;
            layerTotalSignal.push_back(0.0);
            for (lineIterator = layerIterator->second.begin(); \
                 lineIterator != layerIterator->second.end(); ++lineIterator)
            {
                ratePointer = lineIterator->second.find("rate");
                if (ratePointer == lineIterator->second.end())
                {
                    tmpString = "Keyword <rate> not found!!!";
                    std::cout << tmpString << std::endl;
                    throw std::invalid_argument(tmpString);
                }
                layerTotalSignal[layerTotalSignal.size() - 1] += ratePointer->second;
            }
            totalSignal += layerTotalSignal[layerTotalSignal.size() - 1];
            // Now we already have area (provided) and ratio (dividing by totalSignal).
            // We can therefore calculate the signal keeping the proper ratios.
            for (lineIterator = layerIterator->second.begin(); \
                 lineIterator != layerIterator->second.end(); ++lineIterator)
            {
                ratePointer = lineIterator->second.find("rate");
                area = c_it->second * (ratePointer->second / totalSignal);
                ratePointer = lineIterator->second.find("energy");
                if (ratePointer == lineIterator->second.end())
                {
                    tmpString = "Keyword <energy> not found!!!";
                    std::cout << tmpString << std::endl;
                    throw std::invalid_argument(tmpString);
                }
                position = ratePointer->second;
                fwhm = Math::getFWHM(position, noise, fano, quantum);
                for (i = 0; i < nChannels; i++)
                {
                    spectrum[i] += Math::hypermet(energy[i], \
                                                  area, position, fwhm, \
                                                  shortTailArea, shortTailSlope, \
                                                  longTailArea, longTailSlope, stepHeight);
                }
            }
        }
    }
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
