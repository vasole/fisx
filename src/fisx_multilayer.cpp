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
#include <cmath>
#include <stdexcept>
//#include <iostream>
#include <sstream>
#include <iomanip>

namespace fisx
{

std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                XRF::getMultilayerFluorescence(const std::vector<std::string> & elementList,
                                               const Elements & elementsLibrary, \
                                               const std::vector<int> & layerList, \
                                               const std::vector<std::string> &  familyList, \
                                               const int & secondary, \
                                               const int & useGeometricEfficiency,
                                               const int & useMassFractions, \
                                               const double & secondaryCalculationLimit)
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
    std::vector<Layer>::size_type jLayer;
    std::vector<Layer>::size_type bLayer;
    Detector detector = this->configuration.getDetector();
    std::string msg;
    const double PI = acos(-1.0);
    const double & alphaIn = this->configuration.getAlphaIn();
    const double & alphaOut = this->configuration.getAlphaOut();
    std::vector<double> geometricEfficiency;
    double sinAlphaIn = sin(alphaIn*(PI/180.));
    double sinAlphaOut = sin(alphaOut*(PI/180.));
    double tmpDouble;
    const std::vector<double> & energies = actualRays[0];
    std::vector<double> weights;
    std::vector<double> doubleVector;
    std::map<std::string, std::map<std::string, double> > result;
    std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > actualResult;
    std::vector<double> energyThresholdList;

    energyThresholdList.clear();
    // beam is ordered
    // maxEnergy = energies[energies.size() - 1];

    // get the beam after the beam filters
    std::vector<double> muTotal;
    muTotal.resize(energies.size());
    doubleVector.resize(energies.size());
    std::fill(muTotal.begin(), muTotal.end(), 0.0);
    for (iLayer = 0; iLayer < filters.size(); iLayer++)
    {
        doubleVector = filters[iLayer].getTransmission(energies, elementsLibrary);
        for (iRay = 0; iRay < energies.size(); iRay++)
        {
            actualRays[1][iRay] *= doubleVector[iRay];
        }
    }

    // we can already calculate the geometric efficiency
    geometricEfficiency.resize(sample.size());
    if (useGeometricEfficiency != 0)
    {
        for (iLayer = 0; iLayer < sample.size(); iLayer++)
        {
            geometricEfficiency[iLayer] = this->getGeometricEfficiency(iLayer);
        }
    }
    else
    {
        for (iLayer = 0; iLayer < sample.size(); iLayer++)
        {
            geometricEfficiency[iLayer] = 1.0;
        }
    }

    std::vector<std::vector<double> >sampleLayerEnergies;
    std::vector<std::vector<std::string> > sampleLayerEnergyNames;
    std::vector<std::vector<double> >sampleLayerRates;
    std::vector<std::vector<double> >sampleLayerMuTotal;
    std::map<std::string, double> sampleLayerComposition;
    std::vector<double>::size_type iLambda;
    std::vector<std::string> sampleLayerFamilies;
    std::vector<std::vector<std::pair<std::string, double> > >sampleLayerPeakFamilies;
    std::vector<std::pair<std::string, double> >::size_type iPeakFamily;
    std::vector<double> sampleLayerDensity;
    std::vector<double> sampleLayerThickness;
    std::vector<double> sampleLayerWeight;
    std::map< std::string, std::map<std::string, double> > escapeRates;

    // * implement a cache
    std::map< std::string, std::map< double, std::map<std::string, std::map<std::string, double> > > > \
                                        excitationFactorsCache;

    int updateEscape;
    updateEscape = 1;

    sampleLayerEnergies.resize(sample.size());
    sampleLayerEnergyNames.resize(sample.size());
    sampleLayerRates.resize(sample.size());
    sampleLayerFamilies.resize(sample.size());
    sampleLayerPeakFamilies.resize(sample.size());
    sampleLayerMuTotal.resize(sample.size());
    sampleLayerDensity.resize(sample.size());
    sampleLayerThickness.resize(sample.size());
    sampleLayerWeight.resize(sample.size());
    escapeRates.clear();

    iRay = energies.size();
    muTotal.resize(sample.size());
    weights.resize(actualRays[1].size());
    double minimumExcitationEnergy = -1.0;
    while (iRay > 0)
    {
        if (minimumExcitationEnergy < 0.0)
        {
            for (std::vector<std::string>::size_type iElement = 0; iElement < energyThresholdList.size(); iElement++)
            {
                if ((energyThresholdList[iElement] < minimumExcitationEnergy) ||
                    (minimumExcitationEnergy < 0.0))
                {
                    minimumExcitationEnergy = energyThresholdList[iElement];
                }
            }
        }
        --iRay;
        if (energies[iRay] < minimumExcitationEnergy)
        {
            // std::cout << "Stopped at Ray " << iRay << std::endl;
            continue;
        }
        weights[iRay] = actualRays[1][iRay];
        tmpDouble = 0.0;
        for(iLayer = 0; iLayer < sample.size(); iLayer++)
        {
            // get muTotal at the incident energy
            layerPtr = &sample[iLayer];
            if (iLayer == 0)
                sampleLayerWeight[iLayer] = 1.0;
            else
                sampleLayerWeight[iLayer] = exp(-tmpDouble);
            muTotal[iLayer] = (*layerPtr).getMassAttenuationCoefficients( \
                                                            energies[iRay], \
                                                            elementsLibrary)["total"];
            // layer thickness and density
            sampleLayerDensity[iLayer] = (*layerPtr).getDensity();
            sampleLayerThickness[iLayer] = (*layerPtr).getThickness();
            tmpDouble += sampleLayerDensity[iLayer] * sampleLayerThickness[iLayer] *\
                         muTotal[iLayer]/sinAlphaIn;
        }

        if (secondary > 0)
        {
            // get the excitation factor for each layer at incident energy
            for(iLayer = 0; iLayer < sample.size(); iLayer++)
            {
                sampleLayerEnergies[iLayer].clear();
                sampleLayerEnergyNames[iLayer].clear();
                sampleLayerRates[iLayer].clear();
                sampleLayerFamilies[iLayer].clear();
                sampleLayerMuTotal[iLayer].clear();
                layerPtr = &sample[iLayer];
                sampleLayerPeakFamilies[iLayer] = (*layerPtr).getPeakFamilies(energies[iRay], elementsLibrary);
                // They are ordered by increasing increasing binding energy
                std::string::size_type iString;
                std::string ele;
                std::string lastEle="dummytext";
                std::string family;
                std::map<std::string, std::map<std::string, double> > tmpResult;
                std::map<std::string, std::map<std::string, double> >::const_iterator c_it;
                std::map<std::string, double> sampleLayerComposition;
                std::map<std::string, double>::const_iterator mapIt;
                std::map<std::string, double>::const_iterator mapIt2;
                sampleLayerComposition = (*layerPtr).getComposition(elementsLibrary);
                for (iPeakFamily = 0 ; iPeakFamily < sampleLayerPeakFamilies[iLayer].size(); iPeakFamily++)
                {
                    iString = sampleLayerPeakFamilies[iLayer][iPeakFamily].first.find(' ');
                    family = sampleLayerPeakFamilies[iLayer][iPeakFamily].first.substr(iString + 1, \
                                    sampleLayerPeakFamilies[iLayer][iPeakFamily].first.size() - iString - 1);
                    ele = sampleLayerPeakFamilies[iLayer][iPeakFamily].first.substr(0, iString);
                    if (ele != lastEle)
                    {
                        // The secondary rates NOT corrected for the beam intensity reaching the layer
                        // The secondary rates NOT corrected for mass of element fraction in layer
                        tmpResult = elementsLibrary.getExcitationFactors(ele, \
                                                                     energies[iRay], \
                                                                     1.0);
                        lastEle = ele;
                    }
                    // and add the energies and rates to the sampleLayerLines
                    for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
                    {
                        // be carefull not to add twice an element
                        if (c_it->first.compare(0, family.length(), family) == 0)
                        {
                            mapIt2 = c_it->second.find("rate");
                            if ((mapIt2->second * sampleLayerComposition[ele]) <= 0.0)
                            {
                                // std::cout << " Lower equal zero " << mapIt2->second << std::endl;
                                continue;
                            }
                            if (secondaryCalculationLimit > 0.0)
                            {
                                if ((mapIt2->second * sampleLayerComposition[ele]) < secondaryCalculationLimit)
                                {
                                    continue;
                                };
                            }
                            mapIt2 = c_it->second.find("energy");
                            if (mapIt2->second < minimumExcitationEnergy)
                            {
                                continue;
                            }
                            sampleLayerEnergies[iLayer].push_back(mapIt2->second);
                            mapIt2 = c_it->second.find("rate");
                            // Store rates already corrected for the beam intensity reaching the layer
                            // Store rates already corrected for the element mass fraction
                            sampleLayerRates[iLayer].push_back(mapIt2->second * sampleLayerComposition[ele] * \
                                                               weights[iRay] * sampleLayerWeight[iLayer]);
                            sampleLayerFamilies[iLayer].push_back(iString);
                            sampleLayerEnergyNames[iLayer].push_back(ele + \
                                                                     " " +\
                                                                     c_it->first);
                        }
                    }
                }
                // We have to add the contribution of coherent scattering
                // We do so by assuming an isotropic emission of the same energy as the
                // incident beam
                sampleLayerEnergyNames[iLayer].push_back("coherent scattering");
                sampleLayerEnergies[iLayer].push_back(energies[iRay]);
                // calculate sample mu total at all those energies
                std::map<std::string, std::vector<double> > tmpStringDoubleVecMap;
                tmpStringDoubleVecMap = (*layerPtr).getMassAttenuationCoefficients( \
                                                                sampleLayerEnergies[iLayer], \
                                                                elementsLibrary);
                sampleLayerMuTotal[iLayer] = tmpStringDoubleVecMap["total"];

                sampleLayerRates[iLayer].push_back((weights[iRay] * sampleLayerWeight[iLayer])*\
                tmpStringDoubleVecMap["coherent"].back());
            }
        }
        // we start calculation
        // mu_1_lambda = Mass attenuation coefficient of iLayer at incident energy
        double mu_1_lambda;
        // mu_1_i = Mass attenuation coefficient of iLayer at fluorescent energy
        double mu_1_i;
        // density and thickness of fluorescent layer
        double density_1;
        double thickness_1;
        // density and thickness of second layer
        double density_2;
        double thickness_2;
        // mu_a_lambda = Mass attenuation coefficient of layers above iLayer at incident energy
        // mu_a_i = Mass attenuation coefficient of layers above iLayer at fluorescent energy
        // mu_b_lambda = Mass attenuation coefficient of layers between iLayer and jLayer at incident energy
        // mu_b_j = Mass attenuation coefficient of layers between iLayer and jLayer at jLayer fluorescent energy j
        // mu_b_j_d_t sum of product mu * density * thickness of layers between iLayer and jLayer at jLayer fluorescent energy j
        double mu_b_j_d_t;
        // mu_2_lambda = Mass attenuation coefficient of jLayer at incident energy
        double mu_2_lambda;
        // mu_2_j = Mass attenuation coefficient of jLayer at jLayer fluorescent energy j
        double mu_2_j;
        // mu_1_j = Mass attenuation coefficient of iLayer at jLayer fluorescent energy j
        double mu_1_j;
        // elementName is the element to be analyzed
        // lineFamily is the family of  elementName X-ray lines to be considered
        double energyThreshold;

        // this line can be moved out of the loop
        std::map<std::string, std::map<std::string, double> > primaryExcitationFactors;
        std::map<std::string, std::map<std::string, double> > tmpExcitationFactors;
        std::map<std::string, std::map<std::string, double> >::const_iterator c_it;
        std::map<std::string, double>::const_iterator mapIt;
        std::map<std::string, double> muTotalFluo;
        double detectionEfficiency;
        double energy;
        std::string key;
        std::string tmpString;
        std::ostringstream tmpStringStream;
        for (std::vector<std::string>::size_type iElement = 0; iElement < elementList.size(); iElement++)
        {
            const std::string & elementName = elementList[iElement];
            const std::string & lineFamily = familyList[iElement];
            int calculationLayer;
            std::string actualLineFamily;
            if (layerList.size() > 1)
                calculationLayer = layerList[iElement];
            else
                calculationLayer = layerList[0];
            actualLineFamily = lineFamily;
            if (lineFamily == "Ka")
            {
                actualLineFamily = "KL";
            }
            if (lineFamily == "Kb")
            {
                // carefull, the actual condition is to start by K and not to be followed by L
                actualLineFamily = "KM";
            }
            if (actualLineFamily == "")
            {
                throw std::runtime_error("All line families case not implemented yet!!!");
            }
            if (iElement < energyThresholdList.size())
            {
                energyThreshold = energyThresholdList[iElement];
            }
            else
            {
                energyThreshold = this->getEnergyThreshold(elementName, \
                                                       actualLineFamily.substr(0, 1), \
                                                       elementsLibrary);
                energyThresholdList.push_back(energyThreshold);
            }

            if (energyThreshold > energies[iRay])
            {
                continue;
            }
            primaryExcitationFactors = elementsLibrary.getExcitationFactors(elementName, \
                                                                        energies[iRay], \
                                                                        weights[iRay]);
            for (iLayer = 0; iLayer < sample.size(); iLayer++)
            {
                double elementMassFractionFactor;
                double elementMassFraction;
                if ((calculationLayer >= 0) && ((iLayer - calculationLayer) != 0))
                {
                    // no need to calculate this layer
                    continue;
                }
                elementMassFractionFactor = 1.0;
                if(true)
                {
                    std::map<std::string, double> sampleLayerComposition;
                    layerPtr = &sample[iLayer];
                    sampleLayerComposition = (*layerPtr).getComposition(elementsLibrary);
                    if (sampleLayerComposition.find(elementName) == sampleLayerComposition.end())
                    {
                        elementMassFraction = 0.0;
                    }
                    else
                    {
                        elementMassFraction = sampleLayerComposition[elementName];
                    }
                }
                if (useMassFractions)
                {
                    elementMassFractionFactor = elementMassFraction;
                }
                // here I should loop for all elements and families
                key = elementName + " " + lineFamily;
                // we need to calculate the layer mass attenuation coefficients at the fluorescent energies
                result.clear();
                if (elementMassFractionFactor == 0.0)
                    continue;
                for (c_it = primaryExcitationFactors.begin(); c_it != primaryExcitationFactors.end(); ++c_it)
                {
                    if ((c_it->first.compare(0, actualLineFamily.length(), actualLineFamily) == 0) || \
                        ((lineFamily == "Kb") && (c_it->first[0] == 'K') && (c_it->first[1] != 'L')))
                    {
                        mapIt = c_it->second.find("factor");
                        if (mapIt == c_it->second.end())
                        {
                            std::cout << "Key <factor> not found in excitation factor" << std::endl;
                        }
                        if (mapIt->second <= 0.0)
                        {
                            // not excited
                            continue;
                        }
                        mapIt = c_it->second.find("energy");
                        energy = mapIt->second;
                        result[c_it->first]["energy"] = energy;
                        mapIt = c_it->second.find("rate");
                        result[c_it->first]["rate"] = mapIt->second;
                        mapIt = c_it->second.find("factor");
                        result[c_it->first]["factor"] = mapIt->second;
                        if (actualResult[key][iLayer].find(c_it->first) == actualResult[key][iLayer].end())
                        {
                            // calculate layer mu total at fluorescent energy
                            // std::cout << "CALCULATING mu_1_i for " << c_it->first << " ";
                            // std::cout << "energy " << energy;
                            result[c_it->first]["mu_1_i"] = \
                                    sample[iLayer].getMassAttenuationCoefficients(energy, \
                                                                                elementsLibrary) ["total"];
                            // calculate detection efficiency of fluorescent energy
                            detectionEfficiency = 1.0;
                            // transmission through upper layers
                            jLayer = iLayer;
                            while (jLayer > 0)
                            {
                                jLayer--;
                                layerPtr = &sample[jLayer];
                                detectionEfficiency *= (*layerPtr).getTransmission(energy, \
                                                                                   elementsLibrary, \
                                                                                   alphaOut);
                            }
                            // transmission through attenuators
                            for (jLayer = 0; jLayer < attenuators.size(); jLayer++)
                            {
                                layerPtr = &attenuators[jLayer];
                                detectionEfficiency *= (*layerPtr).getTransmission(energy, \
                                                                                   elementsLibrary, \
                                                                                   90.0);
                            }

                            // detection efficiency decomposed in geometric and intrinsic
                            detectionEfficiency *= geometricEfficiency[iLayer];

                            if (detector.hasMaterialComposition() || (detector.getMaterialName().size() > 0 ))
                            {
                                if ((detector.getDensity() > 0.0) && (detector.getThickness() > 0.0))
                                {
                                    // calculate intrinsic efficiency
                                    // assuming normal incidence on detector surface
                                    detectionEfficiency *= (1.0 - detector.getTransmission(energy, \
                                                                                elementsLibrary, \
                                                                                90.0));
                                }
                                // calculate escape ratio assuming normal incidence on detector surface
                                escapeRates = detector.getEscape(energy, \
                                                                 elementsLibrary, \
                                                                 elementName + c_it->first, \
                                                                 updateEscape);
                                updateEscape = 0;
                            }


                            result[c_it->first]["energy_threshold"] = energyThreshold;
                            result[c_it->first]["efficiency"] = detectionEfficiency;
                            actualResult[key][iLayer][c_it->first]["efficiency"] = detectionEfficiency;
                            actualResult[key][iLayer][c_it->first]["energy"] = energy;
                            actualResult[key][iLayer][c_it->first]["energy_threshold"] = energyThreshold;
                            actualResult[key][iLayer][c_it->first]["mu_1_i"] = result[c_it->first]["mu_1_i"];
                            actualResult[key][iLayer][c_it->first]["rate"] = 0.0;
                            actualResult[key][iLayer][c_it->first]["primary"] = 0.0;
                            actualResult[key][iLayer][c_it->first]["secondary"] = 0.0;
                        }
                        else
                        {
                            // std::cout << "USING mu_1_i for " << c_it->first << " ";
                            // std::cout << "energy " << energy;
                            result[c_it->first]["efficiency"] = \
                                        actualResult[key][iLayer][c_it->first]["efficiency"];
                            result[c_it->first]["energy"] = actualResult[key][iLayer][c_it->first]["energy"];
                            result[c_it->first]["energy_threshold"] = \
                                                    actualResult[key][iLayer][c_it->first]["energy_threshold"];
                            result[c_it->first]["mu_1_i"] = actualResult[key][iLayer][c_it->first]["mu_1_i"];
                        }
                    }
                }
                if (result.size() == 0)
                {
                    // no need to calculate anything
                    continue;
                }
                // primary
                mu_1_lambda = sample[iLayer].getMassAttenuationCoefficients( \
                                                                energies[iRay], \
                                                                elementsLibrary)["total"];
                density_1 = sample[iLayer].getDensity();
                thickness_1 = sample[iLayer].getThickness();
                for (c_it = result.begin(); c_it != result.end(); ++c_it)
                {
                    mapIt = c_it->second.find("mu_1_i");
                    if (mapIt == c_it->second.end())
                    {
                        throw std::runtime_error("Mass attenuation coefficient not calculated!!!");
                    }
                    mu_1_i = mapIt->second;
                    tmpDouble = (mu_1_lambda / sinAlphaIn) + (mu_1_i / sinAlphaOut);
                    // keep factor for deciding if secondary excitation is to be considered or not
                    tmpDouble = (1.0 - exp( - tmpDouble * density_1 * thickness_1)) / tmpDouble;
                    //result[c_it->first]["criterium"] = tmpDouble;
                    tmpDouble *= (elementMassFractionFactor / sinAlphaIn);
                    result[c_it->first]["primary"] = tmpDouble * \
                                                     primaryExcitationFactors[c_it->first]["rate"] * \
                                                     sampleLayerWeight[iLayer];
                    result[c_it->first]["rate"] = result[c_it->first]["primary"] * \
                                                  result[c_it->first]["efficiency"];
                    result[c_it->first]["secondary"] = 0.0;
                    //std::cout << c_it->first << "efficiency = " << result[c_it->first]["efficiency"] << std::endl;
                    //std::cout << c_it->first << "primary = " << result[c_it->first]["primary"] << std::endl;
                    //std::cout << c_it->first << "energy = " << result[c_it->first]["energy"] << std::endl;
                    //std::cout << c_it->first << "mu_1_i = " << result[c_it->first]["mu_1_i"] << std::endl;
                    /*
                    if ((c_it->first == "KL2") && (iLayer == 0))
                    {
                        std::cout << c_it->first << "efficiency = " << result[c_it->first]["efficiency"] << std::endl;
                        std::cout << c_it->first << "primary = " << result[c_it->first]["primary"] << std::endl;
                        std::cout << c_it->first << "sampleLayerWeight = " << sampleLayerWeight[iLayer] << std::endl;
                        std::cout << c_it->first << "Excitation E = " << energies[iRay] << std::endl;
                        std::cout << c_it->first << "Fluorescence E = " << result[c_it->first]["energy"] << std::endl;
                        std::cout << c_it->first << "mu_1_i = " << mu_1_i << std::endl;
                        std::cout << c_it->first << "mu_1_lambda = " << mu_1_lambda<< std::endl;
                        std::cout << c_it->first << "d * t = " << density_1 * thickness_1 << std::endl;
                        std::cout << c_it->first << "mu_1_lambda/sinALphain = " << mu_1_lambda / sinAlphaIn<< std::endl;
                        std::cout << c_it->first << "mu_1_i/sinALphaOut = " << mu_1_i / sinAlphaOut<< std::endl;
                    }
                    */
                }

                if (secondary > 0)
                {
                    // calculate secondary
                    for (jLayer = 0; jLayer < sample.size(); jLayer++)
                    {
                        if (iLayer == jLayer)
                        {
                            double criteriumMax=secondaryCalculationLimit;
                            // intralayer secondary
                            for(iLambda = 0; iLambda < sampleLayerEnergies[jLayer].size(); iLambda++)
                            {
                                // analogous to incident beam
                                if (energyThreshold > sampleLayerEnergies[jLayer][iLambda])
                                    continue;
                                bool calculate;
                                calculate = true;
                                if (excitationFactorsCache.find(elementName) != excitationFactorsCache.end())
                                {
                                    if (excitationFactorsCache[elementName].find(sampleLayerEnergies[jLayer][iLambda]) \
                                            != excitationFactorsCache[elementName].end())
                                    {
                                        calculate = false;
                                    }
                                }
                                if (calculate)
                                {
                                    excitationFactorsCache[elementName] \
                                                [sampleLayerEnergies[jLayer][iLambda]] = \
                                                        elementsLibrary.getExcitationFactors(elementName, \
                                                        sampleLayerEnergies[jLayer][iLambda], \
                                                        1.0);
                                }
                                tmpExcitationFactors = excitationFactorsCache[elementName] \
                                                        [sampleLayerEnergies[jLayer][iLambda]];
                                double criterium;
                                for (c_it = result.begin(); c_it != result.end(); ++c_it)
                                {
                                    if (tmpExcitationFactors.find(c_it->first) == tmpExcitationFactors.end())
                                    {
                                        continue;
                                    }
                                    mapIt = result[c_it->first].find("mu_1_i");
                                    if (mapIt == result[c_it->first].end())
                                        throw std::runtime_error(" mu_1_i key. Mass attenuation not present???");
                                    criterium = tmpExcitationFactors[c_it->first]["rate"] * \
                                                    sampleLayerRates[jLayer][iLambda];
                                    criterium /= (primaryExcitationFactors[c_it->first]["rate"] * \
                                                     sampleLayerWeight[iLayer]);
                                    if (criterium <= criteriumMax)
                                    {
                                        // std::cout << " criterium = " << criterium;
                                        // std::cout << " criterium Max = " << criteriumMax << std::endl;
                                        continue;
                                    }
                                    mu_1_i = mapIt->second;
                                    tmpDouble = Math::deBoerL0(mu_1_lambda / sinAlphaIn,
                                                               mu_1_i / sinAlphaOut,
                                                               sampleLayerMuTotal[jLayer][iLambda],
                                                               density_1,
                                                               thickness_1);
                                    // Workaround incident angle of 90 degrees and scatter contribution
                                    if ((mu_1_lambda / sinAlphaIn) == sampleLayerMuTotal[jLayer][iLambda])
                                        tmpDouble += Math::deBoerL0(mu_1_i / sinAlphaOut,
                                                               mu_1_lambda / (0.99999*sinAlphaIn),
                                                               sampleLayerMuTotal[jLayer][iLambda],
                                                               density_1,
                                                               thickness_1);
                                    else
                                        tmpDouble += Math::deBoerL0(mu_1_i / sinAlphaOut,
                                                               mu_1_lambda / sinAlphaIn,
                                                               sampleLayerMuTotal[jLayer][iLambda],
                                                               density_1,
                                                               thickness_1);
                                    tmpDouble *= elementMassFractionFactor * (0.5/sinAlphaIn);
                                    tmpDouble *= tmpExcitationFactors[c_it->first]["rate"] * \
                                                    sampleLayerRates[jLayer][iLambda];
                                    tmpStringStream.str(std::string());
                                    tmpStringStream.clear();
                                    tmpStringStream << std::setfill('0') << std::setw(2) << jLayer;
                                    tmpString = sampleLayerEnergyNames[jLayer][iLambda] + " " + tmpStringStream.str();
                                    actualResult[elementName + " " + lineFamily][iLayer][c_it->first][tmpString] = \
                                                                                                    tmpDouble;
                                    result[c_it->first]["secondary"] += tmpDouble;
                                    result[c_it->first]["rate"] += tmpDouble * \
                                                                   result[c_it->first]["efficiency"];
                                    if (criteriumMax > 0.0)
                                    {
                                        if ((tmpDouble/result[c_it->first]["primary"]) < 0.00001)
                                        {
                                            if (criterium > criteriumMax)
                                            {
                                                criteriumMax = criterium;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            if ((sampleLayerWeight[jLayer] / sampleLayerWeight[iLayer]) < 1.0E-4)
                            {
                                // The incident reaching the layer originating the secondary excitation
                                // is relatively much weaker. In the common energy range of XRF interest,
                                // one will certainly recover less than of 50 % of the incident photons
                                // reaching the jLayer originating the secondary excitation "beam".
                                // An additional factor 100 is accounted for because of the possible ratio
                                // of attenuation coefficient of the element at incident and secondary excitation
                                // photon energies
                                continue;
                            }

                            double criterium;
                            // continue;
                            mu_2_lambda = muTotal[jLayer];
                            density_2 = sampleLayerDensity[jLayer];
                            thickness_2 = sampleLayerThickness[jLayer];
                            if (iLayer < jLayer)
                            {
                                double criteriumMax = secondaryCalculationLimit;
                                // interlayer case a)
                                for(iLambda = 0;
                                    iLambda < sampleLayerEnergies[jLayer].size(); \
                                    iLambda++)
                                {
                                    // analogous to incident beam
                                    energy = sampleLayerEnergies[jLayer][iLambda];
                                    if (energyThreshold > energy)
                                        continue;
                                    bool calculate;
                                    calculate = true;
                                    if (excitationFactorsCache.find(elementName) != excitationFactorsCache.end())
                                    {
                                        if (excitationFactorsCache[elementName].find(sampleLayerEnergies[jLayer][iLambda]) \
                                                != excitationFactorsCache[elementName].end())
                                        {
                                            calculate = false;
                                        }
                                    }
                                    if (calculate)
                                    {
                                        excitationFactorsCache[elementName] \
                                                    [sampleLayerEnergies[jLayer][iLambda]] = \
                                                            elementsLibrary.getExcitationFactors(elementName, \
                                                            sampleLayerEnergies[jLayer][iLambda], \
                                                            1.0);
                                    }
                                    tmpExcitationFactors = excitationFactorsCache[elementName] \
                                                            [sampleLayerEnergies[jLayer][iLambda]];
                                    for (c_it = result.begin(); c_it != result.end(); ++c_it)
                                    {
                                        if (tmpExcitationFactors.find(c_it->first) == tmpExcitationFactors.end())
                                        {
                                            // This happens when we look for K lines, but obviously L lines are
                                            // present
                                            //std::cout << "Not considered " << c_it->first ;
                                            //std::cout << " energy = " << sampleLayerEnergies[iLayer][iLambda] << std::endl;
                                            continue;
                                        }
                                        if (tmpExcitationFactors[c_it->first]["rate"] < 1.0e-30)
                                        {
                                            continue;
                                        }
                                        criterium = tmpExcitationFactors[c_it->first]["rate"] * \
                                                    sampleLayerRates[jLayer][iLambda];
                                        criterium /= (primaryExcitationFactors[c_it->first]["rate"] * \
                                                     sampleLayerWeight[iLayer]);
                                        if (criterium <= criteriumMax)
                                            continue;

                                        mapIt = result[c_it->first].find("mu_1_i");
                                        if (mapIt == result[c_it->first].end())
                                            throw std::runtime_error(" mu_1_i key. Mass attenuation not present???");
                                        mu_1_i = mapIt->second;
                                        mu_1_j = \
                                            sample[iLayer].getMassAttenuationCoefficients(energy, \
                                                                                elementsLibrary)["total"];
                                        mu_2_j = sampleLayerMuTotal[jLayer][iLambda];
                                        bLayer = iLayer + 1;
                                        mu_b_j_d_t = 0.0;
                                        while (bLayer < jLayer)
                                        {
                                            mu_b_j_d_t += sampleLayerDensity[bLayer] * \
                                                          sampleLayerThickness[bLayer] * \
                                                          sample[bLayer].getMassAttenuationCoefficients(energy, \
                                                                                    elementsLibrary)["total"];
                                            bLayer++;
                                        }
                                        tmpDouble = std::exp(-mu_1_i * density_1 * thickness_1/sinAlphaOut);
                                        if (tmpDouble < 0.001)
                                            continue;
                                        tmpDouble *= sampleLayerRates[jLayer][iLambda];
                                        if (-(mu_2_lambda/sinAlphaIn) == mu_2_j)
                                            tmpDouble *= Math::deBoerX(mu_2_lambda/(0.99999*sinAlphaIn), \
                                                                  mu_1_i/sinAlphaOut, \
                                                                  density_1 * thickness_1, \
                                                                  density_2 * thickness_2, \
                                                                  mu_1_j, \
                                                                  mu_2_j, \
                                                                  mu_b_j_d_t);
                                        else
                                            tmpDouble *= Math::deBoerX(mu_2_lambda/sinAlphaIn, \
                                                                  mu_1_i/sinAlphaOut, \
                                                                  density_1 * thickness_1, \
                                                                  density_2 * thickness_2, \
                                                                  mu_1_j, \
                                                                  mu_2_j, \
                                                                  mu_b_j_d_t);
                                        tmpDouble *= elementMassFractionFactor * (0.5/sinAlphaIn);
                                        tmpDouble *= tmpExcitationFactors[c_it->first]["rate"];
                                        tmpStringStream.str(std::string());
                                        tmpStringStream.clear();
                                        tmpStringStream << std::setfill('0') << std::setw(2) << jLayer;
                                        tmpString = sampleLayerEnergyNames[jLayer][iLambda] + " " + tmpStringStream.str();
                                        actualResult[elementName + " " + lineFamily][iLayer][c_it->first][tmpString] = \
                                                                                                        tmpDouble;
                                        result[c_it->first]["secondary"] += tmpDouble;
                                        result[c_it->first]["rate"] += tmpDouble * \
                                                                       result[c_it->first]["efficiency"];
                                        if (criteriumMax > 0.0)
                                        {
                                            if ((tmpDouble/result[c_it->first]["primary"]) < 0.00001)
                                            {
                                                if (criterium > criteriumMax)
                                                {
                                                    criteriumMax = criterium;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if (iLayer > jLayer)
                            {
                                double criteriumMax = secondaryCalculationLimit;
                                // interlayer case b)
                                double layerFactor;
                                layerFactor = std::exp(-mu_2_lambda * density_2 * thickness_2/sinAlphaIn);
                                /*
                                This test has been removed in favour of the one dealing with the relative
                                weights of each layer. Indeed, the relative amount of secondary excitation
                                can be huge important. Other point is if we are going to detect something.
                                if (layerFactor < 0.001)
                                {
                                    // No need to calculate anything, the top layer attenuates too much the
                                    // incoming beam
                                    continue;
                                }
                                */
                                for(iLambda = 0;
                                    iLambda < sampleLayerEnergies[jLayer].size(); \
                                    iLambda++)
                                {
                                    // analogous to incident beam
                                    energy = sampleLayerEnergies[jLayer][iLambda];
                                    if (energyThreshold > energy)
                                        continue;
                                    bool calculate;
                                    calculate = true;
                                    if (excitationFactorsCache.find(elementName) != excitationFactorsCache.end())
                                    {
                                        if (excitationFactorsCache[elementName].find(sampleLayerEnergies[jLayer][iLambda]) \
                                                != excitationFactorsCache[elementName].end())
                                        {
                                            calculate = false;
                                        }
                                    }
                                    if (calculate)
                                    {
                                        excitationFactorsCache[elementName] \
                                                    [sampleLayerEnergies[jLayer][iLambda]] = \
                                                            elementsLibrary.getExcitationFactors(elementName, \
                                                            sampleLayerEnergies[jLayer][iLambda], \
                                                            1.0);
                                    }
                                    tmpExcitationFactors = excitationFactorsCache[elementName] \
                                                            [sampleLayerEnergies[jLayer][iLambda]];
                                    for (c_it = result.begin(); c_it != result.end(); ++c_it)
                                    {
                                        if (tmpExcitationFactors.find(c_it->first) == tmpExcitationFactors.end())
                                        {
                                            // This happens when, for instance, we look for K lines, but obviously
                                            // L lines are present
                                            //std::cout << "Not considered " << c_it->first ;
                                            //std::cout << " energy = " << sampleLayerEnergies[iLayer][iLambda] << std::endl;
                                            continue;
                                        }
                                        if (tmpExcitationFactors[c_it->first]["rate"] < 1.0e-30)
                                        {
                                            continue;
                                        }
                                        criterium = tmpExcitationFactors[c_it->first]["rate"] * \
                                                    sampleLayerRates[jLayer][iLambda];
                                        criterium /= (primaryExcitationFactors[c_it->first]["rate"] * \
                                                     sampleLayerWeight[iLayer]);
                                        if (criterium <= criteriumMax)
                                            continue;
                                        mapIt = result[c_it->first].find("mu_1_i");
                                        if (mapIt == result[c_it->first].end())
                                            throw std::runtime_error(" mu_1_i key. Mass attenuation not present???");
                                        mu_1_i = mapIt->second;
                                        mu_1_j = \
                                            sample[iLayer].getMassAttenuationCoefficients(energy, \
                                                                                elementsLibrary)["total"];
                                        mu_2_j = sampleLayerMuTotal[jLayer][iLambda];
                                        bLayer = jLayer + 1;
                                        mu_b_j_d_t = 0.0;
                                        while (bLayer < iLayer)
                                        {
                                            mu_b_j_d_t += sampleLayerDensity[bLayer] * \
                                                          sampleLayerThickness[bLayer] * \
                                                          sample[bLayer].getMassAttenuationCoefficients(energy, \
                                                                                    elementsLibrary)["total"];
                                            bLayer++;
                                        }
                                        tmpDouble = layerFactor * sampleLayerRates[jLayer][iLambda];
                                        if ((mu_2_lambda/sinAlphaIn) == mu_2_j)
                                            tmpDouble *= Math::deBoerX(-mu_2_lambda/(0.99999*sinAlphaIn), \
                                                                  -mu_1_i/sinAlphaOut, \
                                                                  density_1 * thickness_1, \
                                                                  density_2 * thickness_2, \
                                                                  mu_1_j, \
                                                                  mu_2_j, \
                                                                  mu_b_j_d_t);
                                        else
                                            tmpDouble *= Math::deBoerX(-mu_2_lambda/sinAlphaIn, \
                                                                  -mu_1_i/sinAlphaOut, \
                                                                  density_1 * thickness_1, \
                                                                  density_2 * thickness_2, \
                                                                  mu_1_j, \
                                                                  mu_2_j, \
                                                                  mu_b_j_d_t);
                                        tmpDouble *= elementMassFractionFactor * (0.5/sinAlphaIn);
                                        tmpDouble *= tmpExcitationFactors[c_it->first]["rate"];
                                        tmpStringStream.str(std::string());
                                        tmpStringStream.clear();
                                        tmpStringStream << std::setfill('0') << std::setw(2) << jLayer;
                                        tmpString = sampleLayerEnergyNames[jLayer][iLambda] + " " + tmpStringStream.str();
                                        actualResult[elementName + " " + lineFamily][iLayer][c_it->first][tmpString] = \
                                                                                                        tmpDouble;
                                        result[c_it->first]["secondary"] += tmpDouble;
                                        result[c_it->first]["rate"] += tmpDouble * \
                                                                       result[c_it->first]["efficiency"];
                                        if (criteriumMax > 0.0)
                                        {
                                            if ((tmpDouble/result[c_it->first]["primary"]) < 0.00001)
                                            {
                                                if (criterium > criteriumMax)
                                                {
                                                    criteriumMax = criterium;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                // here we are done for the element and the layer
                key = elementName + " " + lineFamily;
                for (c_it = result.begin(); c_it != result.end(); ++c_it)
                {
                    double totalEscape = 0.0;
                    if (detector.hasMaterialComposition() || (detector.getMaterialName().size() > 0 ))
                    {
                        // calculate (if needed) escape ratio
                        escapeRates = detector.getEscape(result[c_it->first]["energy"], \
                                                         elementsLibrary, \
                                                         elementName + c_it->first, \
                                                         updateEscape);
                        if (escapeRates.size())
                        {
                            updateEscape = 0;
                            std::map<std::string, std::map<std::string, double> >::const_iterator c_it2;
     //                       std::map<std::string, double>::const_iterator mapIt;
                            for( c_it2 = escapeRates.begin(); c_it2!= escapeRates.end(); ++c_it2)
                            {
                                tmpString = c_it->first + " "+ c_it2->first;
                                if (actualResult[key][iLayer].find(tmpString) == actualResult[key][iLayer].end())
                                {
                                    mapIt = c_it2->second.find("energy");
                                    if (mapIt == c_it2->second.end())
                                    {
                                        throw std::runtime_error("Missing energy key in escape peak information!");
                                    }
                                    actualResult[key][iLayer][tmpString]["energy"] = mapIt->second;
                                    actualResult[key][iLayer][tmpString]["rate"] = 0.0;
                                    actualResult[key][iLayer][tmpString]["ratio"] = 0.0;
                                    actualResult[key][iLayer][tmpString]["primary"] = 0.0;
                                    actualResult[key][iLayer][tmpString]["secondary"] = 0.0;
                                }
                                mapIt = c_it2->second.find("rate");
                                if (mapIt == c_it2->second.end())
                                {
                                    throw std::runtime_error("Missing rate key in escape peak information!");
                                }
                                totalEscape += mapIt->second;
                                actualResult[key][iLayer][tmpString]["rate"] += mapIt->second * result[c_it->first]["rate"];
                                actualResult[key][iLayer][tmpString]["ratio"] = mapIt->second;
                                // The only meaning of filling "primary" and "secondary" for a escape peak is in order to
                                // be able to evaluate the ratio without having to refer to the actual parent line.
                                actualResult[key][iLayer][tmpString]["primary"] += mapIt->second * result[c_it->first]["primary"];
                                actualResult[key][iLayer][tmpString]["secondary"] += mapIt->second * result[c_it->first]["secondary"];
                            }
                        }
                    }
                    actualResult[key][iLayer][c_it->first]["rate"] += (1.0 - totalEscape) * result[c_it->first]["rate"];
                    // primary and secondary are the same independently of having escape or not.
                    actualResult[key][iLayer][c_it->first]["primary"] += result[c_it->first]["primary"];
                    actualResult[key][iLayer][c_it->first]["secondary"] += result[c_it->first]["secondary"];
                    actualResult[key][iLayer][c_it->first]["massFraction"] = elementMassFraction;
                }
            }
        }
    }
    if (secondary > 1)
    {
        // std::cout << "WARNING: Tertiary excitation under development " << std::endl;
        // approximate tertiary excitation
        // we ignore the case of excitation after double rayleigh/coherent scattering
        std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > >::iterator actualResultIt;
        std::map<std::string, std::map<std::string, double> >::iterator c_it;
        std::map<std::string, std::map<std::string, double> >::iterator it;
        std::map<std::string, std::map<std::string, double> >::iterator it2;
        std::map<std::string, double> contributingKeys;
        std::map<std::string, double>::iterator contributingKeysIt;
        std::string key;
        std::ostringstream tmpStringStream;
        double factorFirst;
        double tertiary;
        std::string ele;
        for (actualResultIt = actualResult.begin(); actualResultIt != actualResult.end(); ++actualResultIt)
        {
            ele = actualResultIt->first.substr(0, actualResultIt->first.find(' '));
            for (iLayer = 0; iLayer < actualResultIt->second.size(); iLayer++)
            {
                for (it = actualResultIt->second[iLayer].begin(); \
                     it != actualResultIt->second[iLayer].end(); ++it)
                {
                    //std::cout << it->first << std::endl;
                    if (it->first.find("esc") != std::string::npos)
                    {
                        // this is a escape line -> Ignore it
                        continue;
                    }
                    if (it->second["massFraction"] < 5.0E-3)
                    {
                        // one should be able to neglect tertiary excitation from elements
                        // with less than 1 % concentration
                        break;
                    }
                    factorFirst = 1.0;
                    if (it->second["primary"] > 0.0)
                    {
                        factorFirst = (it->second["primary"] + it->second["secondary"]) / \
                                     it->second["primary"];
                    }
                    if (factorFirst < 1.01)
                    {
                        // tertiary contribution should already be less than 1 % -> Ignore it
                        continue;
                    }
                    else
                    {
                        tmpStringStream.str(std::string());
                        tmpStringStream.clear();
                        tmpStringStream << std::setfill('0') << std::setw(2) << iLayer;
                        key = ele + " " + it->first + " " + tmpStringStream.str();
                        contributingKeys[key] = factorFirst;
                        // the factor will be the same for all lines starting by KL2, being escape or not
                    }
                }
            }
        }

        for (actualResultIt = actualResult.begin(); actualResultIt != actualResult.end(); ++actualResultIt)
        {
            for (iLayer = 0; iLayer < actualResultIt->second.size(); iLayer++)
            {
                for (it = actualResultIt->second[iLayer].begin(); \
                     it != actualResultIt->second[iLayer].end(); ++it)
                {
                    factorFirst = 1.0;
                    tertiary = 0.0;
                    if (it->second["primary"] > 0.0)
                    {
                        factorFirst = (it->second["primary"] + it->second["secondary"]) / \
                                     it->second["primary"];
                    }
                    if (factorFirst < 1.01)
                    {
                        // It had less than 1 % secondary, we assume tertiary will be even less
                        it->second["tertiary"] = 0.0;
                        continue;
                    }
                    if (it->first.find("esc") != std::string::npos)
                    {
                        // this is a escape line -> the contribution to be used is the one of
                        // the originating line that should already be calculated
                        key = it->first.substr(0, it->first.find(" "));
                        tertiary = (it->second["primary"] + it->second["secondary"]) * \
                                   (actualResultIt->second[iLayer][key]["tertiary"] /  \
                                   (actualResultIt->second[iLayer][key]["primary"] +   \
                                    actualResultIt->second[iLayer][key]["secondary"]));
                    }
                    else
                    {
                        for (contributingKeysIt = contributingKeys.begin(); \
                             contributingKeysIt != contributingKeys.end(); ++contributingKeysIt)
                        {
                            if (it->second.find(contributingKeysIt->first) != it->second.end())
                            {
                                tertiary += it->second[contributingKeysIt->first] * \
                                            (contributingKeysIt->second - 1.0);
                            }
                        }
                    }
                    it->second["tertiary"] = tertiary;
                    // update the total rate to account for primary, secondary and tertiary
                    // rate was equal to (primary + secondary) times a certain efficiency factor
                    // rate = A * (primary + secondary) therefore  now we must update the rate to
                    // account for tertiary rate = A * (primary + secondary + tertiary)
                    it->second["rate"] *= (tertiary + it->second["primary"] + it->second["secondary"]) \
                                          / (it->second["primary"] + it->second["secondary"]);
                }
            }
        }
    }
    this->lastMultilayerFluorescence = actualResult;
    return actualResult;
}

} // namespace fisx
