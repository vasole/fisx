#include "fisx_xrf.h"
#include "fisx_math.h"
#include <cmath>
#include <stdexcept>
//#include <iostream>
#include <sstream>
#include <iomanip>

std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                XRF::getMultilayerFluorescence(const std::vector<std::string> & elementList,
                                               const Elements & elementsLibrary, \
                                               const std::vector<int> & layerList, \
                                               const std::vector<std::string> &  familyList, \
                                               const int & secondary, \
                                               const int & useGeometricEfficiency,
                                               const int & useMassFractions)
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
    const double & detectorDistance = detector.getDistance();
    const double & detectorDiameter = detector.getDiameter();
    double distance;
    std::vector<double> geometricEfficiency;
    const int & referenceLayerIndex = this->configuration.getReferenceLayer();
    double sinAlphaIn = sin(alphaIn*(PI/180.));
    double sinAlphaOut = sin(alphaOut*(PI/180.));
    double tmpDouble;
    std::vector<double> & energies = actualRays[0];
    std::vector<double> weights;
    std::vector<double> doubleVector;
    double maxEnergy;
    std::map<std::string, std::map<std::string, double> > result;
    std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > actualResult;

    // beam is ordered
    maxEnergy = energies[energies.size() - 1];

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
    std::cout << "iRAY = " << iRay << std::endl;
    muTotal.resize(sample.size());
    weights.resize(actualRays[1].size());
    while (iRay > 0)
    {
        --iRay;
        weights[iRay] = actualRays[1][iRay];
        tmpDouble = 0.0;
        for(iLayer = 0; iLayer < sample.size(); iLayer++)
        {
            // get muTotal at the incident energy
            layerPtr = &sample[iLayer];
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
                    ele = sampleLayerPeakFamilies[iLayer][iPeakFamily].first.substr(0, iString);
                    family = sampleLayerPeakFamilies[iLayer][iPeakFamily].first.substr(iString + 1, \
                                    sampleLayerPeakFamilies[iLayer][iPeakFamily].first.size() - iString - 1);
                    // The secondary rates are already corrected for the beam intensity reaching the layer
                    tmpResult = elementsLibrary.getExcitationFactors(ele, \
                                                                     energies[iRay], \
                                                                     weights[iRay] * sampleLayerWeight[iLayer]);
                    // and add the energies and rates to the sampleLayerLines
                    for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
                    {
                        // be carefull not to add twice an element
                        if (c_it->first.compare(0, family.length(), family) == 0)
                        {
                            mapIt2 = c_it->second.find("energy");
                            sampleLayerEnergies[iLayer].push_back(mapIt2->second);
                            mapIt2 = c_it->second.find("rate");
                            sampleLayerRates[iLayer].push_back(mapIt2->second *
                                                               sampleLayerComposition[ele]);
                            sampleLayerFamilies[iLayer].push_back(iString);
                            sampleLayerEnergyNames[iLayer].push_back(ele + \
                                                                     " " +\
                                                                     c_it->first);
                        }
                    }
                }
                // calculate sample mu total at all those energies
                // Doubtfull. It has to be calculated for intermediate layers too.
                sampleLayerMuTotal[iLayer] = (*layerPtr).getMassAttenuationCoefficients( \
                                                                sampleLayerEnergies[iLayer], \
                                                                elementsLibrary)["total"];
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
        // density and thickness of intermediate layers
        double density_b;
        double thickness_b;
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
                actualLineFamily = "KM";
            }
            if (actualLineFamily == "")
            {
                throw std::runtime_error("All line families case not implemented yet!!!");
            }
            energyThreshold = this->getEnergyThreshold(elementName, \
                                                       actualLineFamily.substr(0, 1), \
                                                       elementsLibrary);
            if (energyThreshold > energies[iRay])
            {
                continue;
            }
            primaryExcitationFactors = elementsLibrary.getExcitationFactors(elementName, \
                                                                        energies[iRay], \
                                                                        weights[iRay]);
            for (iLayer = 0; iLayer < sample.size(); iLayer++)
            {
                double elementMassFraction;
                if ((calculationLayer > 0) && (iLayer != calculationLayer))
                {
                    // no need to calculate this layer
                    continue;
                }
                elementMassFraction = 1.0;
                if (useMassFractions)
                {
                    std::map<std::string, double> sampleLayerComposition;
                    layerPtr = &sample[iLayer];
                    sampleLayerComposition = (*layerPtr).getComposition(elementsLibrary);
                    if (sampleLayerComposition.find(elementName) == sampleLayerComposition.end())
                    {
                        // no need to calculate this layer
                        elementMassFraction = 0.0;
                    }
                    else
                    {
                        elementMassFraction = sampleLayerComposition[elementName];
                    }
                }
                // here I should loop for all elements and families
                key = elementName + " " + lineFamily;
                // we need to calculate the layer mass attenuation coefficients at the fluorescent energies
                result.clear();
                for (c_it = primaryExcitationFactors.begin(); c_it != primaryExcitationFactors.end(); ++c_it)
                {
                    if (c_it->first.compare(0, actualLineFamily.length(), actualLineFamily) == 0)
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
                        if (c_it->first.size() == 4)
                        {
                            energyThreshold = this->getEnergyThreshold(elementName, \
                                                                       c_it->first.substr(0, 2), \
                                                                       elementsLibrary);
                        }
                        else
                        {
                            energyThreshold = this->getEnergyThreshold(elementName, \
                                                                       c_it->first.substr(0, 1), \
                                                                       elementsLibrary);
                        }
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
                                // calculate intrinsic efficiency
                                // assuming normal incidence on detector surface
                                detectionEfficiency *= (1.0 - detector.getTransmission(energy, \
                                                                                       elementsLibrary, \
                                                                                       90.0));
                                // calculate escape ratio assuming normal incidence on detector surface
                                escapeRates = detector.getEscape(energy, \
                                                                 elementsLibrary, \
                                                                 c_it->first, \
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
                    tmpDouble = elementMassFraction * \
                        (1.0 - exp( - tmpDouble * density_1 * thickness_1)) / (tmpDouble * sinAlphaIn);
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
                    if (false && (c_it->first == "KL2") && (iLayer == 0))
                    {
                        std::cout << c_it->first << "efficiency = " << result[c_it->first]["efficiency"] << std::endl;
                        std::cout << c_it->first << "primary = " << result[c_it->first]["primary"] << std::endl;
                        std::cout << c_it->first << "sampleLayerWeight = " << sampleLayerWeight[iLayer] << std::endl;
                        std::cout << c_it->first << "Excitation E = " << energies[iRay] << std::endl;
                        std::cout << c_it->first << "FLuoresnce E = " << result[c_it->first]["energy"] << std::endl;
                        std::cout << c_it->first << "mu_1_i = " << mu_1_i << std::endl;
                        std::cout << c_it->first << "mu_1_lambda = " << mu_1_lambda<< std::endl;
                        std::cout << c_it->first << "d * t = " << density_1 * thickness_1 << std::endl;
                        std::cout << c_it->first << "mu_1_lambda/sinALphain = " << mu_1_lambda / sinAlphaIn<< std::endl;
                        std::cout << c_it->first << "mu_1_i/sinALphaOut = " << mu_1_i / sinAlphaOut<< std::endl;
                    }
                }

                /*
                // initialize output in case nothing comes out (element not excited)
                // this is to be done at the very end
                if (actualResult[key].find(iLayer) == actualResult[key].end())
                {
                    // element family was not excited;
                    continue;
                }
                if (actualResult[key].find(iLayer) == actualResult.end())
                {
                    // element family has not been excited (and it will never will) because
                    // we start by the highest energy
                    // we need to fill some defaults
                    std::map<std::string, double> tmpMap;
                    tmpMap = elementsLibrary.getElement(elementName).getEmittedXRayLines(1000.);
                    for (mapIt = tmpMap.begin(); mapIt != tmpMap.end(); ++mapIt)
                    {
                        if (mapIt->first.compare(0, actualLineFamily.length(), actualLineFamily) == 0)
                        {
                            actualResult[key][iLayer][mapIt->first]["energy"] = mapIt->second;
                            actualResult[key][iLayer][mapIt->first]["primary"] = 0.0;
                            actualResult[key][iLayer][mapIt->first]["rate"] = 0.0;
                        }
                    }
                }
                */
                if (secondary > 0)
                {
                    // calculate secondary
                    for (jLayer = 0; jLayer < sample.size(); jLayer++)
                    {
                        if (iLayer == jLayer)
                        {
                            // intralayer secondary
                            for(iLambda = 0; iLambda < sampleLayerEnergies[jLayer].size(); iLambda++)
                            {
                                // analogous to incident beam
                                tmpExcitationFactors = elementsLibrary.getExcitationFactors(elementName, \
                                            sampleLayerEnergies[jLayer][iLambda], \
                                            sampleLayerRates[jLayer][iLambda]);
                                for (c_it = result.begin(); c_it != result.end(); ++c_it)
                                {
                                    if (tmpExcitationFactors.find(c_it->first) == tmpExcitationFactors.end())
                                    {
                                        continue;
                                    }
                                    // I could put a higher limit here
                                    if (tmpExcitationFactors[c_it->first]["rate"] <= 1.0e-30)
                                        continue;
                                    mapIt = result[c_it->first].find("mu_1_i");
                                    if (mapIt == result[c_it->first].end())
                                        throw std::runtime_error(" mu_1_i key. Mass attenuation noy present???");
                                    mu_1_i = mapIt->second;
                                    tmpDouble = Math::deBoerL0(mu_1_lambda / sinAlphaIn,
                                                               mu_1_i / sinAlphaOut,
                                                               sampleLayerMuTotal[jLayer][iLambda],
                                                               density_1,
                                                               thickness_1);
                                    tmpDouble += Math::deBoerL0(mu_1_i / sinAlphaOut,
                                                               mu_1_lambda / sinAlphaIn,
                                                               sampleLayerMuTotal[jLayer][iLambda],
                                                               density_1,
                                                               thickness_1);
                                    tmpDouble *= elementMassFraction * (0.5/sinAlphaIn);
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
                                    /*
                                    if ((sampleLayerEnergyNames[iLayer][iLambda] == "Fe KL2") || \
                                        (sampleLayerEnergyNames[iLayer][iLambda] == "Fe KL3"))
                                    {
                                        std::cout << c_it->first << " FROM ";
                                        std::cout << sampleLayerEnergyNames[iLayer][iLambda] << std::endl;
                                        std::cout << "Enhancement = ";
                                        std::cout << tmpDouble/result[c_it->first]["primary"] << std::endl;
                                        std::cout << "mu1 " << mu_1_lambda / sinAlphaIn << std::endl;
                                        std::cout << "mu2 " << mu_1_i / sinAlphaOut << std::endl;
                                        std::cout << "muj " << sampleLayerMuTotal[iLayer][iLambda] << std::endl;
                                        std::cout << "L0 = " << Math::deBoerL0(mu_1_lambda / sinAlphaIn,
                                                               mu_1_i / sinAlphaOut,
                                                               sampleLayerMuTotal[iLayer][iLambda],
                                                               density_1,
                                                               thickness_1) << std::endl;
                                        std::cout << "mu1 " << mu_1_i / sinAlphaOut << std::endl;
                                        std::cout << "mu2 " << mu_1_lambda / sinAlphaIn << std::endl;
                                        std::cout << "muj " << sampleLayerMuTotal[iLayer][iLambda] << std::endl;
                                        std::cout << "L0 = " << Math::deBoerL0(mu_1_i / sinAlphaOut,
                                                               mu_1_lambda / sinAlphaIn,
                                                               sampleLayerMuTotal[iLayer][iLambda],
                                                               density_1,
                                                               thickness_1) << std::endl;
                                    }
                                    */
                                }
                            }
                        }
                        else
                        {
                            // continue;
                            mu_2_lambda = muTotal[jLayer];
                            density_2 = sampleLayerDensity[jLayer];
                            thickness_2 = sampleLayerThickness[jLayer];
                            if (iLayer < jLayer)
                            {
                                // interlayer case a)
                                for(iLambda = 0;
                                    iLambda < sampleLayerEnergies[jLayer].size(); \
                                    iLambda++)
                                {
                                    // analogous to incident beam
                                    energy = sampleLayerEnergies[jLayer][iLambda];
                                    tmpExcitationFactors = elementsLibrary.getExcitationFactors( \
                                                            elementName, \
                                                            energy, \
                                                            sampleLayerRates[jLayer][iLambda]);
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
                                        tmpDouble = Math::deBoerX(mu_2_lambda/sinAlphaIn, \
                                                                  mu_1_i/sinAlphaOut, \
                                                                  density_1 * thickness_1, \
                                                                  density_2 * thickness_2, \
                                                                  mu_1_j, \
                                                                  mu_2_j, \
                                                                  mu_b_j_d_t);
                                        tmpDouble *= elementMassFraction * (0.5/sinAlphaIn);
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
                                    }
                                }
                            }
                            if (iLayer > jLayer)
                            {
                                // interlayer case b)
                                for(iLambda = 0;
                                    iLambda < sampleLayerEnergies[jLayer].size(); \
                                    iLambda++)
                                {
                                    // analogous to incident beam
                                    energy = sampleLayerEnergies[jLayer][iLambda];
                                    tmpExcitationFactors = elementsLibrary.getExcitationFactors( \
                                                            elementName, \
                                                            energy, \
                                                            sampleLayerRates[jLayer][iLambda]);
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
                                        tmpDouble = Math::deBoerX(-mu_2_lambda/sinAlphaIn, \
                                                                  -mu_1_i/sinAlphaOut, \
                                                                  density_1 * thickness_1, \
                                                                  density_2 * thickness_2, \
                                                                  mu_1_j, \
                                                                  mu_2_j, \
                                                                  mu_b_j_d_t);
                                        tmpDouble *= elementMassFraction * (0.5/sinAlphaIn);
                                        tmpDouble *= tmpExcitationFactors[c_it->first]["rate"];
                                        tmpStringStream.str(std::string());
                                        tmpStringStream.clear();
                                        tmpStringStream << std::setfill('0') << std::setw(2) << iLayer;
                                        tmpString = sampleLayerEnergyNames[jLayer][iLambda] + " " + tmpStringStream.str();
                                        actualResult[elementName + " " + lineFamily][iLayer][c_it->first][tmpString] = \
                                                                                                        tmpDouble;
                                        result[c_it->first]["secondary"] += tmpDouble;
                                        result[c_it->first]["rate"] += tmpDouble * \
                                                                       result[c_it->first]["efficiency"];
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
                        escapeRates = detector.getEscape(energy, \
                                                         elementsLibrary, \
                                                         c_it->first, \
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
                                actualResult[key][iLayer][tmpString]["primary"] += mapIt->second * result[c_it->first]["primary"];
                                actualResult[key][iLayer][tmpString]["secondary"] += mapIt->second * result[c_it->first]["secondary"];
                            }
                        }
                    }
                    actualResult[key][iLayer][c_it->first]["rate"] += (1.0 - totalEscape) * result[c_it->first]["rate"];
                    actualResult[key][iLayer][c_it->first]["primary"] += result[c_it->first]["primary"];
                    actualResult[key][iLayer][c_it->first]["secondary"] += result[c_it->first]["secondary"];
                }
            }
        }
    }
    this->lastMultilayerFluorescence = actualResult;
    return actualResult;
}

