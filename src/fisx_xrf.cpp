#include "fisx_xrf.h"
#include "fisx_math.h"
#include "fisx_simpleini.h"
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

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
    const double & detectorDistance = detector.getDistance();
    const double & detectorDiameter = detector.getDiameter();
    double distance;
    double geometricEfficiency;
    const int & referenceLayerIndex = this->configuration.getReferenceLayer();
    double sinAlphaIn = sin(alphaIn*(PI/180.));
    double sinAlphaOut = sin(alphaOut*(PI/180.));
    double tmpDouble;
    std::string tmpString;

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
    double maxEnergy;

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
    double geometricEfficiency;
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
                const int & useGeometricEfficiency)
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
    const Element & element = elementsLibrary.getElement(elementName);
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
    std::cout << "element name = " << elementName << " " << lineFamily << std::endl;;

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
        primaryExcitationFactors = elementsLibrary.getExcitationFactors(elementName, \
                                                                    energies[iRay], \
                                                                    weights[iRay]);
        for (iLayer = 0; iLayer < sample.size(); iLayer++)
        {
            // here I should loop for all elements and families
            key = elementName + " " + lineFamily;
            // we need to calculate the layer mass attenuation coefficients at the fluorescent energies
            result.clear();
            for (c_it = primaryExcitationFactors.begin(); c_it != primaryExcitationFactors.end(); ++c_it)
            {
                if (c_it->first.compare(0, lineFamily.length(), lineFamily) == 0)
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
                tmpDouble = (1.0 - exp( - tmpDouble * density_1 * thickness_1)) / (tmpDouble * sinAlphaIn);
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
                if ((c_it->first == "KL2") && (iLayer == 0))
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
                    if (mapIt->first.compare(0, lineFamily.length(), lineFamily) == 0)
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
                                tmpDouble *= (0.5/sinAlphaIn);
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
                                    tmpDouble *= (0.5/sinAlphaIn);
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
                                    tmpDouble *= (0.5/sinAlphaIn);
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
                        }
                    }
                }
                actualResult[key][iLayer][c_it->first]["rate"] += (1.0 - totalEscape) * result[c_it->first]["rate"];
                actualResult[key][iLayer][c_it->first]["primary"] += result[c_it->first]["primary"];
                actualResult[key][iLayer][c_it->first]["secondary"] += result[c_it->first]["secondary"];
            }
        }
    }
    this->lastMultilayerFluorescence = actualResult;
    return actualResult;
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
                        std::cout << tmpString << std::cout;
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
                        std::cout << tmpString << std::cout;
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
                std::cout << tmpString << std::cout;
                throw std::invalid_argument(tmpString);
            }
            // TODO: Deal with Ka, Kb, L, L1, L2, L3, ...
            tmpString = tmpStringVector[0] + " " + tmpStringVector[1];
            emissionRatiosPointer = emissionRatios.find(tmpString);
            if (emissionRatiosPointer == emissionRatios.end())
            {
                tmpString = "Undefined emission ratios for element " + tmpStringVector[0] +\
                            " family " + tmpStringVector[1];
                std::cout << tmpString << std::cout;
                throw std::invalid_argument(tmpString);
            }
            // Emission ratios has the form "Cr K" but we have received peakFamily can have the form "Cr K index"
            // We have to to normalize the signal from that element, family and layer to unit area,
            // and multiply by the supplied area. This could have been already done ...
            layerIterator = emissionRatiosPointer->second.find(layerIndex);
            if (layerIterator == emissionRatiosPointer->second.end())
            {
                tmpString = "I do not have information for layer number " + tmpStringVector[2];
                std::cout << tmpString << std::cout;
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
                    std::cout << tmpString << std::cout;
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
                    std::cout << tmpString << std::cout;
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
