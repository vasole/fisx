#include "fisx_xrf.h"
#include <math.h>
#include <stdexcept>
#include <iostream>

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
    this->configuration.setGeometry(alphaIn, alphaOut, scatteringAngle);
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

std::map<std::string, std::map<std::string, double> > XRF::getFluorescence(const std::string elementName, \
                const Elements & elementsLibrary, const int & sampleLayerIndex, \
                const std::string & lineFamily, const int & secondary)
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
    std::map<std::string, std::map<std::string, double> >  result;
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
    distance = 0.0;
    for (iLayer = 0; iLayer < sampleLayerIndex; iLayer++)
    {
        layerPtr = &sample[iLayer];
        if (iLayer < referenceLayerIndex)
        {
            // TODO improve distance calculation
            distance += (*layerPtr).getThickness() / sinAlphaIn;
        }
        doubleVector = (*layerPtr).getTransmission(energies, \
                                            elementsLibrary, alphaIn);
        for (iRay = 0; iRay < energies.size(); iRay++)
        {
            weights[iRay] *= doubleVector[iRay];
        }
    }

    // we can already calculate the geometric efficiency
    if (detectorDiameter > 0.0)
    {
        distance += detectorDistance;
        // calculate geometric efficiency 0.5 * (1 - cos theta)
        geometricEfficiency = 0.5 * (1.0 - (detectorDiameter / sqrt(pow(distance, 2) + pow(detectorDiameter, 2))));
    }
    else
    {
        geometricEfficiency = 1.0 ;
    }

    // we have reached the layer we are interesed on
    // calculate its total mass attenuation coefficient at each incident energy
    std::fill(muTotal.begin(), muTotal.end(), 0.0);
    layerPtr = &sample[sampleLayerIndex];
    if ((*layerPtr).hasMaterialComposition())
    {
        const Material & material = (*layerPtr).getMaterial();
        doubleVector = elementsLibrary.getMassAttenuationCoefficients( \
                            material.getComposition(), energies)["total"];
    }
    else
    {
       doubleVector = elementsLibrary.getMassAttenuationCoefficients(\
                            (*layerPtr).getMaterialName(), energies)["total"];
    }
    for (iRay = 0; iRay < energies.size(); iRay++)
    {
        muTotal[iRay] = doubleVector[iRay] /sinAlphaIn;
    }

    std::map<std::string, std::map<std::string, double> > tmpResult;
    std::map<std::string, std::map<std::string, double> >::const_iterator c_it;
    std::map<std::string, double>::const_iterator mapIt;
    std::map<std::string, double> muTotalFluo;
    std::map<std::string, double> detectionEfficiency;

    iRay = energies.size();
    while (iRay > 0)
    {
        --iRay;
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
                detectionEfficiency[c_it->first] = 1.0;

                // transmission through upper layers
                iLayer = sampleLayerIndex;
                while (iLayer > 0)
                {
                    --iLayer;
                    layerPtr = &sample[sampleLayerIndex];
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

                // detection efficienty decomposed in geometric and intrinsic
                if (detectorDiameter > 0.0)
                {
                    // apply geometric efficiency 0.5 * (1 - cos theta)
                    detectionEfficiency[c_it->first] *= geometricEfficiency;

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
            tmpDouble = (1.0 - exp(- tmpDouble * \
                                   sample[sampleLayerIndex].getDensity() * \
                                   sample[sampleLayerIndex].getThickness())) / (tmpDouble * sinAlphaIn);
            mapIt = c_it->second.find("rate");
            actualResult[c_it->first]["rate"] += mapIt->second * tmpDouble * \
                                                detectionEfficiency[c_it->first];
        }
    }
    return actualResult;
}
