#include "fisx_xrf.h"
#include <math.h>
#include <stdexcept>
#include <iostream>

XRF::XRF()
{
    // initialize geometry with default parameters
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

void XRF::setBeamFilters(const std::vector<Layer> &layers)
{
    this->recentBeam = true;
    this->configuration.setBeamFilters(layers);
}

void XRF::setSample(const std::vector<Layer> & layers, const int & referenceLayer)
{
    this->configuration.setSample(layers, referenceLayer);
}

void XRF::setAttenuators(const std::vector<Layer> & attenuators)
{
    this->configuration.setAttenuators(attenuators);
}

void XRF::setDetector(const Detector & detector)
{
    this->configuration.setDetector(detector);
}


/* This implementation tries to optimize things but becomes difficult to understand
std::map< std::string, std::map< std::string, std::map<std::string, std::map<std::string, double> > > >\
                XRF::getExpectedPrimaryEmission(const std::vector<std::string> & elementList,
                                                const Elements & elements)
{
    // get the actual beam reaching the sample
    const Beam & beam = this->configuration.getBeam();
    std::vector<std::vector<double> >actualRays = beam.getBeamAsDoubleVectors();
    std::vector<double> & energies = actualRays[0];
    std::vector<double> weights;
    std::vector<double> muTotal;
    std::vector<double> tmpDouble;
    //const std::vector<Layer> & attenuators = this->configuration.getAttenuators();
    std::vector<Layer> attenuators;
    std::vector<Layer>::size_type iLayer;
    Material material;
    std::vector<Ray>::size_type iRay;
    std::vector<std::string>::size_type iString;
    double argument;

    muTotal.resize(energies.size());
    std::fill(muTotal.begin(), muTotal.end(), 0.0);
    attenuators = this->configuration.getBeamFilters();
    for (iLayer = 0; iLayer < attenuators.size(); iLayer++)
    {
        if (attenuators[iLayer].hasMaterialComposition())
        {
            material = attenuators[iLayer].getMaterial();
            tmpDouble = elements.getMassAttenuationCoefficients( \
                                material.getComposition(), energies)["total"];
        }
        else
        {
           tmpDouble = elements.getMassAttenuationCoefficients(\
                                attenuators[iLayer].getMaterialName(), energies)["total"];
        }
        for (iRay = 0; iRay < tmpDouble.size(); iRay++)
        {
            muTotal[iRay] += (tmpDouble[iRay] * \
                               attenuators[iLayer].density * \
                               attenuators[iLayer].thickness);
        }
    }

    for(iRay = 0; iRay < actualRays.size(); iRay++)
    {
        actualRays[1][iRay] *= exp(- muTotal[iRay]) ;
    }

    // we have got the actual beam reaching the sample

    // For each layer
    // for each element
    // for each energy
    // calculate the number of primary vacancies (it is function of energy and weight)
    // calculate the propagation of those vacancies
    // get the emitted lines from the final vacancy distribution

    // std::map<std::string, double> initialVacancies;
    // std::map<std::string, std::map<std::string, double> > transitions;
    // A transition is specified by a key (IUPAC notation), a rate and an energy
    std::map<std::string, std::vector< std::map<std::string, std::map<std::string, double> > > > transitions;

    for (iString = 0; iString < elementList.size(); iString++)
    {
        const Element & elementObject = elements.getElement(elementList[iString]);
        // calculate the inital vacancy distribution for each energy
        transitions[elementList[iString]].resize(energies.size());
        for (iRay = 0; iRay < energies.size(); iRay++)
        {
            //initialVacancies = \
            //    elementObject.getInitialPhotoelectricVacancyDistribution(energies[iRay], cascade, useFluorescenceYield);
            // and the emitted x-rays energies and ratios already multiplied by the fluorescence yield*
            // taking into account the cascade
            transitions[elementList[iString]][iRay] = elementObject.getXRayLinesFromVacancyDistribution( \
                        elementObject.getInitialPhotoelectricVacancyDistribution(energies[iRay]), 1, true);
        }
    }


    // copy the rays reaching the sample into other vector that will be updated
    // as we go through the sample
    weights = actualRays[1];
    const std::vector<Layer> & sample = this->configuration.getSample();
    // map[element][sample_layer_number][transition][energy and rate]
    std::map<std::string, std::vector<std::map<std::string, std::map<std::string, double> > > > emittedRays;

    // create the needed space
    for (iString = 0; iString < elementList.size(); iString++)
    {
       emittedRays[elementList[iString]].resize(sample.size());
    }

    for (iLayer = 0; iLayer < sample.size(); iLayer++)
    {
        // get the layer total mass attenuation coefficients
        if (sample[iLayer].hasMaterialComposition())
        {
            material = sample[iLayer].getMaterial();
            muTotal = elements.getMassAttenuationCoefficients( \
                                material.getComposition(), energies)["total"];
        }
        else
        {
            muTotal = elements.getMassAttenuationCoefficients(\
                                attenuators[iLayer].getMaterialName(), energies)["total"];
        }
        // get the contribution for each element
        for (iString = 0; iString < elementList.size(); iString++)
        {
            const Element & elementObject = elements.getElement(elementList[iString]);
            std::vector<std::map<std::string, map<std::string, double> > > transitionsVector;
            std::map<std::string, map<std::string, double> > >::transitions;
            std::map<std::string, map<std::string, double> > >::iterator transitionsIt;
            transitionsVector = elements.getExcitationFactors(elementList[iString], energies, weights);
            // the rays were ordered by ascending energy, so the last one contains all the emitted lines
            std::map<std::string, double> muFluoTotal;
            transitions = transitionsVector[transitionsVector.size() - 1];
            tmpDouble.resize(transitions.size());
            for (transitionsIt = transitions.begin(); transitionsIt != transitionsIt.end(); ++transitionsIt)
            {
                transitionsIt->second;
            }

            // We have the total mass attenuation coefficient at each incident energy
            // We have the excitation factor at each emitted energy
            // We miss: - the sample total mass attenuation coefficients at each energy
            //          - the transmission through attenuators (and upper layers)
            //          - the detection efficiency

            for (iRay=0; iRay < energies.size(); iRay++)
            {
            }

            //elementObject.getShellConstants()

        }
        // update weights
        tmpDouble = sample[iLayer].getTransmission(energies, elements);
        for (iRay = 0; iRay < tmpDouble.size(); iRay++)
        {
            weights[iRay] *= tmpDouble[iRay];
        }
    }
}

*/

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
    const Detector & detector = this->configuration.getDetector();
    const Element & element = elementsLibrary.getElement(elementName);
    std::string msg;
    std::map<std::string, std::map<std::string, double> >  result;
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
    for (std::vector<Layer>::size_type iLayer = 0; iLayer < filters.size(); iLayer++)
    {
        layerPtr = &filters[iLayer];
        doubleVector = (*layerPtr).getTransmission(energies, elementsLibrary);
        for (iRay = 0; iRay < energies.size(); iRay++)
        {
            actualRays[1][iRay] *= doubleVector[iRay];
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
            distance += (*layerPtr).thickness / sinAlphaIn;
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
    std::map<std::string, std::map<std::string, double> > actualResult;
    std::map<std::string, std::map<std::string, double> >::const_iterator c_it;
    std::map<std::string, double>::const_iterator mapIt;
    std::map<std::string, double> muTotalFluo;
    std::map<std::string, double> detectionEfficiency;

    iRay = energies.size();
    while (iRay > 0)
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
            // calculate the transmission in the way back for each energy (it will be the same for each energy)
            // in the sample upper layers
            // in the attenuators (DO NOT FORGET the funny factor) Because of that it would be better that
            // the layer takes care of the calculation
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
        }
    }
}
