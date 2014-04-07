#include "fisx_xrf.h"
#include "fisx_math.h"
#include <cmath>
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
            distance += (*layerPtr).getThickness() / sinAlphaOut;
        }
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
        distance += detectorDistance;
        // calculate geometric efficiency 0.5 * (1 - cos theta)
        geometricEfficiency = 0.5 * (1.0 - (distance / sqrt(pow(distance, 2) + pow(0.5 * detectorDiameter, 2))));
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
    for (iLayer = 0; iLayer < sampleLayerIndex; iLayer++)
    {
        layerPtr = &sample[iLayer];
        if (iLayer < referenceLayerIndex)
        {
            distance += (*layerPtr).getThickness() / sinAlphaOut;
        }
    }
    // we can calculate the geometric efficiency for the given layer
    // calculate geometric efficiency 0.5 * (1 - cos theta)
    return (0.5 * (1.0 - (distance / sqrt(pow(distance, 2) + pow(0.5 * detectorDiameter, 2)))));
}

std::map<std::string, std::map<std::string, double> > XRF::getMultilayerFluorescence( \
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
    std::vector<double> & energies = actualRays[0];
    std::vector<double> weights;
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


    std::vector<std::vector<double> >sampleLayerEnergies;
    std::vector<std::vector<double> >sampleLayerRates;
    std::vector<std::vector<double> >sampleLayerMuTotal;
    std::map<std::string, double> sampleLayerComposition;
    std::vector<double>::size_type iLambda;
    std::vector<std::string> sampleLayerFamilies;
    std::vector<std::vector<std::pair<std::string, double> > >sampleLayerPeakFamilies;
    std::vector<std::pair<std::string, double> >::size_type iPeakFamily;
    std::vector<double> sampleLayerDensity;
    std::vector<double> sampleLayerThickness;

    sampleLayerEnergies.resize(sample.size());
    sampleLayerRates.resize(sample.size());
    sampleLayerMuTotal.resize(sample.size());
    sampleLayerDensity.resize(sample.size());
    sampleLayerThickness.resize(sample.size());

    iRay = energies.size() - 1;
    double weight;
    muTotal.resize(sample.size());
    while (iRay > 0)
    {
        --iRay;
        weights[iRay] = actualRays[1][iRay];
        // get the excitation factor for each layer at incident energy
        for(iLayer = 0; iLayer < sample.size(); iLayer++)
        {
            sampleLayerEnergies[iLayer].clear();
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
                tmpResult = elementsLibrary.getExcitationFactors(ele, energies[iRay],  weights[iRay]);
                // and add the energies and rates to the sampleLayerLines
                for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
                {
                    mapIt2 = c_it->second.find("energy");
                    sampleLayerEnergies[iLayer].push_back(mapIt2->second);
                    mapIt2 = c_it->second.find("rate");
                    sampleLayerRates[iLayer].push_back(mapIt2->second *
                                                       sampleLayerComposition[ele]);
                    sampleLayerFamilies[iLayer].push_back(iString);
                }
            }
            // calculate sample mu total at all those energies
            // Doubtfull. It has to be calculated for intermediate layers too.
            sampleLayerMuTotal[iLayer] = (*layerPtr).getMassAttenuationCoefficients( \
                                                            sampleLayerEnergies[iLayer], \
                                                            elementsLibrary)["total"];
            // get muTotal at the incident energy
            muTotal[iLayer] = (*layerPtr).getMassAttenuationCoefficients( \
                                                            energies[iRay], \
                                                            elementsLibrary)["total"];
            // layer thickness and density
            sampleLayerDensity[iLayer] = (*layerPtr).getDensity();
            sampleLayerThickness[iLayer] = (*layerPtr).getThickness();
        }

        // we start calculation by the bottom layer
        iLayer = sample.size();
        while(iLayer > 0)
        {
            --iLayer;
            // calculate excitation from top to bottom (case b in de Boer's article)
            for (jLayer = 0; jLayer < iLayer; jLayer++)
            {
                // get the mu, density and thickness b terms for each emission from j
                //mu_d_t_b = 0.0;
                for (int i = jLayer + 1; i < jLayer; i++)
                {
                    // for each fitted element and family
                    // get binding energy of family
                    // calculate contribution from all the emitted energies
                    /*
                    for (int iEnergy=0; iEnergy < sampleLayerEnergies[i]; iEnergy++)
                    {
                    //mu_d_t_b += sampleLayerDensity[i] * sampleLayerThickness[i] * \
                                    sampleLayerMuTotal[iEnergy];

                    }
                    */
                }
            }
        }


            /*
            layerPtr = &sample[iLayer];
            sampleLayerComposition[iLayer] = (*layerPtr).getComposition();
            for (mapIt = sampleLayerComposition.begin(); \
                 mapIt != sampleLayerComposition.end(); ++mapIt)
            {
                // get excitation factors for each element
                tmpResult = elementsLibrary.getExcitationFactors(mapIt->first,
                                                    energies[iRay], weights[iRay]);
                // and add the energies and rates to the sampleLayerLines
                for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
                {
                    mapIt2 = c_it->second.find("energy");
                    sampleLayerEnergies[iLayer].push_back(mapIt2->second);
                    mapIt2 = c_it->second.find("rate");
                    sampleLayerRates[iLayer].push_back(mapIt2->second * mapIt->second);
                }

                // the layer
                sampleLayerMuTotal[iLayer] = (*layerPtr).getMassAttenuationCoefficients( \
                                                            sampleLayerEnergies[iLayer], \
                                                            elementsLibrary)["total"];
            }
            */
    }
}
