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
#include "fisx_xrfconfig.h"
#include "fisx_simpleini.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>

namespace fisx
{

XRFConfig::XRFConfig()
{
    this->setGeometry(45.0, 45.0, 90.);
}

void XRFConfig::setGeometry(const double & alphaIn, const double & alphaOut, const double & scatteringAngle)
{
    this->alphaIn = alphaIn;
    this->alphaOut = alphaOut;
    this->scatteringAngle = scatteringAngle;
}

void XRFConfig::readConfigurationFromFile(const std::string & fileName)
{
    SimpleIni iniFile = SimpleIni(fileName);
    std::map<std::string, std::string> sectionContents;
    std::string key;
    std::string content;
    std::locale loc;
    std::map<std::string, std::vector<double> > mapDoubles;
    std::map<std::string, std::vector<std::string> > mapStrings;
    std::vector<std::string> stringVector;
    std::vector<std::string>::size_type iStringVector;
    std::vector<double> doubleVector;
    std::vector<int> intVector, flagVector;
    std::vector<int>::size_type iIntVector;
    std::map<std::string, std::string>::const_iterator c_it;
    std::vector<Material>::size_type iMaterial;
    std::vector<std::string> splitString;
    Material material;
    Layer layer;
    long counter;
    bool multilayerSample;
    double value;


    // find out if it is a fix or a PyMca configuration file
    /* This not needed since a fisx file format is not defined yet
    bool fisxFile;
    sectionContents.clear();
    sectionContents = iniFile.readSection("fisx", false);
    fisxFile = true;
    if(!sectionContents.size())
    {
        fisxFile = false;
    }
    if (fisxFile)
    {
        std::cout << "Not implemented" << std::endl;
        throw std::invalid_argument("Only PyMca configuration files supported.");
    }
    */

    // Assume is a PyMca generated file.
    // TODO: Still to find out if it is a fit output file or a configuration file
    // Assume it is configuration file.
    // In case of fit file, the configuration is under [result.config]
    // GET BEAM
    sectionContents.clear();
    sectionContents = iniFile.readSection("fit", false);
    if(!sectionContents.size())
    {
        sectionContents = iniFile.readSection("result.config", false);
        if(!sectionContents.size())
        {
            throw std::invalid_argument("File not recognized as a PyMca configuration file.");
        }
        std::cout << "fit result file" << std::endl;
    }
    else
    {
        // In case of PyMca.ini file, the configuration is under [Fit.Configuration]
        ;
    }
    mapDoubles.clear();
    content = sectionContents["energy"];
    iniFile.parseStringAsMultipleValues(content, mapDoubles["energy"], -666.0);
    content = sectionContents["energyweight"];
    iniFile.parseStringAsMultipleValues(content, mapDoubles["weight"], -1.0);
    content = sectionContents["energyscatter"];
    iniFile.parseStringAsMultipleValues(content, intVector, -1);
    content = sectionContents["energyflag"];
    iniFile.parseStringAsMultipleValues(content, flagVector, 0);

    /*
    std::cout << "Passed" << std::endl;
    std::cout << "Energy size = " << mapDoubles["energy"].size() << std::endl;
    std::cout << "weight size = " << mapDoubles["weight"].size() << std::endl;
    std::cout << "scatter = " << intVector.size() << std::endl;
    std::cout << "falg = " << flagVector.size() << std::endl;
    */

    if (mapDoubles["weight"].size() == 0)
    {
        mapDoubles["weight"].resize(mapDoubles["energy"].size());
        std::fill(mapDoubles["weight"].begin(), mapDoubles["weight"].end(), 1.0);
    }
    if (intVector.size() == 0)
    {
        intVector.resize(mapDoubles["energy"].size());
        std::fill(intVector.begin(), intVector.end(), 1.0);
    }
    if (flagVector.size() == 0)
    {
        flagVector.resize(mapDoubles["energy"].size());
        std::fill(flagVector.begin(), flagVector.end(), 1.0);
    }

    counter = 0;
    iIntVector = flagVector.size();
    while(iIntVector > 0)
    {
        iIntVector--;
        if ((flagVector[iIntVector] > 0) && (mapDoubles["energy"][iIntVector] != -666.0))
        {
            if(mapDoubles["energy"][iIntVector] <= 0.0)
            {
                throw std::invalid_argument("Negative excitation beam photon energy");
            }
            if(mapDoubles["weight"][iIntVector] < 0.0)
            {
                throw std::invalid_argument("Negative excitation beam photon weight");
            }
            if(intVector[iIntVector] < 0)
            {
                std::cout << "WARNING: " << "Negative characteristic flag. ";
                std::cout << "Assuming not a characteristic photon energy." << std::endl;
                intVector[iIntVector] = 0;
            }
            counter++;
        }
        else
        {
            // index not to be considered
            mapDoubles["energy"].erase(mapDoubles["energy"].begin() + iIntVector);
            mapDoubles["weight"].erase(mapDoubles["weight"].begin() + iIntVector);
            intVector.erase(intVector.begin() + iIntVector);
        }
    }
    this->setBeam(mapDoubles["energy"], mapDoubles["weight"], intVector);
    // GET THE MATERIALS
    iniFile.getSubsections("Materials", stringVector, false);
    this->materials.clear();
    if (stringVector.size())
    {
        std::string comment;
        double density;
        double thickness;
        std::vector<std::string> compoundList;
        std::vector<double> compoundFractions;
        std::string key;
        std::string::size_type j;
        // Materials found
        for (iStringVector = 0; iStringVector < stringVector.size(); iStringVector++)
        {
            sectionContents.clear();
            sectionContents = iniFile.readSection(stringVector[iStringVector], true);
            splitString.clear();
            iniFile.parseStringAsMultipleValues(stringVector[iStringVector],
                                            splitString,
                                            std::string(),
                                            '.');
            for (c_it = sectionContents.begin(); c_it != sectionContents.end(); ++c_it)
            {
                key = c_it->first;
                for (j = 0; j < key.size(); j++)
                {
                    key[j] = std::toupper(key[j], loc);
                }
                if (key == "DENSITY")
                {
                    iniFile.parseStringAsSingleValue(c_it->second, density, -1.0);
                }
                else if (key == "THICKNESS")
                {
                    iniFile.parseStringAsSingleValue(c_it->second, thickness, -1.0);
                }
                else if (key == "COMPOUNDLIST")
                {
                    iniFile.parseStringAsMultipleValues(c_it->second, \
                                                        compoundList, std::string());
                }
                else if ((key == "COMPOUNDFRACTION") || (key == "COMPOUNDFRACTIONS"))
                {
                    iniFile.parseStringAsMultipleValues(c_it->second, \
                                                        compoundFractions, -1.0);
                }
                else
                {
                    comment = c_it->second;
                }
            }
            material = Material(splitString[1], density, thickness, comment);
            material.setComposition(compoundList, compoundFractions);
            this->materials.push_back(material);
        }
    }

    // GET BEAM FILTERS AND ATTENUATORS
    sectionContents.clear();
    sectionContents = iniFile.readSection("attenuators", false);
    mapDoubles.clear();
    doubleVector.clear();
    stringVector.clear();
    this->beamFilters.clear();
    this->attenuators.clear();
    this->sample.clear();
    this->detector = Detector();
    multilayerSample = false;
    for (c_it = sectionContents.begin(); c_it != sectionContents.end(); ++c_it)
    {
        // std::cout << c_it->first << " " << c_it->second << std::endl;
        content = c_it->second;
        iniFile.parseStringAsMultipleValues(content, doubleVector, -1.0);
        iniFile.parseStringAsMultipleValues(content, stringVector, std::string());
        if (doubleVector.size() == 0.0)
        {
            std::cout << "WARNING: Empty line in attenuators section. Offending key is: "<< std::endl;
            std::cout << "<" << c_it->first << ">" << std::endl;
            continue;
        }
        if (doubleVector[0] > 0.0)
        {
            if (c_it->first.substr(0, c_it->first.size() - 1) == "BeamFilter")
            {
                // BeamFilter0 = 0, -, 0.0, 0.0, 1.0
                // std::cout << "BEAMFILTER" << std::endl;
                layer = Layer(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                layer.setMaterial(stringVector[1]);
                for(iMaterial = 0; iMaterial < this->materials.size(); iMaterial++)
                {
                    if(this->materials[iMaterial].getName() == stringVector[1])
                    {
                        layer.setMaterial(this->materials[iMaterial]);
                    }
                }
                this->beamFilters.push_back(layer);
            }
            else
            {
                // atmosphere = 0, -, 0.0, 0.0, 1.0
                // Matrix = 0, MULTILAYER, 0.0, 0.0, 45.0, 45.0, 0, 90.0
                if (stringVector.size() == 8 )
                {
                    // Matrix
                    if (doubleVector[6] > 0.0)
                    {
                        this->setGeometry(doubleVector[4], doubleVector[5], doubleVector[7]);
                    }
                    else
                    {
                        this->setGeometry(doubleVector[4], doubleVector[5], doubleVector[4] + doubleVector[5]);
                    }
                    if (stringVector[1] == "MULTILAYER")
                    {
                        multilayerSample = true;
                    }
                    else
                    {
                        // funny factor is not set for the sample
                        layer = Layer(c_it->first, doubleVector[2], doubleVector[3], 1.0);
                        layer.setMaterial(stringVector[1]);
                        for(iMaterial = 0; iMaterial < this->materials.size(); iMaterial++)
                        {
                            if(this->materials[iMaterial].getName() == stringVector[1])
                            {
                                layer.setMaterial(this->materials[iMaterial]);
                            }
                        }
                        this->sample.push_back(layer);
                    }
                }
                else
                {
                    if (c_it->first.substr(0, 8) == "Detector")
                    {
                        // DETECTOR
                        // std::cout << "DETECTOR " << std::endl;
                        this->detector = Detector(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                        this->detector.setMaterial(stringVector[1]);
                        for(iMaterial = 0; iMaterial < this->materials.size(); iMaterial++)
                        {
                            if(this->materials[iMaterial].getName() == stringVector[1])
                            {
                                detector.setMaterial(this->materials[iMaterial]);
                            }
                        }
                    }
                    else
                    {
                        // Attenuator
                        // std::cout << "ATTENUATOR " << std::endl;
                        layer = Layer(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                        layer.setMaterial(stringVector[1]);
                        for(iMaterial = 0; iMaterial < this->materials.size(); iMaterial++)
                        {
                            if(this->materials[iMaterial].getName() == stringVector[1])
                            {
                                layer.setMaterial(this->materials[iMaterial]);
                            }
                        }
                        this->attenuators.push_back(layer);
                    }
                }
            }
        }
    }
    // for the time being it is not yet in the file
    this->referenceLayer = 0;
    // GET MULTILAYER SAMPLE IF NEEDED
    if (multilayerSample)
    {
        sectionContents.clear();
        sectionContents = iniFile.readSection("multilayer", false);
        for (c_it = sectionContents.begin(); c_it != sectionContents.end(); ++c_it)
        {
            // std::cout << c_it->first << " " << c_it->second << std::endl;
            content = c_it->second;
            iniFile.parseStringAsMultipleValues(content, doubleVector, -1.0);
            iniFile.parseStringAsMultipleValues(content, stringVector, std::string());
            if (doubleVector.size() == 0.0)
            {
                std::cout << "WARNING: Empty line in multilayer section. Offending key is: "<< std::endl;
                std::cout << "<" << c_it->first << ">" << std::endl;
                continue;
            }
            if (doubleVector[0] > 0.0)
            {
                    //BeamFilter0 = 0, -, 0.0, 0.0, 1.0
                    layer = Layer(c_it->first, doubleVector[2], doubleVector[3]);
                    layer.setMaterial(stringVector[1]);
                    for(iMaterial = 0; iMaterial < this->materials.size(); iMaterial++)
                    {
                        if(this->materials[iMaterial].getName() == stringVector[1])
                        {
                            layer.setMaterial(this->materials[iMaterial]);
                        }
                    }
                    this->sample.push_back(layer);
            }
        }
    }

    // CONCENTATIONS SETUP
    sectionContents.clear();
    sectionContents = iniFile.readSection("concentrations", false);
    for (c_it = sectionContents.begin(); c_it != sectionContents.end(); ++c_it)
    {
        key = c_it->first;
        iniFile.toUpper(key);
        iniFile.parseStringAsSingleValue(c_it->second, value, -1.0);
        if (key == "DISTANCE")
        {
            this->detector.setDistance(value);
        }
        if (key == "AREA")
        {
            this->detector.setActiveArea(value);
        }
    }
}

void XRFConfig::setBeam(const std::vector<double> & energy, \
                        const std::vector<double> & weight, \
                        const std::vector<int> & characteristic, \
                        const std::vector<double> & divergency)
{
    this->beam.setBeam(energy, weight, characteristic, divergency);
}

void XRFConfig::setBeam(const double & energy, const double & divergency)
{
    this->beam.setBeam(energy, divergency);
}

void XRFConfig::setBeam(const Beam & beam)
{
    this->beam = beam;
}

const Beam & XRFConfig::getBeam() const
{
    return this->beam;
}


void XRFConfig::setBeamFilters(const std::vector<Layer> & filters)
{
    this->beamFilters = filters;
}

void XRFConfig::setSample(const std::vector<Layer> & layers, const int & referenceLayer)
{
    if (referenceLayer >= (int) layers.size())
    {
        throw std::invalid_argument("Reference layer must be smaller than number of layers");
    }
    this->sample = layers;
    this->referenceLayer = referenceLayer;
}

void XRFConfig::setAttenuators(const std::vector<Layer> & attenuators)
{
    this->attenuators = attenuators;
}

void XRFConfig::setDetector(const Detector & detector)
{
    this->detector = detector;
}

std::ostream& operator<< (std::ostream& o, XRFConfig const& config)
{
    std::vector<Layer>::size_type i;
    o << "BEAM" << std::endl;
    o << config.beam << std::endl;
    o << "BEAM FILTERS" << std::endl;
    for(i = 0; i < config.beamFilters.size(); i++)
    {
        o << config.beamFilters[i] << std::endl;
    }
    o << "SAMPLE" << std::endl;
    for(i = 0; i < config.sample.size(); i++)
    {
        o << config.sample[i] << std::endl;
    }
    o << "ATTENUATORS" << std::endl;
    for(i = 0; i < config.attenuators.size(); i++)
    {
        o << config.attenuators[i] << std::endl;
    }
    o << "DETECTOR" << std::endl;
    o << config.detector << std::endl;

    o << "GEOMETRY" << std::endl;
    o << "Alpha In(deg): "<< config.getAlphaIn() << std::endl;
    o << "Alpha In(deg): "<< config.getAlphaOut() << std::endl;
    return o;
}

} // namespace fisx
