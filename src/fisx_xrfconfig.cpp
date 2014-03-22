#include "fisx_xrfconfig.h"
#include "fisx_simpleini.h"
#include <stdexcept>
#include <iostream>
#include <algorithm>

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
    std::vector<double>::size_type iDoubleVector;
    std::vector<int> intVector, flagVector;
    std::vector<int>::size_type iIntVector;
    std::map<std::string, std::string>::const_iterator c_it;
    Material material;
    Layer layer;
    bool fisxFile;
    long counter;
    bool multilayerSample;
    double value;

    // find out if it is a fix or a PyMca configuration file
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
        return;
    }
    // Assume is a PyMca generated file.
    // TODO: Still to find out if it is a fir output file or a configuration file
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
            throw std::invalid_argument("File not recognized as a fisx or PyMca configuration file.");
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
            if(mapDoubles["weight"][iIntVector] <= 0.0)
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
                material = Material(stringVector[iStringVector], density, thickness, comment);
                material.setComposition(compoundList, compoundFractions);
                this->materials.push_back(material);
            }
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
            if (c_it->first.substr(0, 9) == "BeamFilter")
            {
                //BeamFilter0 = 0, -, 0.0, 0.0, 1.0
                layer = Layer(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                layer.setMaterial(stringVector[1]);
                this->beamFilters.push_back(layer);
            }
            else
            {
                // atmosphere = 0, -, 0.0, 0.0, 1.0
                // Matrix = 0, MULTILAYER, 0.0, 0.0, 45.0, 45.0, 0, 90.0
                if (stringVector.size() == 8 )
                {
                    // Matrix
                    this->setGeometry(doubleVector[5], doubleVector[6], doubleVector[7]);
                    if (stringVector[1] == "MULTILAYER")
                    {
                        std::cout << "MULTILAYER NOT PARSED YET" << std::endl;
                        multilayerSample = true;
                    }
                    else
                    {
                        layer = Layer(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                        layer.setMaterial(stringVector[1]);
                        this->sample.push_back(layer);
                    }
                }
                else
                {
                    if (stringVector[1] == "Detector")
                    {
                        // DETECTOR
                        std::cout << "DETECTOR " << std::endl;
                        this->detector = Detector(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                        this->detector.setMaterial(stringVector[1]);
                    }
                    else
                    {
                        // Attenuator
                        std::cout << "ATTENUATOR " << std::endl;
                        layer = Layer(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                        layer.setMaterial(stringVector[1]);
                        this->attenuators.push_back(layer);
                    }
                }
            }
        }
    }

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
                    layer = Layer(c_it->first, doubleVector[2], doubleVector[3], doubleVector[4]);
                    layer.setMaterial(stringVector[1]);
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
