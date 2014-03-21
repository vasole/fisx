#include "fisx_xrfconfig.h"
#include "fisx_simpleini.h"
#include <stdexcept>
#include <iostream>

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

void XRFConfig::importConfigurationFromFile(const std::string & fileName)
{
    SimpleIni iniFile = SimpleIni(fileName);
    std::map<std::string, std::string> sectionContents;
    std::string content;
    std::map<std::string, std::vector<double> > mapDoubles;
    std::map<std::string, std::vector<std::string> > mapStrings;
    std::vector<int> intVector, flagVector;
    std::vector<int>::size_type iIntVector;
    bool fisxFile;
    long counter;

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
    // Assume it is configuration file
    sectionContents.clear();
    sectionContents = iniFile.readSection("fit", false);
    if(!sectionContents.size())
    {
        throw std::invalid_argument("File not recognized as a fisx configuration file.");
    }
    mapDoubles.clear();
    content = sectionContents["energy"];
    iniFile.parseStringAsMultipleValues(content, mapDoubles["energy"], -1.0);
    content = sectionContents["energyweight"];
    iniFile.parseStringAsMultipleValues(content, mapDoubles["weight"], -1.0);
    content = sectionContents["energyscatter"];
    iniFile.parseStringAsMultipleValues(content, intVector, -1);
    content = sectionContents["energyflag"];
    iniFile.parseStringAsMultipleValues(content, flagVector, 0);
    counter = 0;
    for(iIntVector = 0; iIntVector < flagVector.size(); iIntVector++)
    {
        if (flagVector[iIntVector] > 0)
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
}

void XRFConfig::setBeam(const std::vector<double> & energy, \
                        const std::vector<double> & weight, \
                        const std::vector<int> & characteristic, \
                        const std::vector<double> & divergency)
{
    this->beam.setBeam(energy, weight, characteristic, divergency);
}
