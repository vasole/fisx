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
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include "fisx_elements.h"

namespace fisx
{


const std::string Elements::defaultDataDir()
{
    char *path;

    path = std::getenv("FISX_DATA_DIR");
    if (path != NULL)
    {
        return std::string(path);
    }
    else
    {
        return FISX_DATA_DIR;
    }
}

Elements::Elements(std::string epdl97Directory, std::string bindingEnergiesFileName, std::string crossSectionsFile)
{
    // this is to simplify an initialization equivalent to that of PyMca:
    // - use EPDL97 to calculate partial photoelectric cross sections
    // - use a different binding energies file
    // - use different mass atenuation coefficients (i.e. XCOM)
    this->initialize(epdl97Directory, bindingEnergiesFileName);
    if (crossSectionsFile.size())
    {
        this->setMassAttenuationCoefficientsFile(crossSectionsFile);
    }
}

Elements::Elements(std::string epdl97Directory)
{
    // pure EPDL97 initialization
    if (epdl97Directory.size() < 1)
    {
        this->initialize(Elements::defaultDataDir(), "");
    }
    else
    {
        this->initialize(epdl97Directory, "");
    }
}

Elements::Elements(std::string directoryName, short pymca)
{
    std::string BINDING_ENERGIES="EADL97_BindingEnergies.dat";
    std::string joinSymbol;
    std::string bindingEnergies;
    std::string crossSections;

#ifdef _WIN32
    joinSymbol = "\\";
#elif _WIN64
    joinSymbol = "\\";
#else
    joinSymbol = "//";
#endif
    // check if directoryName already ends with the joinSymbol
    if (directoryName.substr(directoryName.size() - 1, 1) == joinSymbol)
    {
        joinSymbol = "";
    }

    if (pymca)
    {
        // binding energies
        bindingEnergies = directoryName + joinSymbol + "BindingEnergies.dat";
        const std::string K_SHELL_CONSTANTS_FILE = "KShellConstants.dat";
        const std::string L_SHELL_CONSTANTS_FILE = "LShellConstants.dat";
        const std::string M_SHELL_CONSTANTS_FILE = "MShellConstants.dat";
        const std::string K_SHELL_RADIATIVE_FILE = "KShellRates.dat";
        const std::string L_SHELL_RADIATIVE_FILE = "LShellRates.dat";
        const std::string M_SHELL_RADIATIVE_FILE = "MShellRates.dat";
        const std::string XCOM_FILE = "XCOM_CrossSections.dat";
        std::string xcom = directoryName + joinSymbol + XCOM_FILE;
        std::string tmpFileName;
        this->initialize(directoryName, bindingEnergies);
        this->setMassAttenuationCoefficientsFile(xcom);
        tmpFileName = directoryName + joinSymbol + K_SHELL_CONSTANTS_FILE;
        this->setShellConstantsFile("K", tmpFileName);
        tmpFileName = directoryName + joinSymbol + L_SHELL_CONSTANTS_FILE;
        this->setShellConstantsFile("L", tmpFileName);
        tmpFileName = directoryName + joinSymbol + M_SHELL_CONSTANTS_FILE;
        this->setShellConstantsFile("M", tmpFileName);
        tmpFileName = directoryName + joinSymbol + K_SHELL_RADIATIVE_FILE;
        this->setShellRadiativeTransitionsFile("K", tmpFileName);
        tmpFileName = directoryName + joinSymbol + L_SHELL_RADIATIVE_FILE;
        this->setShellRadiativeTransitionsFile("L", tmpFileName);
        tmpFileName = directoryName + joinSymbol + M_SHELL_RADIATIVE_FILE;
        this->setShellRadiativeTransitionsFile("M", tmpFileName);
    }
    else
    {
        // binding energies
        bindingEnergies = directoryName + joinSymbol + BINDING_ENERGIES;
        this->initialize(directoryName, "");
    }
}

void Elements::initialize(std::string epdl97Directory, std::string bindingEnergiesFile)
{
#include "fisx_defaultelementsinfo.h"
    const std::string K_SHELL_CONSTANTS_FILE = "EADL97_KShellConstants.dat";
    const std::string L_SHELL_CONSTANTS_FILE = "EADL97_LShellConstants.dat";
    const std::string M_SHELL_CONSTANTS_FILE = "EADL97_MShellConstants.dat";
    const std::string K_SHELL_RADIATIVE_FILE = "EADL97_KShellRadiativeRates.dat";
    const std::string L_SHELL_RADIATIVE_FILE = "EADL97_LShellRadiativeRates.dat";
    const std::string M_SHELL_RADIATIVE_FILE = "EADL97_MShellRadiativeRates.dat";
    const std::string K_SHELL_NONRADIATIVE_FILE = "EADL97_KShellNonradiativeRates.dat";
    const std::string L_SHELL_NONRADIATIVE_FILE = "EADL97_LShellNonradiativeRates.dat";
    const std::string M_SHELL_NONRADIATIVE_FILE = "EADL97_MShellNonradiativeRates.dat";
    std::string shell;
    int i, j;
    std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};
    std::string symbol;
    int    atomicNumber;
    std::vector<double> energyGrid;
    std::map<std::string, std::vector<double> > massAttenuationCoefficients;
    std::string joinSymbol;
    std::string filename;

    // Indicate we are going to configure everything
    this->shellConstantsFile["K"] = "";
    this->shellConstantsFile["L"] = "";
    this->shellConstantsFile["M"] = "";
    this->shellRadiativeTransitionsFile["K"] = "";
    this->shellRadiativeTransitionsFile["L"] = "";
    this->shellRadiativeTransitionsFile["M"] = "";
    this->shellNonradiativeTransitionsFile["K"] = "";
    this->shellNonradiativeTransitionsFile["L"] = "";
    this->shellNonradiativeTransitionsFile["M"] = "";

    // initialize EPDL97
    this->epdl97.setDataDirectory(epdl97Directory);

    // read a different binding energies file if requested
    if (bindingEnergiesFile.size() > 0)
    {
        this->epdl97.loadBindingEnergies(bindingEnergiesFile);
    }

    // TODO Read the XCOM cross sections instead of EPDL ones?
    // file XCOM_CrossSections.dat

    // create the default set of elements
    this->elementList.resize(N_PREDEFINED_ELEMENTS);
    for (i = 0; i < N_PREDEFINED_ELEMENTS; i++)
    {
        symbol = defaultElementsInfo[i].symbol;
        atomicNumber = defaultElementsInfo[i].z;
        this->elementList[i] = Element(symbol, atomicNumber);
        this->elementList[i].setAtomicMass(defaultElementsInfo[i].atomicMass);
        this->elementList[i].setBindingEnergies(epdl97.getBindingEnergies(atomicNumber));
        massAttenuationCoefficients = epdl97.getMassAttenuationCoefficients(atomicNumber);
        this->elementList[i].setMassAttenuationCoefficients(massAttenuationCoefficients["energy"],\
                                       massAttenuationCoefficients["photoelectric"],              \
                                       massAttenuationCoefficients["coherent"],                   \
                                       massAttenuationCoefficients["compton"],                    \
                                       massAttenuationCoefficients["pair"]);
        for (j = 0; j < 10; j++)
        {
            shell = shellList[j];
            this->elementList[i].setPartialPhotoelectricMassAttenuationCoefficients(shell,     \
                                                        massAttenuationCoefficients["energy"], \
                                                        massAttenuationCoefficients[shell]);
        }
        this->elementDict[symbol] = i;
    }

#ifdef _WIN32
    joinSymbol = "\\";
#elif _WIN64
    joinSymbol = "\\";
#else
    joinSymbol = "//";
#endif

    // check if epdl97Directory already ends with the joinSymbol
    if (epdl97Directory.substr(epdl97Directory.size() - 1, 1) == joinSymbol)
    {
        joinSymbol = "";
    }

    // Load the shell constants
    filename = epdl97Directory + joinSymbol + K_SHELL_CONSTANTS_FILE;
    this->setShellConstantsFile("K", filename);
    filename = epdl97Directory + joinSymbol + L_SHELL_CONSTANTS_FILE;
    this->setShellConstantsFile("L", filename);
    filename = epdl97Directory + joinSymbol + M_SHELL_CONSTANTS_FILE;
    this->setShellConstantsFile("M", filename);

    // Load the non-radiative transition ratios (alwys from epdl97)
    filename = epdl97Directory + joinSymbol + K_SHELL_NONRADIATIVE_FILE;
    this->setShellNonradiativeTransitionsFile("K", filename);
    filename = epdl97Directory + joinSymbol + L_SHELL_NONRADIATIVE_FILE;
    this->setShellNonradiativeTransitionsFile("L", filename);
    filename = epdl97Directory + joinSymbol + M_SHELL_NONRADIATIVE_FILE;
    this->setShellNonradiativeTransitionsFile("M", filename);
    // Load the radiative transition ratios
    filename = epdl97Directory + joinSymbol + K_SHELL_RADIATIVE_FILE;
    this->setShellRadiativeTransitionsFile("K", filename);
    filename = epdl97Directory + joinSymbol + L_SHELL_RADIATIVE_FILE;
    this->setShellRadiativeTransitionsFile("L", filename);
    filename = epdl97Directory + joinSymbol + M_SHELL_RADIATIVE_FILE;
    this->setShellRadiativeTransitionsFile("M", filename);
}

// Element handling
bool Elements::isElementNameDefined(const std::string & elementName) const
{
    if (elementName.size() == 0)
    {
        return false;
    }
    if (this->elementDict.find(elementName) == this->elementDict.end())
    {
        return false;
    }
    return true;
}

std::vector<std::string> Elements::getElementNames()
{
    std::vector<Element>::size_type i;
    std::vector<std::string> result;

    result.resize(this->elementList.size());
    for (i = 0; i < this->elementList.size(); i++)
    {
        result[i] = this->elementList[i].getName();
    }

    return result;
}

const Element & Elements::getElement(const std::string & elementName) const
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i];
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}

Element Elements::getElementCopy(const std::string & elementName)
{
    if (this->isElementNameDefined(elementName))
        return this->elementList[this->elementDict[elementName]];
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}

void Elements::addElement(const Element & element)
{
    std::string name;
    name = element.getName();

    if (this->elementDict.find(name) != this->elementDict.end())
    {
        // an element with that name already exists
        this->elementList[this->elementDict[name]] = element;
    }
    else
    {
        this->elementDict[name] = this->elementList.size();
        this->elementList.push_back(element);
    }
}
// Shell constants
void Elements::setShellConstantsFile(const std::string & mainShellName, \
                                     const std::string & fileName)
{
    SimpleSpecfile sf;
    int    nScans, i;
    std::vector<std::string> tmpLabels;
    std::vector<std::string>::size_type nLabels, j;
    std::vector<std::string> subShells;
    std::vector<std::vector<double> > tmpValues;
    std::vector<std::vector<double> >::size_type n;
    std::map<std::string, double > tmpDict;
    std::string msg;

    if ((mainShellName == "K") || (mainShellName == "L") || (mainShellName == "M"))
    {
        // We have received a valid main shell and not a subshell
        ;
    }
    else
    {
        throw std::invalid_argument("Invalid main shell <" + mainShellName +">");
    }

    sf = SimpleSpecfile(fileName);
    nScans = sf.getNumberOfScans();
    if (mainShellName == "K")
    {
        if (nScans != 1)
        {
            msg = "Number of scans not equal one in K shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("K");
    }
    if (mainShellName == "L")
    {
        if (nScans != 3)
        {
            msg = "Number of scans not equal three in L shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("L1");
        subShells.push_back("L2");
        subShells.push_back("L3");
    }
    if (mainShellName == "M")
    {
        if (nScans != 5)
        {
            msg = "Number of scans not equal five in M shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("M1");
        subShells.push_back("M2");
        subShells.push_back("M3");
        subShells.push_back("M4");
        subShells.push_back("M5");
    }

    for (i = 0; i < nScans ; i++)
    {
        tmpLabels = sf.getScanLabels(i);
        tmpValues = sf.getScanData(i);
        nLabels = tmpLabels.size();
        if (tmpValues[i].size() != nLabels)
        {
            msg = "Number of values does not match number of labels in " + fileName;
            throw std::ios_base::failure(msg);
        }
        tmpDict.clear();
        for (n = 0; n < tmpValues.size(); n++)
        {
            for (j = 0; j < nLabels; j++)
            {
                if (tmpLabels[j] != "Z")
                {
                    if(tmpLabels[j].size() > 4)
                    {
                        // only omegaX, fXY supported
                        if (tmpLabels[j].substr(0, 5) == "omega")
                        {
                            tmpDict["omega"] = tmpValues[n][j];
                        }
                        else
                        {
                            //std::cout << "Label " << tmpLabels[j] << " ignored" << std::endl;
                            ;
                        }
                    }
                    else
                    {
                        tmpDict[tmpLabels[j]] = tmpValues[n][j];
                    }
                }

            }
            this->elementList[n].setShellConstants(subShells[i], tmpDict);
        }
    }
    this->shellConstantsFile[mainShellName] = fileName;
}

void Elements::setShellNonradiativeTransitionsFile(const std::string & mainShellName, \
                                                   const std::string & fileName)
{
    SimpleSpecfile sf;
    int    nScans, i;
    std::vector<std::string> tmpLabels;
    std::vector<std::string>::size_type nLabels, j;
    std::vector<std::string> subShells;
    std::vector<std::vector<double> > tmpValues;
    std::vector<std::vector<double> >::size_type n;
    std::map<std::string, double > tmpDict;
    std::map<std::string, double > bindingE;
    std::string subshell;
    std::string msg;

    if ((mainShellName == "K") || (mainShellName == "L") || (mainShellName == "M"))
    {
        // We have received a valid main shell and not a subshell
        ;
    }
    else
    {
        msg = "Invalid main shell <" + mainShellName +">";
        throw std::invalid_argument(msg);
    }

    sf = SimpleSpecfile(fileName);
    nScans = sf.getNumberOfScans();
    if (mainShellName == "K")
    {
        if (nScans != 1)
        {
            msg = "Number of scans not equal one in K shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("K");
    }
    if (mainShellName == "L")
    {
        if (nScans != 3)
        {
            msg = "Number of scans not equal three in L shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("L1");
        subShells.push_back("L2");
        subShells.push_back("L3");
    }
    if (mainShellName == "M")
    {
        if (nScans != 5)
        {
            msg = "Number of scans not equal five in M shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("M1");
        subShells.push_back("M2");
        subShells.push_back("M3");
        subShells.push_back("M4");
        subShells.push_back("M5");
    }

    for (i = 0; i < nScans ; i++)
    {
        tmpLabels = sf.getScanLabels(i);
        tmpValues = sf.getScanData(i);
        nLabels = tmpLabels.size();
        if (tmpValues[i].size() != nLabels)
        {
            msg = "Number of values does not match number of labels in " + fileName;
            throw std::ios_base::failure(msg);
        }
        tmpDict.clear();
        for (n = 0; n < tmpValues.size(); n++)
        {
            for (j = 0; j < nLabels; j++)
            {
                if (tmpLabels[j] != "Z")
                {
                    tmpDict[tmpLabels[j]] = tmpValues[n][j];
                }

            }
            bindingE = this->elementList[n].getBindingEnergies();
            if(bindingE[subShells[i]] > 0.0)
            {
                this->elementList[n].setNonradiativeTransitions(subShells[i], tmpDict);
            }
        }
    }
    this->shellNonradiativeTransitionsFile[mainShellName] = fileName;
}

void Elements::setShellRadiativeTransitionsFile(const std::string & mainShellName, \
                                                const std::string & fileName)
{
    SimpleSpecfile sf;
    int    nScans, i;
    std::vector<std::string> tmpLabels;
    std::vector<std::string>::size_type nLabels, j;
    std::vector<std::string> subShells;
    std::vector<std::vector<double> > tmpValues;
    std::vector<std::vector<double> >::size_type n;
    std::map<std::string, double > tmpDict;
    std::map<std::string, double > bindingE;
    std::string subshell;
    std::string msg;

    if ((mainShellName == "K") || (mainShellName == "L") || (mainShellName == "M"))
    {
        // We have received a valid main shell and not a subshell
        ;
    }
    else
    {
        msg = "Invalid main shell <" + mainShellName +">";
        throw std::invalid_argument(msg);
    }

    sf = SimpleSpecfile(fileName);
    nScans = sf.getNumberOfScans();
    if (mainShellName == "K")
    {
        if (nScans != 1)
        {
            msg = "Number of scans not equal one in K shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("K");
    }
    if (mainShellName == "L")
    {
        if (nScans != 3)
        {
            msg = "Number of scans not equal three in L shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("L1");
        subShells.push_back("L2");
        subShells.push_back("L3");
    }
    if (mainShellName == "M")
    {
        if (nScans != 5)
        {
            msg = "Number of scans not equal five in M shell constants file " + \
                   fileName;
            throw std::ios_base::failure(msg);
        }
        subShells.push_back("M1");
        subShells.push_back("M2");
        subShells.push_back("M3");
        subShells.push_back("M4");
        subShells.push_back("M5");
    }

    for (i = 0; i < nScans ; i++)
    {
        tmpLabels = sf.getScanLabels(i);
        tmpValues = sf.getScanData(i);
        nLabels = tmpLabels.size();
        if (tmpValues[i].size() != nLabels)
        {
            msg = "Number of values does not match number of labels in " + fileName;
            throw std::ios_base::failure(msg);
        }
        tmpDict.clear();
        for (n = 0; n < tmpValues.size(); n++)
        {
            for (j = 0; j < nLabels; j++)
            {
                if (tmpLabels[j] != "Z")
                {
                    tmpDict[tmpLabels[j]] = tmpValues[n][j];
                }

            }
            bindingE = this->elementList[n].getBindingEnergies();
            if(bindingE[subShells[i]] > 0.0)
            {
                this->elementList[n].setRadiativeTransitions(subShells[i], tmpDict);
            }
        }
    }
    this->shellRadiativeTransitionsFile[mainShellName] = fileName;
}

// Mass attenuation handling
void Elements::setMassAttenuationCoefficientsFile(const std::string & fileName)
{
    int n, nScans;
    SimpleSpecfile sf;
    std::vector<std::string> labels;
    std::vector<std::string>::size_type i, j;
    std::vector<std::vector<double> > data;
    std::vector<double> muEnergy;
    std::vector<double> muCompton;
    std::vector<double> muCoherent;
    std::vector<double> muPair;
    std::vector<double> muPhotoelectric;
    std::string key;

    sf = SimpleSpecfile(fileName);
    nScans = sf.getNumberOfScans();
    if (nScans < 1)
    {
        throw std::ios_base::failure("No scans found in file!");
    }

    for (n = 0; n < nScans; n++)
    {
        labels = sf.getScanLabels(n);
        data = sf.getScanData(n);
        muEnergy.clear();
        muPhotoelectric.clear();
        muCoherent.clear();
        muCompton.clear();
        muPair.clear();
        for (i = 0; i < labels.size(); i++)
        {
            key = this->epdl97.toUpperCaseString(labels[i]);
            if (key.find("ENERGY") != std::string::npos)
            {
                muEnergy.resize(data.size());
                for (j = 0; j < data.size(); j++)
                {
                    muEnergy[j] = data[j][i];
                }
            }
            if (key.find("PHOTO") != std::string::npos)
            {
                muPhotoelectric.resize(data.size());
                for (j = 0; j < data.size(); j++)
                {
                    muPhotoelectric[j] = data[j][i];
                }
            }
            if (key.find("PAIR") != std::string::npos)
            {
                muPair.resize(data.size());
                for (j = 0; j < data.size(); j++)
                {
                    muPair[j] = data[j][i];
                }
            }
            if (key.find("COMPTON") != std::string::npos)
            {
                muCompton.resize(data.size());
                for (j = 0; j < data.size(); j++)
                {
                    muCompton[j] = data[j][i];
                }
            }
            if (key.find("RAYLEIGH") != std::string::npos)
            {
                muCoherent.resize(data.size());
                for (j = 0; j < data.size(); j++)
                {
                    muCoherent[j] = data[j][i];
                }
            }
            if (key.find("COHERENT") != std::string::npos)
            {
                if (key.find("INCOHERENT") == std::string::npos)
                {
                    muCoherent.resize(data.size());
                    for (j = 0; j < data.size(); j++)
                    {
                        muCoherent[j] = data[j][i];
                    }
                }
            }
        }
        // TODO check the sizes of all of them (except pair) are not zero
        // muEnergy.clear();
        // muPhotoelectric.clear();
        // muCoherent.clear();
        // muCompton.clear();
        // muPair.clear();
        this->setMassAttenuationCoefficients(this->elementList[n].getName(), \
                                             muEnergy, \
                                             muPhotoelectric, \
                                             muCoherent, \
                                             muCompton, \
                                             muPair);
    }
}

void Elements::setMassAttenuationCoefficients(const std::string & name,
                                              const std::vector<double> & energy,
                                              const std::vector<double> & photoelectric,
                                              const std::vector<double> & coherent,
                                              const std::vector<double> & compton,
                                              const std::vector<double> & pair)
{
    std::map<std::string, std::vector<double> > massAttenuationCoefficients;
    std::map<std::string, std::vector<double> >::iterator it;
    std::map<std::string, std::pair<double, int> > extractedEdgeEnergies;
    std::vector<std::vector<double>::size_type> shellIndexList;
    std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};
    Element *element;
    std::vector<double> newEnergyGrid;
    std::vector<double>::size_type i, j, shellIndex;
    std::vector<double> oldTotalPhotoelectric;
    std::string msg;
    std::string shell;
    double tmpDouble;
    int atomicNumber, idx;


    if (this->elementDict.find(name) == this->elementDict.end())
    {
        msg = "Name " + name + " not among defined elements";
        throw std::invalid_argument(msg);
    }


    // we have to make sure the partial mass attenuation photoelectric coefficients and the total
    // photoelectric mass attenuation coefficients are self-consistent
    // we are doing this because the only compilation providing subshell photoelectric cross sections
    // is EPDL97

    // we reset the mass attenuation coefficients of the element to the EPDL97 values
    element = &(this->elementList[this->elementDict[name]]);
    atomicNumber = this->elementList[this->elementDict[name]].getAtomicNumber();
    massAttenuationCoefficients = epdl97.getMassAttenuationCoefficients(atomicNumber);

    (*element).setMassAttenuationCoefficients(massAttenuationCoefficients["energy"],          \
                                       massAttenuationCoefficients["photoelectric"],          \
                                       massAttenuationCoefficients["coherent"],               \
                                       massAttenuationCoefficients["compton"],                \
                                       massAttenuationCoefficients["pair"]);

    for (j = 0; j < 10; j++)
    {
        shell = shellList[j];
        (*element).setPartialPhotoelectricMassAttenuationCoefficients(shell,     \
                                                    massAttenuationCoefficients["energy"], \
                                                    massAttenuationCoefficients[shell]);
    }

    // now we calculate the mass attenuation coefficients via EPDL97 at the given energies
    // extended by the EPDL97 ones at the beginning
    // first calculate the appropriate energy grid

    j = 0;
    while(massAttenuationCoefficients["energy"][j] < energy[0])
    {
        j++;
    }
    if (j > 0)
    {
        j--;
    }

    // partial cross sections
    // we extract the edges from supplied mass attenuation coefficients
    // the result contains the indices of the first instance of the energy
    extractedEdgeEnergies = (*element).extractEdgeEnergiesFromMassAttenuationCoefficients(energy, \
                                                                                          photoelectric);

    shellIndexList.clear();
    for (shellIndex = 0; shellIndex < 10; shellIndex++)
    {
        shell = shellList[shellIndex];
        if (extractedEdgeEnergies.find(shell) != extractedEdgeEnergies.end())
        {
            // we have the index prior to starting to excite that shell but that index corresponds to
            // the old, non extended, energy grid
            idx = extractedEdgeEnergies[shell].second;
            shellIndexList.push_back(idx);
        }
    }

    newEnergyGrid.resize(j + energy.size());
    for (i = 0; i < j; i++)
    {
        newEnergyGrid[i] = massAttenuationCoefficients["energy"][i];
    }

    for (i = j; i < newEnergyGrid.size() ; i++)
    {
        newEnergyGrid[i] = energy[i - j];
        // Check compatibility
        for (shellIndex = 0; shellIndex < 10; shellIndex++)
        {
            shell = shellList[shellIndex];
            if (extractedEdgeEnergies.find(shell) != extractedEdgeEnergies.end())
            {
                if ((i - j - 1) == static_cast<unsigned int>(extractedEdgeEnergies[shell].second))
                {
                    // We are just above an edge according to the input mass attenuation coefficients
                    // We have to check if that energy is compatible with the binding energy.
                    // If the energy was below the binding energy of that shell, we have to make sure the point is considered
                    // (Problem with excitation of Pb at 14 keV when using XCOM cross sections)
                    std::map<std::string, double> tmpMap = this->getBindingEnergies(name);
                    if (extractedEdgeEnergies[shell].first < tmpMap[shell])
                    {
                        // I should interpolate in the supplied grid to get the proper value
                        newEnergyGrid[i] = tmpMap[shell];
                        continue;
                    }
                }
            }
        }
        if (j > 0)
        {
            if (newEnergyGrid[i] < newEnergyGrid[i - 1])
            {
                newEnergyGrid[i] = newEnergyGrid[i - 1] + 0.000001;
            }
        }
    }

    // we have got the extended energy grid
    // recalculate the cross sections on this energy grid

    for (it = massAttenuationCoefficients.begin(); it != massAttenuationCoefficients.end(); ++it)
    {
        massAttenuationCoefficients[it->first].clear();
    }
    massAttenuationCoefficients.clear();
    massAttenuationCoefficients = (*element).getMassAttenuationCoefficients(newEnergyGrid);

    // insert the supplied mass attenuation coefficients
    if (massAttenuationCoefficients["energy"].size() != newEnergyGrid.size())
    {
        msg = "Error setting mass attenuation coefficients  of element " + name;
        throw std::runtime_error(msg);
    }

/*
    // test the energies
    for (i = j; i < newEnergyGrid.size(); i++)
    {
        if (fabs(massAttenuationCoefficients["energy"][i] - energy[i - j]) > 0.010)
        {
            std::cout << "Inconsistent energy grid for element " << name << std::endl;
            std::cout << " Energy = " << massAttenuationCoefficients["energy"][i] << std::endl;
            std::cout << " Input e = " << energy[i - j] << std::endl;
            msg = "Inconsistent energy grid for element " + name;
            //throw std::runtime_error(msg);
        }
    }
*/

    // coherent
    for (i = j; i < newEnergyGrid.size(); i++)
    {
        massAttenuationCoefficients["coherent"][i] = coherent[i - j];
    }

    // compton
    for (i = j; i < newEnergyGrid.size(); i++)
    {
        massAttenuationCoefficients["compton"][i] = compton[i - j];
    }

    // pair
    for (i = j; i < newEnergyGrid.size(); i++)
    {
        massAttenuationCoefficients["pair"][i] = pair[i - j];
    }

    // photoelectric
    // we keep a copy of the old photoelectric
    oldTotalPhotoelectric.resize(newEnergyGrid.size());
    for (i = 0; i < newEnergyGrid.size(); i++)
    {
        oldTotalPhotoelectric[i] = massAttenuationCoefficients["photoelectric"][i];
    }

    for (i = j; i < newEnergyGrid.size(); i++)
    {
        massAttenuationCoefficients["photoelectric"][i] = photoelectric[i - j];
    }


    // partial cross sections
    // we extract the edges from supplied mass attenuation coefficients
    // the result contains the indices of the first instance of the energy
    extractedEdgeEnergies = (*element).extractEdgeEnergiesFromMassAttenuationCoefficients(energy, \
                                                                                          photoelectric);

    if (extractedEdgeEnergies.size() > 0)
    {
        for (shellIndex = 0; shellIndex < 10; shellIndex++)
        {
            shell = shellList[shellIndex];
            if (extractedEdgeEnergies.find(shell) != extractedEdgeEnergies.end())
            {
                // we have the index prior to starting to excite that shell but that index corresponds to
                // the old, non extended, energy grid
                idx = extractedEdgeEnergies[shell].second;
                for (i = 0; i < (j + idx + 1); i++)
                {
                    massAttenuationCoefficients[shell][i] = 0.0;
                }
                for (i = (j + idx + 1); i < newEnergyGrid.size(); i++)
                {
                    tmpDouble = oldTotalPhotoelectric[i];
                    if ( tmpDouble > 0)
                    {
                        tmpDouble = massAttenuationCoefficients[shell][i] / tmpDouble;
                        massAttenuationCoefficients[shell][i] = photoelectric[i - j] * tmpDouble;
                    }
                    else
                    {
                        massAttenuationCoefficients[shell][i] = 0.0;
                    }
                }
            }
        }
    }
    else
    {
        // blindly use EPDL97
        for (shellIndex = 0; shellIndex < 10; shellIndex++)
        {
            shell = shellList[shellIndex];
            for (i = j; i < newEnergyGrid.size(); i++)
            {
                tmpDouble = oldTotalPhotoelectric[i];
                if ( tmpDouble > 0)
                {
                    tmpDouble = massAttenuationCoefficients[shell][i] / tmpDouble;
                    massAttenuationCoefficients[shell][i] = photoelectric[i - j] * tmpDouble;
                }
                else
                {
                    massAttenuationCoefficients[shell][i] = 0.0;
                }
            }
        }
    }

    // finally we are ready to update the mass attenuation coefficients
    (*element).setMassAttenuationCoefficients(massAttenuationCoefficients["energy"],\
                                       massAttenuationCoefficients["photoelectric"],   \
                                       massAttenuationCoefficients["coherent"],        \
                                       massAttenuationCoefficients["compton"],         \
                                       massAttenuationCoefficients["pair"]);

    for (shellIndex = 0; shellIndex < 10; shellIndex++)
    {
        shell = shellList[shellIndex];
        (*element).setPartialPhotoelectricMassAttenuationCoefficients(shell,\
                                            massAttenuationCoefficients["energy"], \
                                            massAttenuationCoefficients[shell]);
    }
}


std::map<std::string, std::vector<double> > Elements::getMassAttenuationCoefficients( \
                                                                const std::string & name) const
{
    std::string msg;
    std::map<std::string, int>::const_iterator it;

    it = this->elementDict.find(name);
    if (it == this->elementDict.end())
    {
        msg = "Name " + name + " not among defined elements";
        throw std::invalid_argument(msg);
    }
    return this->elementList[it->second].getMassAttenuationCoefficients();

}


std::map<std::string, double> Elements::getMassAttenuationCoefficients(std::string name, double energy) const
{
    std::string msg;
    std::map<std::string, double> composition;

    if (this->elementDict.find(name) == this->elementDict.end())
    {
        // composition = this->getCompositionFromFormula(name);
        composition = this->getComposition(name);
        if (composition.size() < 1)
        {
            msg = "Name " + name + " not accepted as element, material or chemical formula";
            throw std::invalid_argument(msg);
        }
        return this->getMassAttenuationCoefficients(composition, energy);
    }
    else
    {
        //elementList.getElement
        const Element & elementReference = this->getElement(name);
        return elementReference.getMassAttenuationCoefficients(energy);
    }
}

std::map<std::string, std::vector<double> > Elements::getMassAttenuationCoefficients( \
                                                        const std::string & name, \
                                                        const std::vector<double> & energy) const
{
    std::string msg;
    std::map<std::string, double> composition;
    std::map<std::string, int>::const_iterator c_it;

    c_it = this->elementDict.find(name);
    if (c_it == this->elementDict.end())
    {
        //composition = this->getCompositionFromFormula(name);
        composition = this->getComposition(name);
        if (composition.size() < 1)
        {
            msg = "Name " + name + " not accepted as element, material or chemical formula";
            throw std::invalid_argument(msg);
        }
        return this->getMassAttenuationCoefficients(composition, energy);
    }
    else
    {
        const Element & elementReference = this->getElement(name);
        return elementReference.getMassAttenuationCoefficients(energy);
    }
}

std::map<std::string, double> Elements::getMassAttenuationCoefficients(\
                                                std::map<std::string, double> inputFormulaDict,\
                                                double inputEnergy) const
{
    std::vector<double> energy;
    std::map<std::string, std::vector<double> >tmpResult;
    std::map<std::string, double> result;

    energy.resize(1);
    energy[0] = inputEnergy;
    tmpResult = this->getMassAttenuationCoefficients(inputFormulaDict, energy);

    result["energy"] = tmpResult["energy"][0];
    result["coherent"] = tmpResult["coherent"][0];
    result["compton"] = tmpResult["compton"][0];
    result["pair"] = tmpResult["pair"][0];
    result["photoelectric"] = tmpResult["photoelectric"][0];
    result["total"] = tmpResult["total"][0];

    return result;
}

std::map<std::string, std::vector<double> > Elements::getMassAttenuationCoefficients(\
                                                std::map<std::string, double> inputFormulaDict,\
                                                                std::vector<double> energy) const
{
    std::string element, msg, name;
    double total, massFraction;
    std::map<std::string, double>::const_iterator c_it;
    std::map<std::string, double> composition;
    std::map<std::string, double> tmpResult;
    std::map<std::string, std::vector<double> > result;
    std::map<std::string, std::vector<double> >::size_type n;
    std::map<std::string, double> elementsDict;
    std::map<std::string, double>::iterator it;
    std::map<std::string , int>::const_iterator mapIterator;

    total = 0.0;
    for (c_it = inputFormulaDict.begin(); c_it != inputFormulaDict.end(); ++c_it)
    {
        massFraction = c_it->second;
        if (massFraction < 0.0)
        {
            msg = "Name " + c_it->first + " has a negative mass fraction!!!";
            throw std::invalid_argument(msg);
        }
        // we may have received formulas ...
        name = c_it->first;
        //composition = this->getCompositionFromFormula(name);
        composition = this->getComposition(name);
        if (composition.size() < 1)
        {
            msg = "Name " + c_it->first + " not understood";
            std::cout << msg << std::endl;
            throw std::invalid_argument(msg);
        }
        for(it = composition.begin(); it != composition.end(); ++it)
        {
            composition[it->first] *= massFraction;
            if (elementsDict.find(it->first) == elementsDict.end())
            {
                elementsDict[it->first] = 0.0;
            }
            elementsDict[it->first] += composition[it->first];
        }
        total += massFraction;
    }

    if (total <= 0.0)
    {
        msg = "Sum of mass fractions is less or equal to 0";
        throw std::invalid_argument(msg);
    }

    result["energy"].resize(energy.size());
    result["coherent"].resize(energy.size());
    result["compton"].resize(energy.size());
    result["pair"].resize(energy.size());
    result["photoelectric"].resize(energy.size());
    result["total"].resize(energy.size());

    for (n = 0; n < result["energy"].size(); n++)
    {
        result["energy"][n] = energy[n];
        result["coherent"][n] = 0.0;
        result["compton"][n] = 0.0;
        result["pair"][n] = 0.0;
        result["photoelectric"][n] = 0.0;
        for (c_it = elementsDict.begin(); c_it != elementsDict.end(); ++c_it)
        {

            massFraction = c_it->second / total;
            element = c_it->first;
            mapIterator = this->elementDict.find(element);
            tmpResult = this->elementList[mapIterator->second].getMassAttenuationCoefficients(energy[n]);
            result["coherent"][n] += tmpResult["coherent"] * massFraction;
            result["compton"][n] += tmpResult["compton"] * massFraction;
            result["pair"][n] += tmpResult["pair"] * massFraction;
            result["photoelectric"][n] += tmpResult["photoelectric"] * massFraction;
        }

        result["total"] [n] = (result["coherent"][n] + result["compton"][n]) + \
                               result["pair"][n] + result["photoelectric"][n];
    }
    return result;
}


// Materials handling
Material Elements::createMaterial(const std::string & name, const double & density, \
                              const double & thickness, const std::string & comment)
{
    std::string msg;
    Material material;
    std::map<std::string, double> composition;

    material.initialize(name, density, thickness, comment);
    composition = this->getCompositionFromFormula(name);
    if (composition.size() > 0)
    {
        material.setComposition(composition);
    }
    return material;
}

void Elements::addMaterial(const std::string & name, const double & density, \
                              const double & thickness, const std::string & comment, \
                              const int & errorOnReplace)
{
    std::string msg;
    Material material;
    std::map<std::string, double> composition;

    if (this->getMaterialIndexFromName(name) < this->materialList.size())
    {
        if (errorOnReplace)
        {
            msg = "Elements::createMaterial. Already existing material: " +  name;
            throw std::invalid_argument(msg);
        }
        else
        {
            this->removeMaterial(name);
        }
    }
    material.initialize(name, density, thickness, comment);
    this->materialList.push_back(material);

    // Try to set the composition from the name
    composition = this->getCompositionFromFormula(name);
    if (composition.size() > 0)
    {
        this->setMaterialComposition(name, composition);
    }
}

void Elements::setMaterialComposition(const std::string & materialName, \
                                      const std::vector<std::string> & names,\
                                      const std::vector<double> & amounts)
{
    std::string msg;
    std::vector<Material>::size_type i;

    i = this->getMaterialIndexFromName(materialName);
    if (i == this->materialList.size())
    {
        msg = "Elements::setMaterialComposition. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    this->materialList[i].setComposition(names, amounts);
}

void Elements::setMaterialComposition(const std::string & materialName, \
                                const std::map<std::string, double> & composition)
{
    std::string msg;
    std::vector<Material>::size_type i;

    i = this->getMaterialIndexFromName(materialName);
    if (i >= this->materialList.size())
    {
        msg = "Elements::setMaterialComposition. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    this->materialList[i].setComposition(composition);
}

const Material & Elements::getMaterial(const std::string & materialName)
{
    std::string msg;
    std::vector<Material>::size_type i;

    i = this->getMaterialIndexFromName(materialName);
    if (i >= this->materialList.size())
    {
        msg = "Elements::getMaterial. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    return this->materialList[i];
}

std::vector<std::string> Elements::getMaterialNames()
{
    std::vector<Material>::size_type i;
    std::vector<std::string> result;

    result.resize(this->materialList.size());
    for (i = 0; i < this->materialList.size(); i++)
    {
        result[i] = this->materialList[i].getName();
    }

    return result;
}

Material Elements::getMaterialCopy(const std::string & materialName)
{
    std::string msg;
    std::vector<Material>::size_type i;

    i = this->getMaterialIndexFromName(materialName);
    if (i >= this->materialList.size())
    {
        msg = "Elements::getMaterial. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    return this->materialList[i];
}


void Elements::addMaterial(const Material & material, const int & errorOnReplace)
{
    std::string msg;
    std::string materialName;
    std::vector<Material>::size_type i;


    materialName = material.getName();

    i = this->getMaterialIndexFromName(materialName);
    if (i < this->materialList.size())
    {
        if (errorOnReplace)
        {
            msg = "Elements::addMaterial. Already existing material: " +  materialName;
            throw std::invalid_argument(msg);
        }
        else
        {
            // try to use the same position
            this->materialList[i] = material;
        }
    }
    else
    {
        this->materialList.push_back(material);
    }
    // TODO: Make sure the material can be interpreted in terms of the supplied composition
    // because the composition can include other materials.
    // If made that way, the internal list of materials will be "clean" of references from
    // one material to another. Is it needed?
}

void Elements::removeMaterials()
{
    this->materialList.clear();
}

std::map<std::string, double> Elements::getComposition(const std::string & name) const
{
    std::string msg;
    std::map<std::string, double> result;
    std::map<std::string, double> tmpResult;
    std::map<std::string, double> composition;
    std::map<std::string, double>::const_iterator c_it, c_it2;
    std::map<std::string , int>::const_iterator matIterator;
    std::vector<Material>::size_type i;
    double total;

    // check if name is a valid element or formula
    result = this->getCompositionFromFormula(name);
    if (result.size() > 0)
    {
        return result;
    }

    // check if it is a material
    i = this->getMaterialIndexFromName(name);
    if (i == this->materialList.size())
    {
        // result at this point must be empty, we can send it back.
        return result;
    }

    // make sure we give back the elemental material composition
    tmpResult = this->materialList[i].getComposition();
    if (tmpResult.size() < 1)
    {
        // throw an exception because we have an undefined material in the list of materials
        msg = "Material " + name + " with empty or non-valid composition";
        throw std::invalid_argument(msg);
    }

    // make sure (again ?) that material is normalized
    total = 0;
    for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
    {
        total += c_it->second;
    }
    for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
    {
        tmpResult[c_it->first] /= total;
    }

    // now, a material can be made of formulas, elements or other materials.
    for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
    {
        composition = this->getComposition(c_it->first);
        if (composition.size() < 1)
        {
            return composition;
        }
        for (c_it2 = composition.begin(); c_it2 != composition.end(); ++c_it2)
        {
            if (result.find(c_it2->first) == result.end())
            {
                result[c_it2->first] = 0.0;
            }
            result[c_it2->first] += composition[c_it2->first] * tmpResult[c_it->first];
        }
    }
    return result;
}

std::map<std::string, double> Elements::getComposition(const std::string & name, \
                                                 const std::vector<Material> & additionalMaterials) const
{
    std::vector<Material>::size_type i;
    std::map<std::string, double> result;
    std::map<std::string, double> tmpResult;
    std::map<std::string, double> composition;
    std::map<std::string, double>::const_iterator c_it, c_it2;
    double total;

    // check if name is a valid element or formula
    result = this->getCompositionFromFormula(name);
    if (result.size() > 0)
    {
        return result;
    }

    // We are dealing with a material
    if (additionalMaterials.size() == 0)
    {
        // The additional materials are not going to be needed
        return this->getComposition(name);
    }

    // check if it is a material of the supplied list
    for (i = 0; i < additionalMaterials.size(); i++)
    {
        if(additionalMaterials[i].getName() == name)
        {
            break;
        }
    }
    if (i == additionalMaterials.size())
    {
        // The additional materials not needed in this iteration
        tmpResult =  this->getComposition(name);
    }
    else
    {
        tmpResult =  additionalMaterials[i].getComposition();
    }
    if (tmpResult.size() < 1)
    {
        // throw an exception because we have an undefined material in the list of materials
        std::string msg;
        msg = "Material " + name + " with empty or non-valid composition";
        throw std::invalid_argument(msg);
    }
    // We are dealing with a material that can be defined (or redefined) using the
    // supplied list and a material can be made of formulas, elements or other materials.
    for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
    {
        composition = this->getComposition(c_it->first, additionalMaterials);
        if (composition.size() < 1)
        {
            return composition;
        }
        for (c_it2 = composition.begin(); c_it2 != composition.end(); ++c_it2)
        {
            if (result.find(c_it2->first) == result.end())
            {
                result[c_it2->first] = 0.0;
            }
            result[c_it2->first] += composition[c_it2->first] * tmpResult[c_it->first];
        }
    }
    // make sure (again ?) that material is normalized
    total = 0;
    for (c_it = result.begin(); c_it != result.end(); ++c_it)
    {
        total += c_it->second;
    }
    std::cout << "total = " << total << std::endl;
    for (c_it = result.begin(); c_it != result.end(); ++c_it)
    {
        result[c_it->first] /= total;
    }
    return result;
}

std::map<std::string, double> Elements::getCompositionFromFormula(const std::string & formula) const
{
    std::map<std::string, double> parsedFormula;
    std::map<std::string, double>::iterator it;
    std::string name, msg;
    double total;
    std::map<std::string , int>::const_iterator matIterator;
    // TODO: Still to multiply by Atomic Weight!!!!

    parsedFormula = this->parseFormula(formula);
    if (parsedFormula.size() < 1)
    {
        // Formula not understood
        return parsedFormula;
    }
    total = 0.0;
    for (it = parsedFormula.begin(); it != parsedFormula.end(); ++it)
    {
        name = it->first;
        matIterator = this->elementDict.find(name);
        if (matIterator == this->elementDict.end())
        {
            // msg = "Name " + name + " not among defined elements";
            // std::cout << msg << std::endl;
            parsedFormula.clear();
            return parsedFormula;
        }
        // multiply the number of atoms by the atomic weight
        parsedFormula[name] *= this->elementList[matIterator->second].getAtomicMass();
        total += parsedFormula[name];
    }
    for (it = parsedFormula.begin(); it != parsedFormula.end(); ++it)
    {
        parsedFormula[it->first] /= total;
    }
    return parsedFormula;
}

std::vector<std::map<std::string, std::map<std::string, double> > >Elements::getExcitationFactors( \
                            const std::string & element,
                            const std::vector<double> & energies,
                            const std::vector<double> & weights) const
{
    const Element & elementObject = this->getElement(element);
    std::map<std::string, std::map<std::string, double> > result;

    // for the time being no other way to produce vacancies than photoelectric effect
    return elementObject.getPhotoelectricExcitationFactors(energies, weights);
}

std::map<std::string, std::map<std::string, double> > Elements::getExcitationFactors( \
                            const std::string & element,
                            const double & energy,
                            const double & weight) const
{
    std::vector<double> energies;
    std::vector<double> weights;

    energies.push_back(energy);
    weights.push_back(weight);
    return this->getExcitationFactors(element, energies, weights)[0];
}

std::map<std::string, double> Elements::parseFormula(const std::string & formula) const
{
    std::map<std::string, double> composition;
    std::map<std::string, double> tmpComposition;
    std::map<std::string, double>::iterator it;
    std::string::size_type i, p1, p2, length;
    std::string::const_iterator c_it;
    std::vector<int> openParenthesis;
    std::vector<int> closeParenthesis;
    std::vector<int>::size_type nParenthesis;
    std::string newFormula;
    bool parsingKey;
    std::string lastKey;
    std::string lastNumber;
    std::vector<std::string> keys;
    std::vector<double> numbers;
    bool spacesPresent;
    double factor;


    //std::cout << "Received formula = " << formula << std::endl;
    composition.clear();

    if (formula.size() < 1)
        return composition;

    if (islower(formula[0]) || isdigit(formula[0])) {
        return composition;
    }

    // look for parenthesis
    parsingKey = false;
    spacesPresent = false;
    for(i=0; i < formula.size(); i++)
    {
        if (formula[i] =='(')
        {
            openParenthesis.push_back(i);
        }
        if (formula[i] ==')')
        {
            closeParenthesis.push_back(i);
        }
        if (formula[i] == ' ')
        {
            spacesPresent = true;
        }
    }
    // if number of opened and closed do not match
    // return an empty map (or should I throw an exception
    // for just this case?
    if (openParenthesis.size() != closeParenthesis.size())
    {
        //std::cout << "RECEIVED Formula NOT passed parenthesis check" << std::endl;
        return composition;
    }

    if (spacesPresent)
    {
        //std::cout << "RECEIVED Formula NOT passed spaces check" << std::endl;
        return composition;
    }

    nParenthesis = openParenthesis.size();
    if ( nParenthesis > 0)
    {
        if (nParenthesis > 1)
        {
            p1 = 0;
            for (i = 0; i < nParenthesis; ++i)
            {
                if (openParenthesis[i] < closeParenthesis[0])
                    p1 = openParenthesis[i];
            }
        }
        else
        {
            p1 = openParenthesis[0];
        }
        p2 = closeParenthesis[0];
        length = p2 - p1 - 1;
        if(length < 1)
        {
            // empty substring
            return composition;
        }
        tmpComposition = this->parseFormula(formula.substr(p1 + 1, length));
        if (tmpComposition.size() < 1)
        {
            return tmpComposition;
        }
        // get the numbers assoiated to the parenthesis
        if ((p2+1) == formula.size())
        {
            // the number associated  to that parenthesis is 1
            factor = 1.0;
        }
        else
        {
            i = p2 + 1;
            while( i < formula.size())
            {
                if (isupper(formula[i]) || islower(formula[i]) || (formula[i] == '(') || (formula[i] == ')'))
                {
                    break;
                }
                i += 1;
            }
            if (!this->stringToDouble(formula.substr(p2+1, i - (p2+1)), factor))
            {
                return composition;
            };
        }
        /* better generate a new formula than to try to deal with everything here
        for(it = tmpComposition.begin(); it != tmpComposition.end(); ++it)
        {
            tmpComposition[it->first] *= factor;
        }
        */
        if (p1 == 0)
        {
            newFormula = "";
        }
        else
        {
            newFormula = formula.substr(0, p1);
        }
        for(it = tmpComposition.begin(); it != tmpComposition.end(); ++it)
        {
            newFormula += it->first + this->toString(it->second * factor);
        }
        if (i < (formula.size()-1))
        {
            newFormula += formula.substr(i, formula.size() - i);
        }
        return this->parseFormula(newFormula);
    }
    else
    {
        // when we arrive here we have a string without parenthesis
        if (!isupper(formula[0]))
        {
            // invalid first character;
            return composition;
        }
        i = 0;
        lastKey = "";
        while (i < formula.size())
        {
            if (isupper(formula[i]))
            {
                if(lastKey.size() > 0)
                {
                    keys.push_back(lastKey);
                    if (lastNumber.size() > 0)
                    {
                        if (this->stringToDouble(lastNumber, factor))
                            numbers.push_back(factor);
                        else
                            return composition;
                    }
                    else
                    {
                        numbers.push_back(1.0);
                    }
                }
                lastNumber = "";
                parsingKey = true;
                lastKey = formula.substr(i, 1);
            }
            else
            {
                if (islower(formula[i]))
                {
                    lastKey += formula.substr(i, 1);
                }
                else
                {
                    // parsing number
                    parsingKey = false;
                    lastNumber += formula.substr(i, 1);
                }
            }
            i += 1;
        }
        if (i == formula.size())
        {
            if (parsingKey)
            {
                keys.push_back(lastKey);
                numbers.push_back(1.0);
            }
            else
            {
                keys.push_back(lastKey);
                if (this->stringToDouble(lastNumber, factor))
                    numbers.push_back(factor);
                else
                    return composition;
            }
        }
        if (keys.size() != numbers.size())
        {
            return composition;
        }
        for (i = 0; i < keys.size(); ++ i)
        {
            if (composition.find(keys[i]) == composition.end())
            {
                composition[keys[i]] = numbers[i];
            }
            else
            {
                composition[keys[i]] += numbers[i];
            }
        }
    }
    return composition;
}

const std::string & Elements::getShellConstantsFile(const std::string & mainShellName) const
{
    std::map<std::string, std::string>::const_iterator c_it;

    c_it = this->shellConstantsFile.find(mainShellName);
    if (c_it == this->shellConstantsFile.end())
    {
        throw std::invalid_argument("Invalid main shell. It should be K, L or M");
    }
    return c_it->second;
}

const std::string & Elements::getShellRadiativeTransitionsFile( \
                                                const std::string & mainShellName) const
{
    std::map<std::string, std::string>::const_iterator c_it;

    c_it = this->shellRadiativeTransitionsFile.find(mainShellName);
    if (c_it == this->shellRadiativeTransitionsFile.end())
    {
        throw std::invalid_argument("Invalid main shell. It should be K, L or M");
    }
    return c_it->second;
}

const std::string & Elements::getShellNonradiativeTransitionsFile( \
                                                const std::string & mainShellName) const
{
    std::map<std::string, std::string>::const_iterator c_it;

    c_it = this->shellNonradiativeTransitionsFile.find(mainShellName);
    if (c_it == this->shellNonradiativeTransitionsFile.end())
    {
        throw std::invalid_argument("Invalid main shell. It should be K, L or M");
    }
    return c_it->second;
}

std::vector<std::pair<std::string, double> > Elements::getPeakFamilies( \
                            const std::vector<std::string> & elementList, const double & energy) const
{
    std::map<std::string, double>::const_iterator c_it;
    std::vector<std::string>::size_type i, j;
    std::vector<std::string> shells;
    std::vector<std::pair<std::string, double> >result;

    result.clear();
    for (i = 0; i < elementList.size(); i++)
    {
        shells = getElement(elementList[i]).getExcitedShells(energy);
        if (shells.size())
        {
            const std::map<std::string, double> & bindingEnergies = \
                        getElement(elementList[i]).getBindingEnergies();
            for (j = 0; j < shells.size(); j++)
            {
                c_it = bindingEnergies.find(shells[j]);
                if ((shells[j][0] == 'K') || (shells[j][0] == 'L') || (shells[j][0] == 'M'))
                {
                    if (getElement(elementList[i]).getShellConstants(shells[j])["omega"] > 0.0)
                    {
                        result.push_back(std::make_pair(elementList[i] + " " + shells[j], c_it->second));
                    }
                }
            }
        }
    }
    if(result.size())
        std::sort(result.begin(), result.end(), sortVectorOfExcited());

    // now fill an easier to wrap vector
    return result;
};


const std::map<std::string, double> & Elements::getBindingEnergies(const std::string & element) const
{
    return getElement(element).getBindingEnergies();
}

std::vector<std::pair<std::string, double> > Elements::getPeakFamilies( \
                                            const std::string & name, const double & energy) const
{
    std::map<std::string, double> composition;
    std::map<std::string, double>::const_iterator c_it;
    std::vector<std::string> names;
    std::string msg;

    composition = this->getComposition(name);
    if (composition.size() < 1)
    {
        msg = "Name " + name + " not accepted as element, material or chemical formula";
        throw std::invalid_argument(msg);
    }
    names.clear();
    for (c_it = composition.begin(); c_it != composition.end(); ++c_it)
    {
        names.push_back(c_it->first);
    }
    return this->getPeakFamilies(names, energy);
}


bool Elements::stringToDouble(const std::string& str, double& number)
{

    std::istringstream i(str);

    if (!(i >> number))
    {
        // Number conversion failed
        return false;
    }
    return true;
}

std::string Elements::toString(const double& number)
{
    std::ostringstream ss;
    ss << number;
    std::string s(ss.str());
    return s;
}

std::map<std::string, double> Elements::getEmittedXRayLines(const std::string & elementName, const double & energy) const
{
    return this->getElement(elementName).getEmittedXRayLines(energy);
}

const std::map<std::string, double> & Elements::getNonradiativeTransitions(const std::string & elementName, \
                                                                     const std::string & subshell) const
{
    return this->getElement(elementName).getNonradiativeTransitions(subshell);
}

const std::map<std::string, double> & Elements::getRadiativeTransitions(const std::string & elementName, \
                                                                  const std::string & subshell) const
{
    return this->getElement(elementName).getRadiativeTransitions(subshell);
}

std::map<std::string, double> Elements::getShellConstants(const std::string & elementName, \
                                                    const std::string & subshell) const
{
    return this->getElement(elementName).getShellConstants(subshell);
}



std::map<std::string,std::map<std::string, double> > Elements::getEscape( \
                                        const std::map<std::string, double> & composition,
                                        const double & energy, \
                                        const double & energyThreshold, \
                                        const double & intensityThreshold, \
                                        const int & nThreshold , \
                                        const double & alphaIn , \
                                        const double & thickness) const
{
    std::map<std::string,  std::map<std::string, double> > result;
    std::map<std::string, std::map<std::string, double> > excitationFactors;
    std::map<std::string, double>::const_iterator c_it;
    std::string element;
    double massFraction;
    double muIncident;
    double muFluorescence;
    double sinAlphaIn;
    double tmpDouble;
    std::string tmpString;
    double intrinsicEfficiency;

    if (this->isEscapeCacheCompatible(composition,energyThreshold,intensityThreshold, \
                                        nThreshold , alphaIn, thickness))
    {
        std::map< double, std::map<std::string,std::map<std::string, double> > >::const_iterator it;
        it = this->detectorEscapeCache.find(energy);
        if ( it != this->detectorEscapeCache.end())
        {
            // std::cout << "USING CACHE" <<  energy << std::endl;
            return it->second;
        }
    }

    if (alphaIn == 90.)
        sinAlphaIn = 1.0;
    else
    {
        sinAlphaIn = std::sin(alphaIn * (3.141592653589793/180.));
        if (sinAlphaIn < 0.0)
        {
            sinAlphaIn = - sinAlphaIn;
        }
    }
    muIncident = this->getMassAttenuationCoefficients(composition, energy)["total"];
    result.clear();

    if (thickness > 0.0)
    {
        intrinsicEfficiency = (1.0 - std::exp(-muIncident * thickness / sinAlphaIn));
    }
    else
    {
        // this is just to get rid of a not very clever non-initialized warning because
        // if thickness is not bigger than 0.0 intrinsicEfficiency will not be used in
        // any case.
        intrinsicEfficiency = 1.0;
    }
    for (c_it = composition.begin(); c_it != composition.end(); c_it++)
    {
        element = c_it->first;
        massFraction = c_it->second;
        excitationFactors = this->getExcitationFactors(element, energy, massFraction);
        std::map<std::string, std::map<std::string, double> >::const_iterator it;
        std::map<std::string, double>::const_iterator mapIt;
        // factor uses mu["shell photoelectric"] / mu["total photoelectric"]
        // double factor;
        // rate uses mu["shell photoelectric"]
        // It is the product factor * mu["total photoelectric"]
        double rate;
        double fluorescentEnergy;
        for ( it = excitationFactors.begin(); it != excitationFactors.end(); ++it)
        {
            // mapIt = it->second.find("factor");
            // factor = mapIt->second;
            mapIt = it->second.find("rate");
            rate = mapIt->second;
            mapIt = it->second.find("energy");
            fluorescentEnergy = mapIt->second;
            muFluorescence = this->getMassAttenuationCoefficients(composition, fluorescentEnergy)["total"];
            tmpDouble = sinAlphaIn * (muFluorescence / muIncident);
            tmpString = element + "_" + it->first + "esc";
            rate *= (0.5 /  muIncident) * ( 1.0 - tmpDouble * std::log( 1 + 1.0 / tmpDouble));
            if (rate > intensityThreshold)
            {
                if (thickness > 0.0)
                {
                    // This is to give the escape peak rate per detected photon
                    // and not per incident photon.
                    // It is not a correct approximation, but this avoids the calculation
                    // of exponential integral functions.
                    rate /= intrinsicEfficiency;
                }
                result[tmpString] ["rate"] = rate;
                result[tmpString] ["energy"] = energy - fluorescentEnergy;
            }
        }
    }
    return result;
}

const std::vector<Material>::size_type Elements::getMaterialIndexFromName(const std::string & name) const
{
    std::vector<Material>::size_type i;

    this->materialList.begin();
    for (i = 0; i < this->materialList.size(); i++)
    {
        if(this->materialList[i].getName() == name)
        {
            return i;
        }
    }
    return this->materialList.size();
}


void Elements::removeMaterial(const std::string & name)
{
    std::string msg;
    std::vector<Material>::size_type i;
    i = this->getMaterialIndexFromName(name);
    if ( i >= this->materialList.size())
    {
        msg = "Elements::setMaterialComposition. Non existing material: " +  name;
        throw std::invalid_argument(msg);
    }
    this->materialList.erase(this->materialList.begin() + i);
}

void Elements::setElementCascadeCacheEnabled(const std::string & elementName, const int & flag)
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        this->elementList[i].setCascadeCacheEnabled(flag);
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}
int Elements::isElementCascadeCacheFilled(const std::string & elementName) const
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].isCascadeCacheFilled();
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}
void Elements::fillElementCascadeCache(const std::string & elementName)
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].fillCascadeCache();
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}
void Elements::emptyElementCascadeCache(const std::string & elementName)
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].emptyCascadeCache();
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}

void Elements::fillCache(const std::string & elementName, const std::vector< double> & energy)
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].fillCache(energy);
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}

void Elements::updateCache(const std::string & elementName, const std::vector< double> & energy)
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].updateCache(energy);
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}


void Elements::setCacheEnabled(const std::string & elementName, const int & flag)
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].setCacheEnabled(flag);
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}

void Elements::clearCache(const std::string & elementName)
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].clearCache();
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}

void Elements::clearEscapeCache(void)
{
    this->detectorEscapeCache.clear();
}


bool Elements::isEscapeCacheCompatible(\
                                        const std::map<std::string, double> & composition,
                                        const double & energyThreshold, \
                                        const double & intensityThreshold, \
                                        const int & nThreshold , \
                                        const double & alphaIn , \
                                        const double & thickness) const
{
    // std::cout << 0 << std::endl;
    if (this->isEscapeCacheEnabled() == 0)
    {
        return false;
    }
    if (this->detectorEscapeCache.size() > 0)
    {
        // std::cout << 1 << std::endl;
        if (energyThreshold == this->detectorEnergyThresholdUsedInCache)
        {
            // std::cout << 2 << std::endl;
            if (intensityThreshold == this->detectorIntensityThresholdUsedInCache)
            {
                // std::cout << 3 << std::endl;
                // std::cout << " input " << nThreshold << " cached " << this->detectorNThresholdUsedInCache << std::endl;
                if(nThreshold == this->detectorNThresholdUsedInCache)
                {
                    // std::cout << 4 << std::endl;
                    if (alphaIn == this->detectorAlphaInUsedInCache)
                    {
                        //std::cout << 5 << std::endl;
                        //std::cout << " input " << thickness << " cached " << this->detectorThicknessUsedInCache << std::endl;
                        if (thickness == this->detectorThicknessUsedInCache)
                        {
                            //std::cout << 6 << std::endl;
                            // we have to compare the composition
                            if (composition.size() == this->detectorCompositionUsedInCache.size())
                            {
                                // std::cout << 7 << std::endl;
                                if (std::equal(composition.begin(), composition.end(), this->detectorCompositionUsedInCache.begin()))
                                {
                                    // std::cout << 8 << std::endl;
                                    return true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return false;
}

void Elements::updateEscapeCache(\
                                        const std::map<std::string, double> & composition,
                                        const std::vector<double> & energyList, \
                                        const double & energyThreshold, \
                                        const double & intensityThreshold, \
                                        const int & nThreshold , \
                                        const double & alphaIn , \
                                        const double & thickness)
{
    std::vector<double>::size_type i;
    double energy;

    if (this->isEscapeCacheEnabled() == 0)
    {
        std::cout << "WARNING: Filling escape cache when escape cache is disabled" << std::endl;
    }

    if (!this->isEscapeCacheCompatible(composition, energyThreshold, \
                                        intensityThreshold, nThreshold, alphaIn, thickness))
    {
        // not compatible, clear the detector cache
        this->clearEscapeCache();
    }
    for (i = 0; i < energyList.size(); ++i)
    {
        energy = energyList[i];
        if (this->detectorEscapeCache.find(energy) == this->detectorEscapeCache.end())
        {
            // std::cout << "filling energy " << energy << std::endl;
            this->detectorEscapeCache[energy] = this->getEscape(composition, energy, energyThreshold, \
                                                                intensityThreshold, nThreshold, alphaIn, thickness);
        }
        this->detectorCompositionUsedInCache = composition;
        this->detectorEnergyThresholdUsedInCache = energyThreshold;
        this->detectorIntensityThresholdUsedInCache = intensityThreshold;
        this->detectorNThresholdUsedInCache = nThreshold;
        // std::cout << "Storing cache " << nThreshold << std::endl;
        this->detectorAlphaInUsedInCache = alphaIn;
        this->detectorThicknessUsedInCache = thickness;
    }
}


const int Elements::isCacheEnabled(const std::string & elementName) const
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].isCacheEnabled();
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}

int Elements::getCacheSize(const std::string & elementName) const
{
    std::map<std::string, int>::const_iterator it;
    int i;
    if (this->isElementNameDefined(elementName))
    {
        it = this->elementDict.find(elementName);
        i = it->second;
        return this->elementList[i].getCacheSize();
    }
    else
        throw std::invalid_argument("Invalid element: " + elementName);
}


} // namespace fisx
