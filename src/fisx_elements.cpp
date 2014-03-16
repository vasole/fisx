#include <iostream>
#include <math.h>
#include <stdexcept>
#include <cctype>
#include <sstream>
#include "fisx_elements.h"

Elements::Elements(std::string epdl97Directory, std::string bindingEnergiesFileName, std::string crossSectionsFile)
{
    // this is to simplify an initialization equivalent to that of PyMca:
    // - use EPDL97 to calculate partial photoelectric cross sections
    // - use a different binding energies file
    // - use different mass atenuation coefficients (i.e. XCOM)
    this->initialize(epdl97Directory, bindingEnergiesFileName);
    this->setMassAttenuationCoefficientsFile(crossSectionsFile);
}

Elements::Elements(std::string epdl97Directory, std::string bindingEnergiesFileName)
{
    this->initialize(epdl97Directory, bindingEnergiesFileName);
}

Elements::Elements(std::string directoryName)
{
    // pure EPDL97 initialization
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

    // binding energies
    bindingEnergies = directoryName + joinSymbol + BINDING_ENERGIES;

    this->initialize(directoryName, "");
}


void Elements::initialize(std::string epdl97Directory, std::string bindingEnergiesFile)
{
#include "fisx_defaultelementsinfo.h"
    std::string shell;
    int i, j;
    std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};
    std::string symbol;
    int    atomicNumber;
    std::vector<double> energyGrid;
    std::map<std::string, std::vector<double> > massAttenuationCoefficients;
    std::string joinSymbol;
    std::string filename;

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
    this->setShellConstantsFile(filename, "K");
    filename = epdl97Directory + joinSymbol + L_SHELL_CONSTANTS_FILE;
    this->setShellConstantsFile(filename, "L");
    filename = epdl97Directory + joinSymbol + M_SHELL_CONSTANTS_FILE;
    this->setShellConstantsFile(filename, "M");
    return;
    filename = epdl97Directory + joinSymbol + K_SHELL_RADIATIVE_FILE;
    this->setShellRadiativeTransitionsFile(filename, "K");
    filename = epdl97Directory + joinSymbol + K_SHELL_NONRADIATIVE_FILE;
    this->setShellNonradiativeTransitionsFile(filename, "K");

    // Load the L shell
    filename = epdl97Directory + joinSymbol + L_SHELL_RADIATIVE_FILE;
    this->setShellRadiativeTransitionsFile(filename, "L");
    filename = epdl97Directory + joinSymbol + L_SHELL_NONRADIATIVE_FILE;
    this->setShellNonradiativeTransitionsFile(filename, "L");

    // Load the M shell
    filename = epdl97Directory + joinSymbol + M_SHELL_RADIATIVE_FILE;
    this->setShellRadiativeTransitionsFile(filename, "M");
    filename = epdl97Directory + joinSymbol + M_SHELL_NONRADIATIVE_FILE;
    this->setShellNonradiativeTransitionsFile(filename, "M");
}

// Element handling
bool Elements::isElementNameDefined(const std::string & elementName)
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

const Element & Elements::getElement(const std::string & elementName)
{
    if (this->isElementNameDefined(elementName))
        return this->elementList[this->elementDict[elementName]];
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

void Elements::addElement(Element & element)
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
void Elements::setShellConstantsFile(std::string fileName, std::string mainShellName)
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
}

void Elements::setShellRadiativeTransitionsFile(std::string fileName, std::string mainShellName)
{

}

void Elements::setShellNonradiativeTransitionsFile(std::string fileName, std::string mainShellName)
{

}

// Mass attenuation handling
void Elements::setMassAttenuationCoefficientsFile(std::string fileName)
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

void Elements::setMassAttenuationCoefficients(std::string name,
                                              std::vector<double> energy,
                                              std::vector<double> photoelectric,
                                              std::vector<double> coherent,
                                              std::vector<double> compton,
                                              std::vector<double> pair)
{
    std::map<std::string, std::vector<double> > massAttenuationCoefficients;
    std::map<std::string, std::vector<double> >::iterator it;
    std::map<std::string, std::pair<double, int> > extractedEdgeEnergies;
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

    newEnergyGrid.resize(j + energy.size());
    for (i = 0; i < j; i++)
    {
        newEnergyGrid[i] = massAttenuationCoefficients["energy"][i];
    }
    for (i=j; i < newEnergyGrid.size() ; i++)
    {
        newEnergyGrid[i] = energy[i - j];
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

    // test the energies
    for (i = j; i < newEnergyGrid.size(); i++)
    {
        if (fabs(massAttenuationCoefficients["energy"][i] - energy[i - j]) > 1.0e-10)
        {
            msg = "Inconsistent energy grid for element " + name;
            throw std::runtime_error(msg);
        }
    }

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


std::map<std::string, std::vector<double> > Elements::getMassAttenuationCoefficients(std::string name)
{
    std::string msg;
    if (this->elementDict.find(name) == this->elementDict.end())
    {
        msg = "Name " + name + " not among defined elements";
        throw std::invalid_argument(msg);
    }
    return this->elementList[this->elementDict[name]].getMassAttenuationCoefficients();

}


std::map<std::string, double> Elements::getMassAttenuationCoefficients(std::string name, double energy)
{
    std::string msg;
    std::map<std::string, double> composition;

    if (this->elementDict.find(name) == this->elementDict.end())
    {
        composition = this->getCompositionFromFormula(name);
        if (composition.size() < 1)
        {
            msg = "Name " + name + " not accepted as element or chemical formula";
            throw std::invalid_argument(msg);
        }
        return this->getMassAttenuationCoefficients(composition, energy);
    }
    else
    {
        return this->elementList[this->elementDict[name]].getMassAttenuationCoefficients(energy);
    }
}

std::map<std::string, std::vector<double> > Elements::getMassAttenuationCoefficients(std::string name, \
                                                                                     std::vector<double> energy)
{
    std::string msg;
    std::map<std::string, double> composition;

    if (this->elementDict.find(name) == this->elementDict.end())
    {
        composition = this->getCompositionFromFormula(name);
        if (composition.size() < 1)
        {
            msg = "Name " + name + " not accepted as element or chemical formula";
            throw std::invalid_argument(msg);
        }
        return this->getMassAttenuationCoefficients(composition, energy);
    }
    else
    {
        return this->elementList[this->elementDict[name]].getMassAttenuationCoefficients(energy);
    }
}

std::map<std::string, double> Elements::getMassAttenuationCoefficients(\
                                                std::map<std::string, double> inputFormulaDict,\
                                                double inputEnergy)
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
                                                                std::vector<double> energy)
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
        composition = this->getCompositionFromFormula(name);
        if (composition.size() < 1)
        {
            msg = "Name " + c_it->first + " not understood";
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
            tmpResult = this->elementList[this->elementDict[element]].getMassAttenuationCoefficients(energy[n]);
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
                              const double & thickness, const std::string & comment)
{
    std::string msg;
    Material material;
    std::map<std::string, double> composition;


    if (this->materialDict.find(name) != this->materialDict.end())
    {
        msg = "Elements::createMaterial. Already existing material: " +  name;
        throw std::invalid_argument(msg);
    }
    material.initialize(name, density, thickness, comment);
    this->materialList.resize(this->materialDict.size() + 1);
    this->materialList[this->materialDict.size()] = material;
    this->materialDict[name] = this->materialList.size() - 1;

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
    if (this->materialDict.find(materialName) == this->materialDict.end())
    {
        msg = "Elements::setMaterialComposition. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    this->materialList[this->materialDict[materialName]].setComposition(names, amounts);
}

void Elements::setMaterialComposition(const std::string & materialName, \
                                const std::map<std::string, double> & composition)
{
    std::string msg;
    if (this->materialDict.find(materialName) == this->materialDict.end())
    {
        msg = "Elements::setMaterialComposition. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    this->materialList[this->materialDict[materialName]].setComposition(composition);
}

const Material & Elements::getMaterial(const std::string & materialName)
{
    std::string msg;
    if (this->materialDict.find(materialName) == this->materialDict.end())
    {
        msg = "Elements::getMaterial. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    return this->materialList[this->materialDict[materialName]];
}

Material Elements::getMaterialCopy(const std::string & materialName)
{
    std::string msg;
    if (this->materialDict.find(materialName) == this->materialDict.end())
    {
        msg = "Elements::getMaterial. Non existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }
    return this->materialList[this->materialDict[materialName]];
}


void Elements::addMaterial(Material material)
{
    std::string msg;
    std::string materialName;

    materialName = material.getName();

    if (this->materialDict.find(materialName) != this->materialDict.end())
    {
        msg = "Elements::addMaterial. Already existing material: " +  materialName;
        throw std::invalid_argument(msg);
    }

    this->materialDict[materialName] = this->materialList.size();
    // TODO: Make sure the material can be interpreted in terms of the supplied composition
    // because the composition can include other materials.
    // If made that way, the internal list of materials will be "clean" of references from
    // one material to another. Is it needed?
    this->materialList.push_back(material);
}

std::map<std::string, double> Elements::getComposition(const std::string & name)
{
    std::string msg;
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

    // check if it is a material
    if (this->materialDict.find(name) == this->materialDict.end())
    {
        // result at this point must be empty, we can send it back.
        return result;
    }

    // make sure we give back the elemental material composition
    tmpResult = this->materialList[this->materialDict[name]].getComposition();
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


std::map<std::string, double> Elements::getCompositionFromFormula(const std::string & formula)
{
    std::map<std::string, double> parsedFormula;
    std::map<std::string, double>::iterator it;
    std::string name, msg;
    double total;
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
        if (this->elementDict.find(name) == this->elementDict.end())
        {
            // msg = "Name " + name + " not among defined elements";
            // std::cout << msg << std::endl;
            parsedFormula.clear();
            return parsedFormula;
        }
        // multiply the number of atoms by the atomic weight
        parsedFormula[name] *= this->elementList[this->elementDict[name]].getAtomicMass();
        total += parsedFormula[name];
    }
    for (it = parsedFormula.begin(); it != parsedFormula.end(); ++it)
    {
        parsedFormula[it->first] /= total;
    }
    return parsedFormula;
}


std::map<std::string, double> Elements::parseFormula(const std::string & formula)
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
            if (!this->StringToDouble(formula.substr(p2+1, i - (p2+1)), factor))
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
                        if (this->StringToDouble(lastNumber, factor))
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
                if (this->StringToDouble(lastNumber, factor))
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

bool Elements::StringToDouble(const std::string& str, double& number)
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
