#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2016 European Synchrotron Radiation Facility
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
#include "fisx_epdl97.h"
#include "fisx_simplespecfile.h"
#include <stdexcept>
#include <iostream>
#include <math.h>
// #include <ctime>

namespace fisx
{

EPDL97::EPDL97()
{
    this->initialized = false;
    this->bindingEnergiesFile = "Unknown";
    this->crossSectionsFile = "Unknown";
    this->bindingEnergy.clear();
    this->muInputLabels.clear();
    this->muInputValues.clear();
    this->muLabelToIndex.clear();
    this->muEnergy.clear();
}

EPDL97::EPDL97(std::string directoryName)
{
    this->bindingEnergy.clear();
    this->muInputLabels.clear();
    this->muInputValues.clear();
    this->muLabelToIndex.clear();
    this->muEnergy.clear();
    // No check on directoryName.
    this->setDataDirectory(directoryName);
}

void EPDL97::setDataDirectory(std::string directoryName)
{
    this->bindingEnergy.clear();
    this->muInputLabels.clear();
    this->muInputValues.clear();
    this->muLabelToIndex.clear();
    this->muEnergy.clear();
    this->initialized = false;
    this->bindingEnergiesFile = "Unknown";
    this->crossSectionsFile = "Unknown";
    this->loadData(directoryName);

}
void EPDL97::loadData(std::string directoryName)
{
    std::string BINDING_ENERGIES="EADL97_BindingEnergies.dat";
    std::string CROSS_SECTIONS="EPDL97_CrossSections.dat";
    std::string joinSymbol;
    std::string filename;

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

    // Load the binding energies
    filename = directoryName + joinSymbol + BINDING_ENERGIES;
    this->loadBindingEnergies(filename);

    // Load the cross sections
    filename = directoryName + joinSymbol + CROSS_SECTIONS;
//     clock_t startTime = clock();
    this->loadCrossSections(filename);
//     std::cout << "CROSS SECTIONS " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << std::endl;

    //
    this->directoryName = directoryName;

    //everything went fine
    this->initialized = true;
}

// Binding energies
void EPDL97::loadBindingEnergies(std::string fileName)
{
    SimpleSpecfile sf;
    int    nScans;

    std::vector<std::string> tmpLabels;
    std::vector<std::string>::size_type nLabels, n;

    std::vector<std::vector<double> > tmpValues;
    std::vector<std::map<std::string, double> >::size_type i;
    std::string key;
    std::string msg;

    sf = SimpleSpecfile(fileName);
    nScans = sf.getNumberOfScans();
    if (nScans != 1)
    {
        msg = "EPDL97: Number of scans not equal one in binding energies file " + \
               fileName;
        throw std::ios_base::failure(msg);
    }

    tmpLabels = sf.getScanLabels(0);
    tmpValues = sf.getScanData(0);
    nLabels = tmpLabels.size();
    if (tmpValues[0].size() != nLabels)
    {
        std::cout << fileName;
        std::cout << " nLabels = " << nLabels;
        std::cout << " nValues = " << tmpValues[0].size();
        throw std::ios_base::failure("EPDL97: Number of values does not match number of labels");
    }

    this->bindingEnergy.resize(tmpValues.size());

    for (i = 0; i < this->bindingEnergy.size(); i++)
    {
        for (n = 1; n < nLabels; n++)
        {
            if (tmpLabels[n].substr(0, 1) == "K")
            {
                // single character
                key = "K";
            }
            else
            {
                if (tmpLabels[n].size() < 3)
                {
                    // two characters
                    key = tmpLabels[n].substr(0, 2);
                }
                else
                {
                    if (tmpLabels[n].substr(3, 1) == "(")
                    {
                        // three characters
                        key = tmpLabels[n].substr(0, 3);
                    }
                    else
                    {
                        // two characters
                        key = tmpLabels[n].substr(0, 2);
                    }
                }
            }
            this->bindingEnergy[i][key] = tmpValues[i][n];
        }
    }
    this->bindingEnergiesFile = fileName;
    // TODO, for a complete initialization some attenuation coefficients dhould be there
    this->initialized = true;
}

void EPDL97::loadCrossSections(std::string fileName)
{
    SimpleSpecfile sf;
    int    nScans;
    int i, n;
    std::vector<std::string> tmpLabels;
    std::vector<std::string>::size_type nLabels;
    std::vector<std::vector<double> >::size_type j;
    std::map<std::string, int> labelMap;
    std::string tmpString, key;
    std::string interestingLabels[15] = {"energy", "compton", "coherent", "photoelectric", "total",\
                                    "K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};
    std::vector<double>    *pVec;

    sf = SimpleSpecfile(fileName);
    nScans = sf.getNumberOfScans();
    if (nScans < 99)
    {
        std::cout << "Reading filename " << fileName << " Number of scans = " << nScans << std::endl;
        throw std::ios_base::failure("EPDL97: Not enough scans in cross sections file");
    }

    //get the labels from the first scan
    this->muInputLabels = sf.getScanLabels(0);
    nLabels = this->muInputLabels.size();

    if(nLabels < 15)
    {
        throw std::ios_base::failure("EPDL97: Not enough labels in cross sections file");
    }

    //identify the label indices containing the relevant information for us:
    //energy, coherent, photoelectric, compton, K, L1, L2, L3, M1, M2, M3, M4, M5 and total cross sections
    for (i = 0; i < 15; i++)
    {
        this->muLabelToIndex[interestingLabels[i]] = -1;
        key = this->toUpperCaseString(interestingLabels[i]);
        if (key == "ALL OTHER")
        {
            key = "ALLOTHER";
        }
        n = 0;
        while ((this->muLabelToIndex[interestingLabels[i]] == -1) && (n < (int) nLabels))
        {
            tmpString = this->toUpperCaseString(this->muInputLabels[n]);
            if (tmpString.find(key) != std::string::npos)
            {
                if (key == "COHERENT")
                {
                    // distinguish coherent and incoherent
                    if (tmpString.find("INCOHERENT") == std::string::npos)
                    {
                        this->muLabelToIndex[interestingLabels[i]] = n;
                    }
                }
                else
                {
                    if (key.size() < 3)
                    {
                        if(key.substr(0, key.size()) == tmpString.substr(0, key.size()))
                        {
                            this->muLabelToIndex[interestingLabels[i]] = n;
                        }
                    }
                    else
                    {
                        this->muLabelToIndex[interestingLabels[i]] = n;
                    }
                }
            }
            n++;
        }
        if(this->muLabelToIndex[interestingLabels[i]] == -1)
        {
            throw std::ios_base::failure("Cannot find mandatory cross section in file!");
        }
    }

    // we have got all the mapping, now read the needed data
    this->muInputValues.resize(nScans);
    this->muEnergy.resize(nScans);

    for (i = 0; i < nScans; i++)
    {
        tmpLabels = sf.getScanLabels(i);
        // check if the labels are the same (they should)
        if (tmpLabels.size() != nLabels)
        {
            std::cout << "Scan " << i <<\
                " does not have the same amount of labels as " << (i + 1) << std::endl;
            throw std::length_error("EPDL97: All scans do not have the same number of labels");
        }
        for (j = 0; j < nLabels; j++)
        {
            if (tmpLabels[j] != this->muInputLabels[j])
            {
                throw std::length_error("EPDL97: All scans do not have the same labels");
            }
        }
        // read the data
        this->muInputValues[i] = sf.getScanData(i);
        this->muEnergy[i].resize(this->muInputValues[i].size());
        for (j = 0; j < this->muInputValues[i].size(); j++)
        {
            pVec = &(this->muInputValues[i][j]);
            this->muEnergy[i][j] = (*pVec)[this->muLabelToIndex["energy"]];
            // recalculate the photoelectric effect of the non considered shells
            (*pVec)[this->muLabelToIndex["photoelectric"]] = (*pVec)[this->muLabelToIndex["total"]]-\
                                         (*pVec)[this->muLabelToIndex["compton"]]-\
                                         (*pVec)[this->muLabelToIndex["coherent"]];

            if ((*pVec)[this->muLabelToIndex["photoelectric"]] > 0.0)
            {
                if (((*pVec)[this->muLabelToIndex["all other"]] > 0.0) && (i > 17))
                {

                    // there are higher shells than the M5 excited
                    (*pVec)[this->muLabelToIndex["all other"]] = (*pVec)[this->muLabelToIndex["photoelectric"]]-\
                                             (*pVec)[this->muLabelToIndex["K"]]-\
                                             (*pVec)[this->muLabelToIndex["L1"]]-\
                                             (*pVec)[this->muLabelToIndex["L2"]]-\
                                             (*pVec)[this->muLabelToIndex["L3"]]-\
                                             (*pVec)[this->muLabelToIndex["M1"]]-\
                                             (*pVec)[this->muLabelToIndex["M2"]]-\
                                             (*pVec)[this->muLabelToIndex["M3"]]-\
                                             (*pVec)[this->muLabelToIndex["M4"]]-\
                                             (*pVec)[this->muLabelToIndex["M5"]];
                }
                else
                {
                    (*pVec)[this->muLabelToIndex["all other"]] = 0.0;
                }
            }
            else
            {
                // take care of rounding
                (*pVec)[this->muLabelToIndex["photoelectric"]] = 0.0;
            }
            // take care of rounding
            if((*pVec)[this->muLabelToIndex["all other"]] < 0.0)
            {
                (*pVec)[this->muLabelToIndex["all other"]] = 0.0;
            }
        }
    }
    this->crossSectionsFile = fileName;
}

void EPDL97::setBindingEnergies(const int & z, const std::map<std::string, double> & bindingEnergies)
{
    std::map<std::string, double>::iterator it;
    std::string tmpString;

    if(z < 1)
    {
        throw std::runtime_error("EPDL97 Atomic number should be positive");
    }

    this->bindingEnergy[z - 1].clear();
    this->bindingEnergy[z - 1] = bindingEnergies;
}

const std::map<std::string, double> & EPDL97::getBindingEnergies(const int & z)
{
    if(z < 1)
    {
        throw std::runtime_error("EPDL97 Atomic number should be positive");
    }
    if (z >= (int) this->bindingEnergy.size())
    {
        // Repeat for the elements beyond the last one
        return this->bindingEnergy[this->bindingEnergy.size() - 1];
    }
    else
    {
        return this->bindingEnergy[z - 1];
    }
}


std::map<std::string, double> EPDL97::getMassAttenuationCoefficients(const int & z, const double & energy) const
{
    std::pair<long, long> indices;
    long i1, i2, i1w, i2w;
    double A, B, x0, x1, y0, y1, Aw, Bw, x0w, x1w;
    std::map<std::string, int>::const_iterator c_it;
    std::string key;
    int zHelp, idx;
    std::map<std::string, double> result;
    std::map<std::string, double>::const_iterator cStrDoubleIt;
    const std::vector<std::vector<double> > *pVector;

    if(!this->initialized)
    {
        throw std::runtime_error("EPDL97 Mass attenuation coefficients not initialized");
    }

    if(z < 1)
    {
        throw std::runtime_error("EPDL97 Atomic number should be positive");
    }

    zHelp = z - 1;

    if(zHelp > (int) (this->muEnergy.size() - 1))
    {
        // std::cout << "WARNING: Using data from last available element" << std::endl;
        zHelp = (int) (this->muEnergy.size() - 1);
    }


    indices = this->getInterpolationIndices(this->muEnergy[zHelp], energy);

    i1 = indices.first;
    i2 = indices.second;

    x0 = this->muEnergy[zHelp][i1];
    x1 = this->muEnergy[zHelp][i2];
    // std::cout << "EPDL97 i1, i2 " << i1 << " " << i2 <<std::endl;
    // std::cout << "EPDL97 x0, x1 " << x0 << " " << x1 <<std::endl;
    if (energy == x1)
    {
        if ((i2 + 1) < ((int) this->muEnergy[zHelp].size()))
        {
            if (this->muEnergy[zHelp][i2+1] == x1)
            {
                // repeated energy
                i1 = i2;
                i2++;
                x0 = this->muEnergy[zHelp][i1];
                x1 = this->muEnergy[zHelp][i2];
                // std::cout << "RETOUCHED EPDL97 i1, i2 " << i1 << " " << i2 <<std::endl;
                // std::cout << "RETOUCHED EPDL97 x0, x1 " << x0 << " " << x1 <<std::endl;
            }
        }
    }
    result["energy"] = energy;
    if ((i1 == i2) ||((x1 - x0) < 5.E-10))
    {
        // std::cout << "case a" <<std::endl;
        // std::cout << "x0, x1 " << x0 << " " << x1 << " energy = " << energy <<std::endl;
        for (c_it = this->muLabelToIndex.begin(); c_it != this->muLabelToIndex.end(); ++c_it)
        {
            key = c_it->first;
            idx = c_it->second;
            pVector = &(this->muInputValues[zHelp]);
            // std::cout << "key " << key << std::endl;
            if ((key == "total") || (key == "photoelectric") || (key == "energy"))
            {
                continue;
            }
            if ((key == "coherent") || (key == "compton") || (key == "all other"))
            {
                result[key] = (*pVector)[i2][idx];
            }
            else
            {
                if(key.size() < 3)
                {
                    // we are dealing with a shell
                    result[key] = 0.0;
                    cStrDoubleIt = this->bindingEnergy[zHelp].find(key);
                    if (cStrDoubleIt == this->bindingEnergy[zHelp].end())
                    {
                        std::cout << "Key not found " << key << std::endl;
                        throw std::runtime_error("Key not found");
                    }
                    if ((energy >=  cStrDoubleIt->second) && \
                       (cStrDoubleIt->second > 0.0))
                    {
                        // the shell is excited
                        if ((*pVector)[i1][idx] > 0.0)
                        {
                            result[key] = (*pVector)[i1][idx];
                        }
                        else
                        {
                            if (((x1 - x0) < 5.E-10) && ((*pVector)[i2][idx] > 0.0))
                            {
                                result[key] = (*pVector)[i2][idx];
                            }
                             else
                             {
                                // according to the binding energies, the shell is excited, but the
                                // respective mass attenuation is zero. We have to extrapolate
                                i1w = i1;
                                while((*pVector)[i1w][idx] <= 0.0)
                                {
                                    i1w += 1;
                                }
                                i2w = i1w + 1;
                                y0 = (*pVector)[i1w][idx];
                                y1 = (*pVector)[i2w][idx];
                                x0w = this->muEnergy[zHelp][i1w];
                                x1w = this->muEnergy[zHelp][i2w];
                                B = 1.0 / log( x1w / x0w);
                                A = log(x1w/energy) * B;
                                B *= log( energy / x0w);
                                result[key] = exp(A * log(y0) + B * log(y1));
                            }
                        }
                    }
                }
            }
        }
        // std::cout << "case a passed" <<std::endl;
    }
    else
    {
        // y = exp(( log(y0)*log(x1/x) + log(y1)*log(x/x0)) / log(x1/x0))
        B = 1.0 / log( x1 / x0);
        A = log(x1/energy) * B;
        B *= log( energy / x0);

        for (c_it = this->muLabelToIndex.begin(); c_it != this->muLabelToIndex.end(); ++c_it)
        {
            key = c_it->first;
            idx = c_it->second;
            pVector = &(this->muInputValues[zHelp]);
            if ((key == "total") || (key == "photoelectric")|| (key == "energy"))
            {
                // they will be deduced from the others
                continue;
            }
            if ((key == "coherent") || (key == "compton") || (key == "all other"))
            {
                y0 = (*pVector)[i1][idx];
                y1 = (*pVector)[i2][idx];
                if ((y0 > 0.0) && (y1 > 0.0))
                {
                    result[key] = exp(A * log(y0) + B * log(y1));
                }
                else
                {
                    if ((y1 > 0.0) && ((energy - x0) > 1.E-5))
                    {
                        result[key] = exp(B * log(y1));
                    }
                    else
                    {
                        result[key] = 0.0;
                    }
                }
                continue;
            }
            if(key.size() < 3)
            {
                // we are dealing with a shell
                result[key] = 0.0;
                cStrDoubleIt = this->bindingEnergy[zHelp].find(key);
                if (cStrDoubleIt == this->bindingEnergy[zHelp].end())
                {
                    std::cout << "Key not found " << key << std::endl;
                    throw std::runtime_error("Key not found");
                }
                if ((energy >= cStrDoubleIt->second) && \
                    (cStrDoubleIt->second > 0.0))
                {
                    if ((*pVector)[i1][idx] > 0.0)
                    {
                        // usual interpolation case
                        // the shell is excited and the photoelectric coefficient is positive
                        y0 = (*pVector)[i1][idx];
                        y1 = (*pVector)[i2][idx];
                        if ((y0 > 0.0) && (y1 > 0))
                        {
                            result[key] = exp(A * log(y0) + B * log(y1));
                        }
                    }
                    else
                    {
                        // extrapolation case
                        // We are forcing EPDL97 to respect a given set of binding energies
                        i1w = i1;
                        while((*pVector)[i1w][idx] <= 0.0)
                        {
                            i1w += 1;
                        }
                        i2w = i1w + 1;
                        y0 = (*pVector)[i1w][idx];
                        y1 = (*pVector)[i2w][idx];
                        x0w = this->muEnergy[zHelp][i1w];
                        x1w = this->muEnergy[zHelp][i2w];
                        Bw = 1.0 / log( x1w / x0w);
                        Aw = log(x1w/energy) * Bw;
                        Bw *= log( energy / x0w);
                        result[key] = exp(Aw * log(y0) + Bw * log(y1));
                    }
                }
            }
        }
    }
    result["pair"] = 0.0;
    result["photoelectric"] = result["K"] + result["L1"] + result["L2"] + result["L3"] +\
                (result["M1"] + result["M2"] + result["M3"] + result["M4"] + result["M5"] +\
                result["all other"]);
    result["total"] = result["photoelectric"] + result["coherent"] + result["compton"] + result["pair"];
    return result;
}

std::map<std::string, std::vector<double> > EPDL97::getMassAttenuationCoefficients(const int & z, \
                                                const std::vector<double> & energy) const
{
    std::vector<double>::size_type length, i;
    std::map<std::string, double> tmpResult;
    std::map<std::string, std::vector<double> > result;
    std::map<std::string, double>::const_iterator c_it;

    length = energy.size();

    for (i = 0; i < length; i++)
    {
        tmpResult = this->getMassAttenuationCoefficients(z, energy[i]);
        if (i == 0)
        {
            for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
            {
                result[c_it->first].resize(length);
            }
        }
        for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
        {
            result[c_it->first][i] = c_it->second;
        }
    }
    return result;
}


std::map< std::string, std::vector<double> > EPDL97::getMassAttenuationCoefficients(const int & z) const
{

    int i, idx, iMu, nValues;
    std::map<std::string, std::vector<double> > result;
    std::map<std::string, int>::const_iterator c_it;
    std::string key;
    const std::vector<std::vector<double> > *pVector;
    std::vector<double>  tmpVector;
    std::vector<double>  tmpPhotoelectricVector;
    std::vector<double>  nonPhotoelectricVector;

    if(!this->initialized)
    {
        throw std::runtime_error("EPDL97 Mass attenuation coefficients not initialized");
    }

    idx = z - 1;

    if(idx < 0)
    {
        throw std::runtime_error("EPDL97 Atomic number should be positive");
    }

    if(idx >= (int) this->muEnergy.size())
    {
        // return the information of last available element
        idx = (int) (this->muEnergy.size() - 1);
    }

    pVector = &(this->muInputValues[idx]);
    nValues = (int) this->muEnergy[idx].size();
    tmpVector.resize(nValues);
    tmpPhotoelectricVector.resize(nValues);
    nonPhotoelectricVector.resize(nValues);
    std::fill(tmpPhotoelectricVector.begin(), tmpPhotoelectricVector.end(), 0.0);
    std::fill(nonPhotoelectricVector.begin(), nonPhotoelectricVector.end(), 0.0);
    for (c_it = this->muLabelToIndex.begin(); c_it != this->muLabelToIndex.end(); ++c_it)
    {
        key = c_it->first;
        if ((key == "total") || (key == "photoelectric"))
        {
            continue;
        }
        iMu = c_it->second;
        for (i = 0; i < nValues; i++)
        {
            tmpVector[i] = (*pVector)[i][iMu];
            if ((key == "coherent") || (key == "compton") || (key == "pair"))
            {
                nonPhotoelectricVector[i] += tmpVector[i];
            }
            else
            {
                if ((key == "K") || (key == "L1") || (key == "L2") || (key == "L3") || \
                    (key == "M1") || (key == "M2") || (key == "M3") || (key == "M4") || \
                    (key == "M5") || (key == "all other"))
                {
                    tmpPhotoelectricVector[i] += tmpVector[i];
                }
            }
        }
        result[key] = tmpVector;
    }
    std::fill(tmpVector.begin(), tmpVector.end(), 0.0);
    if (this->muLabelToIndex.find("pair") == this->muLabelToIndex.end())
    {
        result["pair"] = tmpVector;
    }
    result["total"].resize(nValues);
    result["photoelectric"] = tmpPhotoelectricVector;
    for (i = 0; i < nValues; i++)
    {
        tmpVector[i] = nonPhotoelectricVector[i] + tmpPhotoelectricVector[i];
    }
    result["total"] = tmpVector;
    return result;
}

std::map<std::string, double> EPDL97::getPhotoelectricWeights(const int & z, \
                                                              const double & energy)
{
    // Given an excitation energy and an optional list of shells to consider
    // gives back the ratio mu(shell, energy)/mu(energy) where mu refers to the photoelectric
    // mass attenuation coefficient.
    // The special shell "all others" refers to all the shells not in the K, L or M groups.
    // For instance, for the K shell, it is the equivalent of (Jk-1)/Jk where Jk is the k jump.
    // The mass attenuation coefficients used in the calculations are also provided.
    std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};
    std::string key;
    int i;
    std::map<std::string, double> tmpResult;
    std::map<std::string, double>::const_iterator c_it;
    std::map<std::string, double> result;
    double muPhotoelectric;

    // get the mass attenuation coefficients
    tmpResult = this->getMassAttenuationCoefficients(z, energy);

    for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
    {
        key = "mu " + c_it->first;
        result[key] = c_it->second;
    }

    muPhotoelectric = tmpResult["photoelectric"];

    if (muPhotoelectric <= 1.0E-10)
    {
        for (i = 0; i < 10; i++)
        {
            result[shellList[i]] = 0.0;
        }
    }
    else
    {
        for (i = 0; i < 10; i++)
        {
            key = shellList[i];
            result[key] = tmpResult[key]/muPhotoelectric;
        }
    }

    return result;
}

std::map<std::string, std::vector<double> > EPDL97::getPhotoelectricWeights(const int & z, \
                                                              const std::vector<double> & energy)
{
    std::vector<double>::size_type length, i;
    std::map<std::string, double> tmpResult;
    std::map<std::string, std::vector<double> > result;
    std::map<std::string, double>::const_iterator c_it;

    length = energy.size();

    for (i = 0; i < length; i++)
    {
        tmpResult = this->getPhotoelectricWeights(z, energy[i]);
        if (i == 0)
        {
            for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
            {
                result[c_it->first].resize(length);
            }
        }
        for (c_it = tmpResult.begin(); c_it != tmpResult.end(); ++c_it)
        {
            result[c_it->first][i] = c_it->second;
        }
    }
    return result;
}


std::string EPDL97::toUpperCaseString(const std::string & str) const
{
    std::string::size_type i;
    std::string converted;
    for(i = 0; i < str.size(); ++i)
        converted += toupper(str[i]);
    return converted;
}

std::pair<long, long> EPDL97::getInterpolationIndices(const std::vector<double> & vec, const double & x) const
{
    static long lastI0 = 0L;
    std::vector<double>::size_type length, iMin, iMax, distance;
    int counter;
    std::pair<long, long> result;

    // try if last point is of any use
    length = vec.size();
    if (lastI0 >= (int) length)
    {
        lastI0 = (long) (length - 1);
    }
    if (x < vec[lastI0])
    {
        iMax = lastI0;
        iMin = 0;
    }
    else
    {
        iMin = lastI0;
        iMax = length - 1;
        // try a point that is close?
        if ((length - iMin) > 21)
        {
            lastI0 = (long) (iMin + 20);
            if (x < vec[lastI0])
            {
                iMax = lastI0;
            }
        }
    }

    counter = 0;
    distance = iMax - iMin;
    while (distance > 1)
    {
        // divide the distance by two
        if (distance > 2)
        {
            distance = distance >> 1;
        }
        else
        {
            distance -= 1;
        }
        lastI0 = (long) (iMin + distance);
        if (x > vec[lastI0])
        {
            iMin = lastI0;
        }
        else
        {
            iMax = lastI0;
        }
        distance = iMax - iMin;
        counter++;
    }
    // std::cout << "Needed " << counter << " iterations " << std::endl;
    result.first = (long) iMin;
    result.second = (long) iMax;
    return result;
}

} // namespace fisx
