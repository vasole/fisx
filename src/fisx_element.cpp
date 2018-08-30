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
#include "fisx_element.h"
#include "fisx_math.h"
#include <iostream>
#include <math.h>
#include <stdexcept>

namespace fisx
{

Element::Element()
{
    // Default constructor required
    this->name = "Unknown";

    // Undefined atomic number
    this->atomicNumber = 0;

    // Unset density
    this->density = 1.0;

    // initialize keys
    this->initPartialPhotoelectricCoefficients();

    // cascade cache
    this->cascadeCacheEnabledFlag = false;

    // calculation cache
    this->setCacheEnabled(0);
}

Element::Element(std::string name, int z = 0)
{
    // No check on element name. Anything can be an element and
    // can have any atomic number
    this->name = name;
    this->setAtomicNumber(z);
    // Unset density
    this->density = 1.0;
    this->initPartialPhotoelectricCoefficients();

    // cascade cache
    this->cascadeCacheEnabledFlag = false;

    // calculation cache
    this->setCacheEnabled(0);
}

void Element::setName(const std::string & name)
{
    this->name = name;
}

std::string Element::getName() const
{
    return  this->name;
}

void Element::setDensity(const double & density)
{
    this->density = density;
}

double Element::getDensity()
{
    return this->density;
}


// Binding energies
void Element::setBindingEnergies(std::map<std::string, double> bindingEnergies)
{
    std::map<std::string, double>::iterator it;
    std::string tmpString;

    // get rid of any shell definition
    this->shellInstance.clear();
    this->bindingEnergy.clear();

    for (it = bindingEnergies.begin(); it != bindingEnergies.end(); ++it)
    {
        this->bindingEnergy[it->first] = it->second;
        tmpString = "";
        if (it->first.size())
        {
            tmpString = it->first.substr(0, 1);
        }
        //std::cout << this->name << " " << it->first << " " << it->second << "tmp = " << tmpString << std::endl;
        if ((tmpString == "K") || (tmpString == "L") || (tmpString == "M"))
        {
            if (this->shellInstance.find(it->first) == this->shellInstance.end())
            {
                // instantiate the appropriate shell
                // Shell constructor will take care of valid shells
                // This uses default constructor ...
                this->shellInstance[it->first] = Shell(it->first);
                /* This does not require default constructor, but MSVC still requires it
                 this->shellInstance.insert\
                    (std::map<std::string, Shell>::value_type(it->first, Shell::Shell(it->first))); */
            }
        }
    }
}

void Element::setBindingEnergies(std::vector<std::string> labels, std::vector<double> bindingEnergies)
{
    std::map<std::string, double> inputData;
    std::vector<std::string>::size_type it;
    std::vector<double>::size_type i;

    if (labels.size() != bindingEnergies.size())
    {
        throw std::invalid_argument("Number of labels does not match number of energies");
    }
    i = 0;
    for (it = 0; it < labels.size(); ++it)
    {
        inputData[labels[it]] = bindingEnergies[i];
        i++;
    }
    this->setBindingEnergies(inputData);
}


const std::map<std::string, double> & Element::getBindingEnergies() const
{
    return this->bindingEnergy;
}


// Mass attenuation coefficients
void Element::setMassAttenuationCoefficients(const std::vector<double> & energies, \
                                        const std::vector<double> & photoelectric, \
                                        const std::vector<double> & coherent, \
                                        const std::vector<double> & compton)
{
    std::vector<double> pair;
    this->setMassAttenuationCoefficients(energies, photoelectric, coherent, compton, pair);
}

void Element::setMassAttenuationCoefficients(const std::vector<double> & energies, \
                                        const std::vector<double> & photoelectric, \
                                        const std::vector<double> & coherent, \
                                        const std::vector<double> & compton, \
                                        const std::vector<double> & pair = std::vector<double>())
{
    std::string msg;
    std::vector<double>::const_iterator c_it;
    std::vector<double>::size_type length, i, pairLength;
    std::map<std::string, std::vector<double> >::iterator mu_it;

    // energies are expected in keV and ordered
    length = energies.size();

    if (photoelectric.size() != length)
    {
        msg = "setMassAttenuationCoefficients: Photoelectric data size not equal to energies data size";
        throw std::invalid_argument(msg);
    }

    if (compton.size() != length)
    {
        msg = "setMassAttenuationCoefficients: Compton data size not equal to energies data size";
        throw std::invalid_argument(msg);
    }

    if (coherent.size() != length)
    {
        msg = "setMassAttenuationCoefficients: Coherent data size not equal to energies data size";
        throw std::invalid_argument(msg);
    }

    pairLength = pair.size();
    if (pairLength > 0)
    {
        if (pairLength != length)
        {
            msg = "setMassAttenuationCoefficients: Pair data size not equal to energies data size";
            throw std::invalid_argument(msg);
        }
    }

    // check energies are supplied in ascending order
    for (i = 0; i < length; i++)
    {
        if (i > 0)
        {
            if (energies[i] < energies[i-1])
            {
                std::cout << "ELEMENT " << this->name << std::endl;
                std::cout << energies[i] << " < " << energies[i-1] << std::endl;
                throw std::invalid_argument("Energies have to be supplied in ascending order");
            }
        }
    }
    this->clearCache();

    if (this->mu.size() > 0)
    {
        for (mu_it = this->mu.begin(); mu_it != this->mu.end(); ++mu_it)
        {
            this->mu[mu_it->first].clear();
        }
        this->mu.clear();
    }

    this->mu["coherent"] = std::vector<double> (coherent);
    this->mu["compton"] = std::vector<double> (compton);
    this->mu["energy"] = std::vector<double> (energies);
    // TODO get rid of this vector?
    this->muEnergy = std::vector<double> (energies);
    if (pairLength > 0)
    {
        this->mu["pair"] = std::vector<double> (pair);
    }
    else
    {
        this->mu["pair"].resize(length);
        for (i = 0; i < length; i++)
        {
            this->mu["pair"][i] = 0.0;
        }
    }
    this->mu["photoelectric"] = std::vector<double> (photoelectric);
    this->mu["total"] = std::vector<double> (coherent);
    for (i = 0; i < length; i++)
    {
        this->mu["total"][i] += this->mu["compton"][i] +\
                                this->mu["pair"][i] + this->mu["photoelectric"][i];
    }
}

void Element::setTotalMassAttenuationCoefficient(const std::vector<double> & energies, \
                                                 const std::vector<double> & total)
{
    // The goal of this method is to provide experimentally measured total mass attenuation coefficients
    // Since the coherent and the Compton cross sections seem to be in good agreement among the different
    // compilations, we assume that the difference total - coherent - compton - pair = photoelectric
    // So, we have to:
    //  0 - check the attenuation coefficents for the individual processes have been entered
    //    1 - check the number of energies is equal to the number of coefficients
    //    2 - merge the energy grids
    //  3 - interpolate compton, coherent and pair at these energies
    throw std::runtime_error("setTotalMassAttenuationCoefficient not implemented yet");
}

const std::map< std::string, std::vector<double> > & Element::getMassAttenuationCoefficients() const
{
    //TODO check initialization
    return this->mu;
}

std::map<std::string, double> Element::getMassAttenuationCoefficients(const double & energy) const
{
    std::pair<long, long> indices;
    long i1, i2;
    double A, B, x0, x1, y0, y1;
    //std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};
    std::map<std::string, std::vector<double> >::const_iterator c_it;
    std::string key;
    std::map< double, std::map< std::string, double> >::const_iterator c_it2;
    std::map<std::string, double> result;

    if (this->muEnergy.size() < 1)
    {
        throw std::runtime_error("Mass attenuation coefficients not initialized yet!");
    }

    if (this->isCacheEnabled())
    {
        c_it2 = this->muCache.find(energy);
        if (c_it2 != this->muCache.end())
        {
            return c_it2->second;
        }
    }

    // TODO: if the partial are not given, use the total photoelectric

    // calculate the partial photoelectric mass attenuation coefficients
    for (c_it = this->muPartialPhotoelectricEnergy.begin();
         c_it != this->muPartialPhotoelectricEnergy.end(); ++c_it)
    {
        if (c_it->second.size() > 0)
        {
            // partial initialized at least for one shell
            //std::cout << " Getting partial photoelectric mass attenuation coefficients " << std::endl;
            result = this->getPartialPhotoelectricMassAttenuationCoefficients(energy);
            break;
        }
    }

    // std::cout << "Calling interpolation" <<std::endl;
    indices = this->getInterpolationIndices(this->muEnergy, energy);

    i1 = indices.first;
    i2 = indices.second;

    x0 = this->muEnergy[i1];
    x1 = this->muEnergy[i2];

    // std::cout << "ELEMENT i1, i2 " << i1 << " " << i2 <<std::endl;
    // std::cout << "ELEMENT x0, x1 " << x0 << " " << x1 <<std::endl;

    if (energy == x1)
    {
        if ((i2 + 1) < ((int) this->muEnergy.size()))
        {
            if (this->muEnergy[i2+1] == x1)
            {
                // repeated energy
                i1 = i2;
                i2++;
                x0 = this->muEnergy[i1];
                x1 = this->muEnergy[i2];
                // std::cout << "RETOUCHED ELEMENT i1, i2 " << i1 << " " << i2 <<std::endl;
                // std::cout << "RETOUCHED ELEMENT x0, x1 " << x0 << " " << x1 <<std::endl;
            }
        }
    }


    result["energy"] = energy;
    if ((i1 == i2) ||((x1 - x0) < 5.E-10))
    {
        // std::cout << "case a" <<std::endl;
        //std::cout << "x0, x1 " << x0 << " " << x1 << " energy = " << energy <<std::endl;
        for (c_it = this->mu.begin(); c_it != this->mu.end(); ++c_it)
        {
            key = c_it->first;
            if ((key == "coherent") || (key == "compton") || (key == "pair"))
            {
                result[key] = c_it->second[i1];
            }
        }
    }
    else
    {
        // y = exp(( log(y0)*log(x1/x) + log(y1)*log(x/x0)) / log(x1/x0))
        // std::cout << "case b" <<std::endl;
        B = 1.0 / log( x1 / x0);
        A = log(x1/energy) * B;
        B *= log( energy / x0);
        for (c_it = this->mu.begin(); c_it != this->mu.end(); ++c_it)
        {
            key = c_it->first;
            //std::cout << "key " << key << std::endl;
            if ((key == "coherent") || (key == "compton") || (key == "pair"))
            {
                // we are left with coherent, compton and pair
                y0 = c_it->second[i1];
                y1 = c_it->second[i2];
                //std::cout << "y0, y1 " << y0 << " " << y1 <<std::endl;
                if ((y0 > 0.0) && (y1 > 0.0))
                {
                    // std::cout << "standard case" <<std::endl;
                    result[key] = exp(A * log(y0) + B * log(y1));
                    // std::cout << "entered value = " << result[key] <<std::endl;
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
            }
        }
    }
    result["photoelectric"] = result["K"] + result["L1"] + result["L2"] + result["L3"] +\
                (result["M1"] + result["M2"] + result["M3"] + result["M4"] + result["M5"] +\
                result["all other"]);

    result["total"] = result["photoelectric"] + result["coherent"] + result["compton"] + result["pair"];
    if (!Math::isFiniteNumber(result["total"]))
    {
        std::cout << "element = " << this->name << std::endl;
        std::cout << "energy = " << energy << std::endl;
        std::cout << "Photo = " << result["photoelectric"] << std::endl;
        std::cout << "coherent = " << result["coherent"] << std::endl;
        std::cout << "compton = " << result["compton"] << std::endl;
        std::cout << "pair = " << result["pair"] << std::endl;
        throw std::runtime_error("Invalid total mass attenuation coefficient");
    }
    return result;
}

std::map<std::string, std::vector<double> > Element::getMassAttenuationCoefficients(\
                                                const std::vector<double> & energy) const
{
    std::vector<double>::size_type length, i;
    std::map<std::string, double> tmpResult;
    std::map<std::string, std::vector<double> > result;
    std::map<std::string, double>::const_iterator c_it;

    length = energy.size();

    for (i = 0; i < length; i++)
    {
        tmpResult = this->getMassAttenuationCoefficients(energy[i]);
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

std::map<std::string, std::pair<double, int> > Element::extractEdgeEnergiesFromMassAttenuationCoefficients()
{
    if(this->mu["photoelectric"].size() < 1)
    {
        throw std::runtime_error("Photoelectric mass attenuation coefficients not initialized");
    }
    return this->extractEdgeEnergiesFromMassAttenuationCoefficients(this->mu["energy"], \
                                                                    this->mu["photoelectric"]);
}

std::map<std::string, std::pair<double, int> > \
        Element::extractEdgeEnergiesFromMassAttenuationCoefficients(const std::vector<double> & muEnergy,\
                                                                const std::vector<double> & muPhotoelectric)
{
    // This function tries to figure out the energies used by the set of photoelectric mass
    // attenuation coefficients in order to be able to calculate the partial photoelectric
    // cross sections of the K, Li and Mi shells.
    // It is based on the supplied mass attenuation coefficents having duplicated energies
    // corresponding at the edge energies. That is the common approach used by the XCOM and
    // the EPDL97 compilations among others.

    std::map<std::string, std::pair<double, int> >result;
    std::map<std::string, std::pair<double, int> >::iterator result_it;
    std::vector<double>::size_type i;
    std::string shellList[16] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5",\
                                "N1", "N2", "N3", "N4", "N5", "N6", "N7"};
    std::string key;
    double energy;
    std::map<double, std::vector<double>::size_type> candidates;
    std::map<double, std::vector<double>::size_type>::const_iterator c_it;

    int    shellIndex;


    if(muPhotoelectric.size() < 1)
    {
        throw std::runtime_error("Empty photoelectric mass attenuation vector");
    }
    if(muEnergy.size() < 1)
    {
        throw std::runtime_error("Empty energy mass attenuation vector");
    }
    if(muEnergy.size() != muPhotoelectric.size())
    {
        throw std::runtime_error("Supplied vectors do not have the same length");
    }

    if (this->bindingEnergy.size() < 1)
    {
        throw std::runtime_error("Binding energies not initialized");
    }

    for(i = 0; i < (muPhotoelectric.size() - 1); i++)
    {
        energy = muEnergy[i];
        if (energy == muEnergy[i + 1])
        {
            // We have a repeated energy
            if (muPhotoelectric[i] < muPhotoelectric[i + 1])
            {
                // we have an energy edge
                candidates[energy] = i;
            }
        }
    }

    if (candidates.size() > 0)
    {
        c_it = candidates.end();
        for (shellIndex = 0; shellIndex < 9; shellIndex++)
        {
            --c_it;
            energy = c_it->first;
            key = shellList[shellIndex];
            if (this->bindingEnergy[key] > 0.0)
                {
                // check the energy is correct within 100 eV
                if(fabs(energy - this->bindingEnergy[key]) < 0.100)
                {
                    result[key].first = energy;
                    result[key].second = (int) c_it->second;
                }
            }
            if (c_it == candidates.begin())
            {
                shellIndex = 9;
            }
        }
    }

    /*
    if (false)
    {
        for (result_it = result.begin(); result_it != result.end(); ++result_it)
        {
            key = result_it->first;
            energy = result_it->second.first;
            std::cout << this->name << " Found shell " << key << " at " << energy << std::endl;
        }
    }
    */
    return result;
}

// Partial photoelectric cross sections

// initialize keys
void Element::initPartialPhotoelectricCoefficients()
{
    std::string photoShells[10] = {"K", "L1", "L2", "L3", \
                                   "M1", "M2", "M3", "M4", "M5", "all other"};
    long i;

    this->clearCache();
    for (i=0; i < 10; i++)
    {
        // This creates (if it does not exist) and clears if not empty
        this->muPartialPhotoelectricEnergy[photoShells[i]].clear();
        this->muPartialPhotoelectricValue[photoShells[i]].clear();
    }
}


void Element::setPartialPhotoelectricMassAttenuationCoefficients(const std::string & shell, \
                                                        const std::vector<double> & energy, \
                                                        const std::vector<double> & partialPhotoelectric)
{
    std::string msg;
    std::vector<double>::size_type i, length;
    double lastEnergy;

    if (this->muPartialPhotoelectricEnergy.find(shell) == this->muPartialPhotoelectricEnergy.end())
    {
        msg = "Shell has to be one of K, L1, L2, L3, M1, M2, M3, M4, M5, all other. Got <" + shell +">";
        throw std::invalid_argument(msg);
    }

    length = energy.size();
    if (length != partialPhotoelectric.size())
    {
        throw std::invalid_argument("Number of energies and of coefficients do not match");
    }

    lastEnergy = 0.0;
    for (i = 0; i < length; i++)
    {
        if (energy[i] < lastEnergy)
        {
            std::cout << "ELEMENT " << this->name << std::endl;
            std::cout << energy[i] << " < " << lastEnergy << std::endl;
            throw std::invalid_argument("Partial photoelectric energies should be in ascending order");
        }
        else
        {
            lastEnergy = energy[i];
        }
    }

    // checks finished, we can go ahead
    this->clearCache();
    this->muPartialPhotoelectricEnergy[shell].clear();
    this->muPartialPhotoelectricValue[shell].clear();

    this->muPartialPhotoelectricEnergy[shell] = std::vector<double>(energy);
    this->muPartialPhotoelectricValue[shell] = std::vector<double>(partialPhotoelectric);
    if (shell != "all other")
    {
        for (i = 1; i < length; i++)
        {
            if (this->muPartialPhotoelectricEnergy[shell][i] < this->bindingEnergy[shell])
            {
                this->muPartialPhotoelectricValue[shell][i] = 0.0;
            }
            else
            {
                // case of repeated values corresponding to an edge where the edge energy is above the set binding energy
                // for instance, exciting lead at 15.19 keV
                if (this->muPartialPhotoelectricEnergy[shell][i] == this->muPartialPhotoelectricEnergy[shell][i - 1])
                {
                    this->muPartialPhotoelectricEnergy[shell][i] += 0.000001;
                    this->muPartialPhotoelectricValue[shell][i - 1] = this->muPartialPhotoelectricValue[shell][i];
                }
            }
        }
    }
    //std::cout << this->muPartialPhotoelectricEnergy[shell][1100] << " " << this->muPartialPhotoelectricValue[shell][1100] << std::endl;
}

std::map<std::string, double> \
    Element::getPartialPhotoelectricMassAttenuationCoefficients(const double & energy) const
{
    std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};
    std::string shell;
    std::map<std::string, double> result;
    std::vector<std::string>::size_type i;
    std::pair<long, long> indices;
    long i1, i2, i1w, i2w;
    double A, B, x0, x1, y0, y1, x0w, x1w;
    std::map<std::string, double >::const_iterator c_itSingle;
    std::map<std::string, std::vector<double> >::const_iterator c_it;
    std::map<std::string, std::vector<double> >::const_iterator y_it;

    if (this->muPartialPhotoelectricEnergy.size() == 0)
    {
        throw std::runtime_error("Partial photoelectric cross sections not initialized");
    }

    // std::cout << " Calculating partials " << std::endl;
    // std::cout << "Entered partials for energy " << energy << std::endl;

    for (i = 0; i < 10; i ++)
    {
        shell = shellList[i];
        // std::cout << "shell " << shell << std::endl;
        result[shell] = 0.0;
        if (shell != "all other")
        {
            c_itSingle = this->bindingEnergy.find(shell);
            if ((c_itSingle->second == 0.0) || (energy < c_itSingle->second))
            {
                continue;
            }
        }
        c_it = this->muPartialPhotoelectricEnergy.find(shell);
        y_it = this->muPartialPhotoelectricValue.find(shell);
        indices = this->getInterpolationIndices(c_it->second, energy);
        i1 = indices.first;
        i2 = indices.second;
        x0 = c_it->second[i1];
        x1 = c_it->second[i2];
        /*
        if (energy == 15.19)
        {

            std::cout << shell << " partials i1, i2 " << i1 << " " << i2 <<std::endl;
            std::cout << " partials x0, x1 " << x0 << " " << x1 <<std::endl;
            std::cout << " values y[i1] " << y_it->second[i1];
            std::cout << " values y[i2] " << y_it->second[i2] << std::endl;
        }
        */
        if (energy == x1)
        {
            if ((i2 + 1) < ((int) c_it->second.size()))
            {
                if (c_it->second[i2+1] == x1)
                {
                    // repeated energy
                    i1 = i2;
                    i2++;
                    x0 = c_it->second[i1];
                    x1 = c_it->second[i2];
                    //std::cout << "RETOUCHED PARTIAL i1, i2 " << i1 << " " << i2 <<std::endl;
                    //std::cout << "RETOUCHED PARTIAL x0, x1 " << x0 << " " << x1 <<std::endl;
                }
            }
        }


        if ((i1 == i2) || ((x1 - x0) < 5.E-10))
        {
            // std::cout << "case a " <<std::endl;
            if (shell == "all other")
            {
                result[shell] = y_it->second[i2];
            }
            else
            {
                y0 = y_it->second[i1];
                if ( y0 > 0.0)
                {
                    result[shell] = y0;
                }
                else
                {
                    y1 = y_it->second[i2];
                    if (((x1 - x0) < 5.E-10) && (y1 > 0.0))
                    {
                        result[shell] = y1;
                    }
                    else
                     {
                        // according to the binding energies, the shell is excited, but the
                        // respective mass attenuation is zero. We have to extrapolate
                        i1w = i1;
                        while(y_it->second[i1w] <= 0.0)
                        {
                            i1w += 1;
                        }
                        i2w = i1w + 1;
                        y0 = y_it->second[i1w];
                        y1 = y_it->second[i2w];
                        x0w = c_it->second[i1w];
                        x1w = c_it->second[i2w];
                        B = 1.0 / log( x1w / x0w);
                        A = log(x1w/energy) * B;
                        B *= log( energy / x0w);
                        result[shell] = exp(A * log(y0) + B * log(y1));
                    }
                }
            }
        }
        else
        {
            // std::cout << "case b " <<std::endl;
            B = 1.0 / log( x1 / x0);
            A = log(x1/energy) * B;
            B *= log( energy / x0);
            y0 = y_it->second[i1];
            y1 = y_it->second[i2];

            if (shell == "all other")
            {
                if ((y0 > 0.0) && (y1 > 0.0))
                {
                    result[shell] = exp(A * log(y0) + B * log(y1));
                }
                else
                {
                    if ((y1 > 0.0) && ((energy - x0) > 1.E-5))
                    {
                        result[shell] = exp(B * log(y1));
                    }
                    else
                    {
                        result[shell] = 0.0;
                    }
                }
            }
            else
            {
                // we are dealing with a shell
                if (y0 > 0.0)
                {
                    // std::cout << "case b1" << std::endl;
                    // usual interpolation case
                    // the shell is excited and the photoelectric coefficient is positive
                    result[shell] = exp(A * log(y0) + B * log(y1));
                }
                else
                {
                    // according to the binding energies, the shell is excited, but the
                    // respective mass attenuation is zero. We have to extrapolate
                    // std::cout << "case b2" << std::endl;
                    i1w = i1;
                    while(y_it->second[i1w] <= 0.0)
                    {
                        i1w += 1;
                    }
                    i2w = i1w + 1;
                    x0w = c_it->second[i1w];
                    x1w = c_it->second[i2w];
                    // check for duplicates
                    if ((x1w - x0w) < 1.0E-10)
                    {
                        i1w += 1;
                        i2w += 1;
                        x0w = c_it->second[i1w];
                        x1w = c_it->second[i2w];
                    }
                    y0 = y_it->second[i1w];
                    y1 = y_it->second[i2w];
                    x0w = c_it->second[i1w];
                    x1w = c_it->second[i2w];
                    result[shell] = exp(log(y0) + (log(y1 / y0) / log(x1w / x0w)) * log(energy/x0w));
                    /*
                    if (energy == 15.19)
                    {
                        std::cout << shell << " partials i1w, i2w " << i1w << " " << i2w <<std::endl;
                        std::cout << " partials x0, x1 " << x0w << " " << x1w <<std::endl;
                        std::cout << " values y[i1] " << y_it->second[i1w];
                        std::cout << " values y[i2] " << y_it->second[i2w] << std::endl;
                        std::cout << " Final " << result[shell] << std::endl;
                    }
                    */
                }
            }
        }
        if (!Math::isFiniteNumber(result[shell]))
        {
            std::cout << "Element " << this->name << std::endl;
            std::cout << "energy " << energy << std::endl;
            std::cout << "shell " << shell << std::endl;
            std::cout << "i1 " << i1 << " i2 " << i2 << std::endl;
            std::cout << "A " << A << " B " << B << std::endl;
            std::cout << "x0 " << x0 << " x1 " << x1 << std::endl;
            std::cout << "y0 " << y0 << " y1 " << y1 << std::endl;
            throw std::runtime_error("Partial photoelectric coefficient is not finite");
        }
    }

    return result;
}

std::vector<std::string> Element::getExcitedShells(const double & energy) const
{
    // get the excited shells from the binding energies
    std::map<std::string, double>::const_iterator c_binding;
    std::string shell;
    std::vector<std::string> result;

    for(c_binding=this->bindingEnergy.begin();\
        c_binding!=this->bindingEnergy.end(); ++c_binding)
    {
        if (c_binding->second > 0.0)
        {
            if (energy > c_binding->second)
            {
                result.push_back(c_binding->first);
            }
        }
    }
    return result;
}

 std::map<std::string, std::vector<double> >Element::getInitialPhotoelectricVacancyDistribution(\
                                                const std::vector<double> & energies) const
{
    std::map<std::string, std::vector<double> > tmpMap;
    std::map<std::string, std::vector<double> > result;
    std::vector<double>::size_type i, j;
    std::string shell;
    double total;
    std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};

    tmpMap = this->getMassAttenuationCoefficients(energies);
    for (i = 0; i < 10; i++)
    {
        shell = shellList[i];
        result[shellList[i]].resize(tmpMap["total"].size());
        for (j = 0; j < tmpMap["total"].size(); j++)
        {
            total = tmpMap["photoelectric"][j];
            if ( total > 0.0)
            {
                result[shell][j] = tmpMap[shell][j] / total;
            }
            else
            {
                result[shell][j] = 0.0;
            }
        }
    }
    return result;
}

std::map<std::string, double> Element::getInitialPhotoelectricVacancyDistribution(\
                                                                const double & energy) const
{
    std::map<std::string, double> tmpMap;
    std::map<std::string, double> result;
    std::vector<double>::size_type i;
    std::string shell;
    double total;
    std::string shellList[10] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "all other"};

    tmpMap = this->getMassAttenuationCoefficients(energy);
    for (i = 0; i < 10; i++)
    {
        shell = shellList[i];
        total = tmpMap["photoelectric"];
        if ( total > 0.0 )
        {
            result[shell] = tmpMap[shell] / total;
        }
        else
        {
            result[shell] = 0.0;
        }
    }
    return result;
}

// Transitions
void Element::setRadiativeTransitions(std::string subshell, \
                                      std::vector<std::string> labels, std::vector<double> values)
{
    std::string msg;
    if (this->bindingEnergy.find(subshell) == this->bindingEnergy.end())
    {
        throw std::invalid_argument("Invalid shell");
    }
    if (this->bindingEnergy[subshell] <= 0.0)
    {
        msg = "Requested shell <" + subshell + "> has non positive binding energy";
        throw std::invalid_argument(msg);
    }
    if (this->shellInstance.find(subshell) == this->shellInstance.end())
    {
        msg = "Requested shell <" + subshell + "> is not a K, L or M subshell";
        throw std::invalid_argument(msg);
    }
    this->shellInstance[subshell].setRadiativeTransitions(labels, values);
}

void Element::setRadiativeTransitions(std::string subshell, std::map<std::string, double> values)
{
    if (this->bindingEnergy.find(subshell) == this->bindingEnergy.end())
    {
        throw std::invalid_argument("Invalid shell");
    }
    if (this->bindingEnergy[subshell] <= 0.0)
    {
        throw std::invalid_argument("Requested shell has non positive binding energy");
    }
    if (this->shellInstance.find(subshell) == this->shellInstance.end())
    {
        throw std::invalid_argument("Requested shell is not a K, L or M subshell");
    }
    this->shellInstance[subshell].setRadiativeTransitions(values);
    this->clearCache();
}

const std::map<std::string, double> & Element::getRadiativeTransitions(const std::string & subshell) const
{
    std::map<std::string, Shell>::const_iterator c_it;
    c_it = this->shellInstance.find(subshell);
    if (c_it == this->shellInstance.end())
    {
        throw std::invalid_argument("Requested shell is not a defined K, L or M subshell");
    }
    return c_it->second.getRadiativeTransitions();
}

void Element::setNonradiativeTransitions(std::string subshell, std::vector<std::string> labels, std::vector<double> values)
{
    if (this->bindingEnergy.find(subshell) == this->bindingEnergy.end())
    {
        throw std::invalid_argument("Invalid shell");
    }
    if (this->bindingEnergy[subshell] <= 0.0)
    {
        throw std::invalid_argument("Requested shell has non positive binding energy");
    }
    if (this->shellInstance.find(subshell) == this->shellInstance.end())
    {
        throw std::invalid_argument("Requested shell is not a K, L or M subshell");
    }
    this->shellInstance[subshell].setNonradiativeTransitions(labels, values);
    this->clearCache();
}

void Element::setNonradiativeTransitions(std::string subshell, std::map<std::string, double> values)
{
    if (this->bindingEnergy.find(subshell) == this->bindingEnergy.end())
    {
        throw std::invalid_argument("Invalid shell");
    }
    if (this->bindingEnergy[subshell] <= 0.0)
    {
        throw std::invalid_argument("Requested shell has non positive binding energy");
    }
    if (this->shellInstance.find(subshell) == this->shellInstance.end())
    {
        throw std::invalid_argument("Requested shell is not a K, L or M subshell");
    }
    this->shellInstance[subshell].setNonradiativeTransitions(values);
    this->clearCache();
}

const std::map<std::string, double> & Element::getNonradiativeTransitions(const std::string & subshell) const
{
    std::map<std::string, Shell>::const_iterator c_it;
    std::string msg;

    c_it = this->shellInstance.find(subshell);
    if ( c_it == this->shellInstance.end())
    {
        msg = "Requested shell <" + subshell + "> is not a defined K, L or M subshell";
        throw std::invalid_argument(msg);
    }

    return c_it->second.getNonradiativeTransitions();
}


void Element::setShellConstants(std::string subshell, std::map<std::string, double> constants)
{
    std::string msg;
    if (this->shellInstance.find(subshell) == this->shellInstance.end())
    {
        msg = "Requested shell <" + subshell + "> is not a defined K, L or M subshell";
        throw std::invalid_argument(msg);
    }
    this->shellInstance[subshell].setShellConstants(constants);
    this->emptyCascadeCache();
    this->clearCache();
}

const std::map<std::string, double> & Element::getFluorescenceRatios(const std::string & subshell) const
{
    const Shell & shell = this->getShell(subshell);
    return shell.getFluorescenceRatios();
}

const std::map<std::string, double> & Element::getAugerRatios(std::string subshell)
{
    Shell shell;
    shell = this->getShell(subshell);
    return shell.getAugerRatios();
}

const std::map<std::string, std::map<std::string, double> > & Element::getCosterKronigRatios(std::string subshell)
{
    Shell shell;
    shell = this->getShell(subshell);
    return shell.getCosterKronigRatios();
}

std::map<std::string, double> Element::getShellConstants(const std::string & subshell) const
{

    std::map<std::string, Shell>::const_iterator it;
    it = this->shellInstance.find(subshell);
    if (it == this->shellInstance.end())
    {
        throw std::invalid_argument("Requested shell is not a defined K, L or M subshell");
    }
    return it->second.getShellConstants();
}

const std::map<std::string, double> & Element::getXRayLines(const std::string & family) const
{
    return this->getFluorescenceRatios(family);
}

std::map<std::string, double> Element::getEmittedXRayLines(const double & energy) const
{
    std::string keys[9] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"};
    std::map<std::string, Shell>::const_iterator shell_it;
    std::map<std::string, double> fluorescenceRatios;
    std::map<std::string, double>::const_iterator c_it;
    std::map<std::string, double> result;
    double rate;
    int i;

    result.clear();
    for (i = 0; i < 9; i++)
    {
        // get all the fluorescence transitions to that shell
        shell_it = this->shellInstance.find(keys[i]);
        if(shell_it == this->shellInstance.end())
        {
            // shell is not defined, we are done
            return result;
        }

        c_it = this->bindingEnergy.find(keys[i]);
        if(c_it == this->bindingEnergy.end())
        {
            std::cout << "Shell defined but energy not set " << keys[i] << std::endl;
            throw std::runtime_error("Shell defined but shell energy not set!");
        }
        if (energy <= c_it->second)
        {
            // we cannot excite that shell
            continue;
        }
        fluorescenceRatios = shell_it->second.getFluorescenceRatios();
        for (c_it = fluorescenceRatios.begin(); c_it != fluorescenceRatios.end(); c_it++)
        {
            rate = shell_it->second.getFluorescenceYield();
            if (rate > 0.0)
            {
                result[c_it->first] = this->getTransitionEnergy(c_it->first);
            }
        }
    }
    return result;
}

double Element::getTransitionEnergy(const std::string & transition) const
{
    std::string fromShell, toShell;
    std::map<std::string, Shell>::const_iterator shell_it;
    std::map<std::string, double>::const_iterator bind_it;
    double energy0, energy1;

    if (transition.size() == 4)
    {
        fromShell = transition.substr(transition.size() - 2 , 2);
        toShell = transition.substr(0, 2);
    }
    else
    {
        if (transition.size() == 3)
        {
            fromShell = transition.substr(transition.size() - 2 , 2);
            toShell = transition.substr(0, 1);
        }
        else
        {
            std::cout << "Fluorescence transition " << transition << std::endl;
            throw std::domain_error("Invalid flurescence transition");
        }
    }

    bind_it = this->bindingEnergy.find(toShell);
    if(bind_it == this->bindingEnergy.end())
    {
        std::cout << "Fluorescence transition " << transition << std::endl;
        throw std::domain_error("Transition to an undefined shell!");
    }
    energy0 = bind_it->second;
    if (energy0 <= 0)
    {
        std::cout << "Fluorescence transition " << transition << std::endl;
        throw std::domain_error("Transition to a shell with 0 binding energy!");
    }
    bind_it = this->bindingEnergy.find(fromShell);
    if (bind_it == this->bindingEnergy.end())
    {
        std::cout << "Fluorescence transition from undefined shell ";
        std::cout << fromShell << std::endl;
        energy1 = 0.0;
    }
    else
    {
        energy1 = bind_it->second;
    }
    if(energy1 <= 0.0)
    {
        if (energy1 < 0.0)
        {
            std::cout << this->name << " " << bind_it->first << " " ;
            std::cout << bind_it->second << std::endl;
            throw std::runtime_error("Negative binding energy!");
        }
        else
        {
#ifndef NDEBUG
            if (0)
            {
                std::cout << "Fluorescence transition from unset energy shell ";
                std::cout << " Element = " << this->name;
                std::cout << "Transition = " << transition << std::endl;
                std::cout << fromShell << "Assuming 3 eV" << std::endl;
            }
#endif
            energy1 = 0.003;
        }
    }
    return (energy0 - energy1);
}

std::map<std::string, double> \
Element::getCascadeModifiedVacancyDistribution(const std::map<std::string, double> & distribution) const
{
    std::map<std::string, double> finalDistribution;
    std::map<std::string, double>::const_iterator c_it;
    std::string keys[9] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"};
    std::vector<std::string>::size_type i;
    std::vector<std::string>::size_type j;
    double rate;
    std::string tmpString;
    //bool cascade = true;
    std::map<std::string, double> transferRatios;
    std::map<std::string, double> fluorescenceRatios;
    std::map<std::string, std::map<std::string, double> >result;
    std::map<std::string, Shell>::const_iterator shell_it;
    std::map<std::string, double>::const_iterator bind_it;

    // get a complete initial distribution of vacancies
    for (i = 0; i < this->shellInstance.size(); i++)
    {
        c_it = distribution.find(keys[i]);
        if (c_it != distribution.end())
        {
            rate = c_it->second;
        }
        else
        {
            rate = 0.0;
        }
        finalDistribution[keys[i]] = rate;
    }

    // update the distribution due to the cascade
    for (i = 0; i < finalDistribution.size(); i++)
    {
        if(finalDistribution[keys[i]] > 0.0)
        {
            // we have initial vacancies in shell i
            // propagate to all the higher shells j
            shell_it = this->shellInstance.find(keys[i]);
            for (j = i + 1; j < finalDistribution.size(); j++)
            {
                rate = 0.0;
                transferRatios.clear();
                transferRatios = shell_it->second.getDirectVacancyTransferRatios(keys[j]);
                for (c_it = transferRatios.begin(); c_it != transferRatios.end(); ++c_it)
                {
                    rate += c_it->second;
                }
                finalDistribution[keys[j]] += (rate * finalDistribution[keys[i]]);
            }
        }
    }
    return finalDistribution;
}

std::map<std::string, std::map<std::string, double> >\
Element::getXRayLinesFromVacancyDistribution(const std::map<std::string, double> & distribution, \
                                             const int & cascade, \
                                             const int & useFluorescenceYield) const
{
    std::map<std::string, double>::const_iterator c_it;
    std::string keys[9] = {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"};
    std::vector<std::string>::size_type i;
    double rate, energy0, energy1;
    std::string tmpString;
    std::map<std::string, double> finalDistribution;
    //bool cascade = true;
    std::map<std::string, double> transferRatios;
    std::map<std::string, double> fluorescenceRatios;
    std::map<std::string, std::map<std::string, double> >result;
    std::map<std::string, Shell>::const_iterator shell_it;
    std::map<std::string, double>::const_iterator bind_it;


    if (cascade != 0)
    {
        if (this->cascadeCacheEnabledFlag &&  useFluorescenceYield && (this->cascadeCache.size() > 0))
        {
            // we are in conditions to use the cached emission
            std::map<std::string, std::map<std::string, \
                                std::map<std::string, double> > >::const_iterator cacheKey;
            result.clear();
            for (std::map<std::string, double>::const_iterator c_it = distribution.begin();
                c_it != distribution.end(); ++c_it)
            {
                if ((c_it->second <= 0.0) || (c_it->first == "all other"))
                {
                    // no vacancies created in that shell
                    continue;
                }
                cacheKey = cascadeCache.find(c_it->first);
                if (cacheKey == cascadeCache.end())
                {
                    std::cout << this->name << " Error processing vacancy on shell " << c_it->first << std::endl;
                    throw std::runtime_error("Vacancy in a shell not present in the cache!");
                    // the above is not good, there can be no emission following a vacancy (rate too low)
                    // it seems to be triggered by As L1 vacancies
                }
                for(std::map<std::string, std::map<std::string, double> >::const_iterator rayIt = \
                    cacheKey->second.begin(); rayIt != cacheKey->second.end(); ++rayIt)
                {
                    if (result.find(rayIt->first) == result.end())
                    {
                        result[rayIt->first] = rayIt->second;
                        result[rayIt->first]["rate"] *= c_it->second;
                    }
                    else
                    {
                        std::map<std::string, double>::const_iterator tmpIterator;
                        tmpIterator = rayIt->second.find("rate");
                        result[rayIt->first] ["rate"] += tmpIterator->second * c_it->second;
                    }
                }
            }
            return result;
        }
        finalDistribution = this->getCascadeModifiedVacancyDistribution(distribution);
    }
    else
    {
        // get a complete initial distribution of vacancies
        for (i = 0; i < this->shellInstance.size(); i++)
        {
            c_it = distribution.find(keys[i]);
            if (c_it != distribution.end())
            {
                rate = c_it->second;
            }
            else
            {
                rate = 0.0;
            }
            finalDistribution[keys[i]] = rate;
        }
    }

    // we have the final vacancy distribution accounting for cascade
    // we just have to generate the dictionary of rates and energies
    // for each transition
    for (i = 0; i < finalDistribution.size(); i++)
    {
        // get all the fluorescence transitions to that shell
        shell_it = this->shellInstance.find(keys[i]);
        fluorescenceRatios = shell_it->second.getFluorescenceRatios();
        for (c_it = fluorescenceRatios.begin(); c_it != fluorescenceRatios.end(); c_it++)
        {
            rate = c_it->second * finalDistribution[keys[i]];
            if (useFluorescenceYield != 0)
            {
                rate *= shell_it->second.getFluorescenceYield();
            }
            if (rate > 0.0)
            {
                result[c_it->first]["rate"] = rate;
                bind_it = this->bindingEnergy.find(keys[i]);
                if(bind_it == this->bindingEnergy.end())
                {
                    std::cout << "Fluorescence transition " << c_it->first << std::endl;
                    throw std::domain_error("Transition to an undefined shell!");
                }
                energy0 = bind_it->second;
                if (energy0 <= 0)
                {
                    std::cout << "Fluorescence transition " << c_it->first << std::endl;
                    throw std::domain_error("Transition to a shell with 0 binding energy!");
                }
                tmpString = c_it->first.substr(c_it->first.size() - 2, 2);
                bind_it = this->bindingEnergy.find(tmpString);
                if (bind_it == this->bindingEnergy.end())
                {
                    std::cout << "Fluorescence transition from undefined shell ";
                    std::cout << tmpString << std::endl;
                    energy1 = 0.0;
                }
                else
                {
                    energy1 = bind_it->second;
                }
                if(energy1 <= 0.0)
                {
                    if (energy1 < 0.0)
                    {
                        std::cout << this->name << " " << bind_it->first << " " ;
                        std::cout << bind_it->second << std::endl;
                        throw std::runtime_error("Negative binding energy!");
                    }
                    else
                    {
#ifndef NDEBUG
                        if (0)
                        {
                        std::cout << "Element = " << this->name << " ";
                        std::cout << "Fluorescence transition " << c_it->first << std::endl;
                        std::cout << "rate = " << rate << std::endl;
                        std::cout << "Fluorescence transition from unset energy shell ";
                        std::cout << " destination shell energy = " << energy0 << std::endl;
                        std::cout << tmpString << "Assuming 3 eV" << std::endl;
                        }
#endif
                        energy1 = 0.003;
                    }
                }
                result[c_it->first]["energy"] = energy0 - energy1;
            }
        }
    }
    return result;
}


const Shell & Element::getShell(const std::string & name) const
{
    std::map<std::string, Shell>::const_iterator it;

    it = this->shellInstance.find(name);
    if (it == this->shellInstance.end())
    {
        std::cout << "Undefined shell " << name << std::endl;
        throw std::invalid_argument("Non defined shell: " + name);
    }
    return it->second;
}


void Element::setAtomicMass(const double & A)
{
    // For the time being only positive numbers accepted
    if (A < 0)
    {
        throw std::invalid_argument("Atomic mass should be positive");
    }
    this->atomicMass = A;
}

const double & Element::getAtomicMass() const
{
    // For the time being only positive numbers accepted
    return this->atomicMass;
}


void Element::setAtomicNumber(const int & z)
{
    // For the time being only positive numbers accepted
    if (z < 1)
    {
        throw std::invalid_argument("Atomic number should be positive");
    }
    this->atomicNumber = z;
}

const int & Element::getAtomicNumber() const
{
    // For the time being only positive numbers accepted
    return this->atomicNumber;
}


std::vector<std::map<std::string, std::map<std::string, double> > >Element::getPhotoelectricExcitationFactors( \
                            const std::vector<double> & energy,
                            const std::vector<double> & weights) const
{
    double weight;
    std::vector<double>::size_type i;
    std::vector<std::map<std::string, std::map<std::string, double> > > result;
    std::map<std::string, std::map<std::string, double> >::iterator it;

    if (energy.size() == 0)
    {
        result.clear();
        return result;
    }
    if (weights.size() == 1)
        weight = weights[0];
    else
        weight = 1.0 / energy.size();

    result.resize(energy.size());
    for(i = 0; i < energy.size(); i++)
    {
        if (weights.size() > 1)
            weight = weights[i];
        result[i] = this->getPhotoelectricExcitationFactors(energy[i], weight);
    }
    return result;
}

std::map<std::string, std::map<std::string, double> > Element::getPhotoelectricExcitationFactors( \
                                                    const double & energy,
                                                    const double & weight) const
{
    std::map<std::string, double>vacancyDistribution;
    std::map<std::string, std::map<std::string, double> > result;
    std::map<std::string, std::map<std::string, double> >::iterator it;
    std::map<double, std::map<std::string, std::map<std::string, double> > >::const_iterator c_it;
    result.clear();

    if (this->isCacheEnabled())
    {
        if (this->excitationFactorsCache.size())
        {
            c_it = this->excitationFactorsCache.find(energy);
            if (c_it != this->excitationFactorsCache.end())
            {
                result = c_it->second;
                for(it = result.begin(); it != result.end(); ++it)
                {
                    it->second["factor"] = it->second["factor"] * weight;
                    it->second["rate"] = it->second["rate"] * weight;
                }
                return result;
            }
        }
    }
    //we have to calculate it
    vacancyDistribution = this->getInitialPhotoelectricVacancyDistribution(energy);
    result = this->getXRayLinesFromVacancyDistribution(vacancyDistribution, 1, 1);
    for(it = result.begin(); it != result.end(); ++it)
    {
        it->second["factor"] = it->second["rate"] * weight;
        it->second["rate"] = it->second["factor"] *\
                        this->getMassAttenuationCoefficients(energy)["photoelectric"];
    }
    return result;
}


std::pair<long, long> Element::getInterpolationIndices(const std::vector<double> & vec, const double & x) const
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

void Element::setCascadeCacheEnabled(const int & flag)
{
    if (flag == 0)
    {
        this->cascadeCacheEnabledFlag = false;
    }
    else
    {
        if (this->cascadeCache.size() < 1)
        {
            this->fillCascadeCache();
        }
        this->cascadeCacheEnabledFlag = true;
    }
}

void Element::fillCascadeCache()
{
    std::map<std::string, Shell>::const_iterator shellIterator;
    bool oldFlag;
    this->cascadeCache.clear();
    // std::cout << "filling cache for element " << this->name << std::endl;
    for (shellIterator = shellInstance.begin(); shellIterator != shellInstance.end(); ++shellIterator)
    {
        std::string subshell = shellIterator->first;
        std::map<std::string, double> tmpDistribution;
        tmpDistribution.clear();
        tmpDistribution[subshell] = 1.0;
        oldFlag = this->cascadeCacheEnabledFlag;
        if (oldFlag)
        {
            // this is needed because otherways it uses the cache when filling the cache :-)
            this->cascadeCacheEnabledFlag = false;
        }
        this->cascadeCache[subshell] = this->getXRayLinesFromVacancyDistribution(tmpDistribution, 1, 1);
        this->cascadeCacheEnabledFlag = oldFlag;
    }
    //std::cout << "cache filled for element " << this->name << std::endl;
}

void Element::emptyCascadeCache()
{
    this->cascadeCache.clear();
}

int Element::isCascadeCacheFilled() const
{
    return this->cascadeCache.size();
}

void Element::setCacheEnabled(const int & flag)
{
    if (flag == 0)
        this->calculationCacheEnabledFlag = false;
    else
        this->calculationCacheEnabledFlag = true;
}

const int Element::isCacheEnabled() const
{
    if (this->calculationCacheEnabledFlag)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

void Element::clearCache()
{
    this->excitationFactorsCache.clear();
    this->muCache.clear();
}

void Element::fillCache(const std::vector<double> & energy)
{
    std::vector<double>::size_type i, maxSize;
    int status;

    this->clearCache();

    if (energy.size() < this->cacheMaximumSize)
    {
        maxSize = energy.size();
    }
    else
    {
        maxSize = this->cacheMaximumSize;
    }
    status = this->isCacheEnabled();
    try
    {
        // This is necessary because clang pushes return value optimization to the limit
        // and creates the key into muCache prior to call getMassAttenuationCoefficients
        this->setCacheEnabled(0);
        for (i = 0; i < maxSize; i++)
        {
            this->muCache[energy[i]] = this->getMassAttenuationCoefficients(energy[i]);
            this->excitationFactorsCache[energy[i]] = this->getPhotoelectricExcitationFactors(energy[i], 1.0);
        }
        this->setCacheEnabled(status);
    }
    catch(...)
    {
        this->setCacheEnabled(status);
        throw;
    }
}

void Element::updateCache(const std::vector< double> & energy)
{
    std::vector<double>::size_type i, eSize;
    int status;

    // This is necessary because clang pushes return value optimization to the limit
    // and creates the key into muCache prior to call getMassAttenuationCoefficients
    status = this->isCacheEnabled();
    this->setCacheEnabled(0);
    eSize = energy.size();
    try
    {
        for (i = 0; i < eSize; i++)
        {
            if (this->muCache.size() < this->cacheMaximumSize)
            {
                if (this->muCache.find(energy[i]) == this->muCache.end())
                {
                    this->muCache[energy[i]] = this->getMassAttenuationCoefficients(energy[i]);
                }
                if (this->excitationFactorsCache.find(energy[i]) == this->excitationFactorsCache.end())
                {
                    this->excitationFactorsCache[energy[i]] = this->getPhotoelectricExcitationFactors(energy[i], 1.0);
                }
            }
        }
        this->setCacheEnabled(status);
    }
    catch(...)
    {
        this->setCacheEnabled(status);
        throw;
    }
    if (this->muCache.size() >= this->cacheMaximumSize)
    {
        std::cout << "Mass attenuation coefficients cache full" << std::endl;
    }
    if (this->excitationFactorsCache.size() >= this->cacheMaximumSize)
    {
        std::cout << "Excitation factors cache full" << std::endl;
    }
}

int Element::getCacheSize() const
{
    return this->muCache.size();
}

} // namespace fisx
