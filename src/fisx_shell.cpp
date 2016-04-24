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
#include <ctype.h>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cstdlib> // needed for atoi
#include "fisx_shell.h"
// #include <iostream>

namespace fisx
{

Shell::Shell()
{
    this->shellConstants["omega"] = 0.0;
}

Shell::Shell(std::string name)
{
    int idx = 0;
    int i, maxShellSubindex;
    int shellMainIndex = -1;
    std::string msg;
    std::string fij;
    std::string digit[10]={"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};

    if (name.size() > 2)
    {
        throw std::invalid_argument("Invalid shell name");
    }
    if (name.substr(0, 1) == "K")
    {
        if (name.size() == 1)
        {
            shellMainIndex = 0;
        }
    }
    if (name.substr(0, 1) == "L")
    {
        if (name.size() == 2)
        {
            shellMainIndex = 1;
        }
    }
    if (name.substr(0, 1) == "M")
    {
        if (name.size() == 2)
        {
            shellMainIndex = 2;
        }
    }

    if(shellMainIndex < 0)
    {
        throw std::invalid_argument("Only shells K, L1, L2, L3, M1, M2, M3, M4 or M5 accepted");
    }

    // omega = 0.0;
    // bindingEnergy = 0.0;
    // ck = NULL;

    if (name.size() == 2)
    {
        if (!this->StringToInteger(name.substr(1, 1), idx))
        {
            msg = "Subshell index must a number but got <" + name.substr(1, 1) + ">";
            throw std::invalid_argument(msg);
        }
        if (idx < 1)
        {
            throw std::invalid_argument("Invalid subshell index");
        }
        if ((shellMainIndex == 1) && (idx > 3))
        {
            std::cout << "INDEX = " << idx << " Obtained from " << name.substr(1, 1) << std::endl;
            throw std::invalid_argument("Incompatible L subshell index");
        }
        if ((shellMainIndex == 2) && (idx > 5))
        {
            throw std::invalid_argument("Incompatible M subshell index");
        }
    }
    this->name = name;
    this->shellMainIndex = shellMainIndex;
    this->subshellIndex = idx;
    this->shellConstants["omega"] = 0.0;

    if (shellMainIndex > 0)
    {
        if (shellMainIndex == 1)
        {
            // Maximum index for L shell is 3
            maxShellSubindex = 3;
        }
        else
        {
            // Maximum index for L shell is 5
            maxShellSubindex = 5;
        }
        for(i = idx + 1; i <= maxShellSubindex; i++)
        {
            fij = "f" + this->name.substr(1, 1) + digit[i];
            this->shellConstants[fij] = 0.0;
        }
    }
}

std::string Shell::toUpperCaseString(const std::string& str)
{
    std::string::size_type i;
    std::string converted;
    for(i = 0; i < str.size(); ++i)
        converted += toupper(str[i]);
    return converted;
}


void Shell::setNonradiativeTransitions(std::vector<std::string> labels, std::vector<double> values)
{
    std::vector<int>::size_type i;

    if (this->nonradiativeTransitions.size() > 0)
    {
        // empty current list
        this->nonradiativeTransitions.clear();
    }
    for (i = 0; i != labels.size(); i++)
    {
        this->nonradiativeTransitions[this->toUpperCaseString(labels[i])] = values[i];
    }
    this->_updateNonradiativeRatios();
}

void Shell::setNonradiativeTransitions(std::map<std::string, double> inputMap)
{
    std::vector<std::string> transitions;
    std::vector<double> values;
    std::vector<double>::size_type iTransitions;
    std::vector<double>::size_type iValues;
    std::map<std::string, double>::const_iterator it;

    transitions.resize(inputMap.size());
    values.resize(inputMap.size());

    iTransitions = 0;
    iValues = 0;
    for(it = inputMap.begin(); it != inputMap.end(); ++it)
    {
        transitions[iTransitions] = it->first;
        values[iValues] = it->second;
        iTransitions++;
        iValues++;
    }
    this->setNonradiativeTransitions(transitions, values);
}

void Shell::setNonradiativeTransitions(const char *c_strings[], const double *values, int nValues)
{
    int i;
    for(i=0; i < nValues; i++)
    {
        this->nonradiativeTransitions[this->toUpperCaseString(std::string(c_strings[i]))] = values[i];
    }
    this->_updateNonradiativeRatios();
}

void Shell::setRadiativeTransitions(std::vector<std::string> labels, std::vector<double> values)
{
    std::string tmpString;
    std::vector<int>::size_type i;
    if (this->radiativeTransitions.size() > 0)
    {
        // empty current list
        this->radiativeTransitions.clear();
    }
    for (i = 0; i != labels.size(); i++)
    {
        tmpString = this->toUpperCaseString(labels[i]);
        if (tmpString.size() == 2)
        {
            // some sets of radiative transitions write KO and KP transitions
            // in those cases PyMca was calculating the energies from KO2 and KP2
            // therefore change the label accordingly
            tmpString += "2";
        }
        if ((tmpString.size() == 5) && (tmpString != "TOTAL"))
        {
            if (isdigit(tmpString.at(3)) && isdigit(tmpString.at(4)))
            {
                // These are the following cases in PyMca:
                // L1P23, L1O45
                // L2P23,
                // L3O45, L3P23, L3P45
                // one should not supply this type of transitions
                // because it will not be possible to derive the associated
                // energy will not be derived from the binding energies
                // PyMca was deriving the energy skipping the last digit.
                // std::cout << "setting key" << tmpString.substr(0, 4) << std::endl;
                this->radiativeTransitions[tmpString.substr(0, 4)] = values[i];
                // other alternative would be to generate two keys with half
                // the rate as being closer to the real stuff. Nevertheless,
                // both transitions would most likely be mixed due to the
                // small difference in energy between them.
            }
            else
            {
                std::cout << "Not a valid transition <" << tmpString << ">" << std::endl;
                tmpString += " is not a valid transition";
                throw std::invalid_argument(tmpString);
            }
        }
        else
        {
            if (tmpString.substr(0, this->name.size()) == this->name)
            {
                this->radiativeTransitions[tmpString] = values[i];
            }
            else
            {
                if(tmpString == "TOTAL")
                {
                    this->radiativeTransitions[tmpString] = values[i];
                }
                else
                {
                    tmpString += " is not a valid transition for shell " + this->name;
                    throw std::invalid_argument(tmpString);
                }
            }
        }
    }
    this->_updateFluorescenceRatios();
}

void Shell::setRadiativeTransitions(std::map<std::string, double> inputMap)
{
    std::vector<std::string> transitions;
    std::vector<double> values;
    std::vector<double>::size_type iTransitions;
    std::vector<double>::size_type iValues;
    std::map<std::string, double>::const_iterator it;

    transitions.resize(inputMap.size());
    values.resize(inputMap.size());

    iTransitions = 0;
    iValues = 0;
    for(it = inputMap.begin(); it != inputMap.end(); ++it)
    {
        transitions[iTransitions] = it->first;
        values[iValues] = it->second;
        iTransitions++;
        iValues++;
    }
    this->setRadiativeTransitions(transitions, values);
}

void Shell::setRadiativeTransitions(const char *c_strings[], const double *values, int nValues)
{
    std::vector<std::string> labels;
    std::vector<double> vValues;
    int i;

    labels.resize(nValues);
    vValues.resize(nValues);

    for(i=0; i < nValues; i++)
    {
        labels[i] = this->toUpperCaseString(std::string(c_strings[i]));
        vValues[i] = values[i];
    }
    this->setRadiativeTransitions(labels, vValues);
}

const std::map<std::string, double> & Shell::getRadiativeTransitions() const
{
    return this->radiativeTransitions;
}

const std::map<std::string, double> & Shell::getNonradiativeTransitions() const
{
    return this->nonradiativeTransitions;
}

const std::map<std::string, double> & Shell::getAugerRatios() const
{
    return this->augerRatios;
}

const std::map<std::string, std::map<std::string, double> > & Shell::getCosterKronigRatios() const
{
    return this->costerKronigRatios;
}
const std::map<std::string, double> & Shell::getFluorescenceRatios() const
{
    return this->fluorescenceRatios;
}


void Shell::_updateFluorescenceRatios()
{
    double total;
    std::string totalLabel = "TOTAL";
    std::map<std::string, double>::iterator iter;

    if (this->fluorescenceRatios.size() > 0)
    {
        this->fluorescenceRatios.clear();
    }

    // calculate the sum of rates
    // Since usually all the rates are supplied, we ignore any possible total provided because
    // we cannot figure out if we deal with level widths, ratios or probabilities.
    total = 0.0 ;
    for (iter = this->radiativeTransitions.begin(); iter != this->radiativeTransitions.end(); ++iter)
    {
        if (iter->first != totalLabel)
        {
            total += iter->second;
        }
    }

    for (iter = this->radiativeTransitions.begin(); iter != this->radiativeTransitions.end(); ++iter)
    {
        // do not calculate unnecessary ratios
        if (iter->second > 0.0)
        {
            if (iter->first != totalLabel)
            {
                this->fluorescenceRatios[iter->first] = iter->second / total;
            }
        }
    }
}

double Shell::getFluorescenceYield() const
{
    std::map<std::string, double>::const_iterator c_it;
    c_it = this->shellConstants.find("omega");
    return c_it->second;
}

void Shell::_updateNonradiativeRatios()
{
    double total;
    double totalAuger;
    double totalCosterKronig;
    double tmpTotal;
    std::string digit[10]={"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
    int i, j;
    int subshellMax;
    std::string totalLabel = "TOTAL";
    std::map<std::string, double>::iterator iter;
    std::string testString;
    std::string fij;


    if (this->augerRatios.size() > 0)
    {
        this->augerRatios.clear();
    }

    if (this->costerKronigRatios.size() > 0)
    {
        this->costerKronigRatios.clear();
    }

    // calculate sum of rates
    total = 0.0 ;
    for (iter = this->nonradiativeTransitions.begin(); iter != this->nonradiativeTransitions.end(); ++iter)
    {
        if (iter->first != totalLabel)
        {
            total += iter->second;
        }
    }

    // since usually not all the transitions are in the file, but only those affecting
    // the L and M shells either as origin or as shell with intervention in the de-excitation
    // process, we have to retrieve the total from the supplied data (if present)
    if (this->nonradiativeTransitions.find(totalLabel) == this->nonradiativeTransitions.end())
    {
        // The supplied data did not include the total, so, we cannot do anything better than
        // to use the sum of the supplied transitions.
        this->nonradiativeTransitions[totalLabel] = total;
    }

    // now we have to separate Auger and different Coster-Kronig transitions
    if (this->name == std::string("K") || \
        this->name == std::string("L3") || \
        this->name == std::string("M5"))
    {
        // no Coster-Kronig transitions
        totalCosterKronig = 0.0;
        subshellMax = 0;
        totalAuger = total;
    }
    else
    {
        totalCosterKronig = 0.0;
        totalAuger = 0.0;
        if (this->shellMainIndex == 1)
        {
            // L1 or L2 shells
            subshellMax = 3;
        }
        else
        {
            // M1, M2, M3 or M4 shells
            subshellMax = 5;
        }

        for (i = this->subshellIndex + 1 ; i < (subshellMax + 1) ; i++)
        {
            testString = this->name;
            testString.append("-");
            testString.append(this->name.substr(0, 1));
            testString.append(digit[i]);
            fij = "f" + this->name.substr(1, 1) + digit[i];
            tmpTotal = 0.0;
            for (iter = this->nonradiativeTransitions.begin(); iter != this->nonradiativeTransitions.end(); ++iter)
            {
                if (iter->second > 0.0)
                {
                    // startswith like function
                    if (iter->first.substr(0, testString.size()) == testString)
                    {
                        // store the pair (key, value)
                        this->costerKronigRatios[fij][iter->first] = iter->second;
                        tmpTotal += iter->second;
                    }
                }
            }
            if (tmpTotal > 0.0)
            {
                for (iter = this->costerKronigRatios[fij].begin();\
                     iter != this->costerKronigRatios[fij].end(); ++iter)
                {
                    this->costerKronigRatios[fij][iter->first] /= tmpTotal;
                }
            }
            // store the total (mainly for debugging purposes)
            this->costerKronigRatios[fij][totalLabel] = tmpTotal;
            totalCosterKronig += tmpTotal;
        }
    }

    // we have the different Coster-Kronig ratios
    // now calculate the Auger ratios
    // we already have the total of the nonradiative transitions
    // and the total of the costerKronig transitions
    totalAuger = this->nonradiativeTransitions[totalLabel] - totalCosterKronig;

    for (iter = this->nonradiativeTransitions.begin(); iter != this->nonradiativeTransitions.end(); ++iter)
    {
        if (iter->second > 0.0)
        {
            // a flag to go out of the loop
            i = 0;
            if (totalCosterKronig > 0.0)
            {
                // an index to iterate all the Coster-Kronig groups
                j = 0;
                while ((i == 0) && j < (subshellMax - this->subshellIndex))
                {
                    // the corresponding fij transition
                    fij = "f" + this->name.substr(1, 1) + digit[(j + this->subshellIndex + 1)];
                    if (this->costerKronigRatios[fij].find(iter->first) != \
                        (this->costerKronigRatios[fij].end()))
                    {
                        // it is a key belonging to a Coster-Kronig transition
                        i = 1;
                    }
                    j++;
                }
            }
            if (i == 0)
            {
                if (iter->first != totalLabel)
                {
                    // it is an Auger transition
                    this->augerRatios[iter->first] = iter->second / totalAuger;
                }
            }
        }
    }
    // store the total auger (for debugging purposes)
    this->augerRatios[totalLabel] = totalAuger;
}

std::map<std::string, double> Shell::getDirectVacancyTransferRatios(const std::string &destination) const
// Return the probabilities of direct transfer of a vacancy to a higher shell following
// an X-ray emission, an Auger transition and Coster-Kronig transitions (if any).
// It multiplies by the fluorescence, auger or Coster-Kronig yields to get the probabilities.
{
    std::map<std::string, double> output;
    std::map<std::string, double>::const_iterator it;
    std::map<std::string, std::map<std::string, double> >::const_iterator ck_it;
    std::string transition;
    double tmpDouble;

    it = this->shellConstants.find("omega");
    if (it->second == 0.0)
    {
        // TODO: Put a less exigent test, omega can be zero for a lot of shells ...
        // std::cout << "WARNING: Probably using unitialized shell constants!!!! " << std::endl;
        ;
    }

    // destination needs two characters as in L2, M3, ...
    if (destination.size() != 2)
    {
        throw std::invalid_argument("Invalid destination subshell. Two characters needed");
    }

    output["fluorescence"] = 0.0;
    output["auger"] = 0.0;

    // build fluorescence transition key
    transition = this->name + destination;
    it = this->fluorescenceRatios.find(transition);
    if (it != this->fluorescenceRatios.end())
    {
        output["fluorescence"] = it->second;
    }

    // Auger pattern Origin-DestinationAny
    transition = this->name + "-" + destination;
    for (it = this->augerRatios.begin(); it != this->augerRatios.end(); ++it)
    {
        if (it->first.substr(0, transition.size()) == transition)
        {
            output["auger"] += it->second;
        }
    }
    // Auger pattern Origin-AnyDestination
    for (it = this->augerRatios.begin(); it != this->augerRatios.end(); ++it)
    {
        if (std::equal(destination.rbegin(), destination.rend(), it->first.rbegin()))
        {
            output["auger"] += it->second;
        }
    }

    // Coster-Kronig
    for (ck_it = this->costerKronigRatios.begin();\
         ck_it != this->costerKronigRatios.end(); ++ck_it)
    {
        output[ck_it->first] = 0.0;
        // Coster-Kronig pattern Origin-DestinationAny
        transition = this->name + "-" + destination;
        for (it = ck_it->second.begin(); it != ck_it->second.end(); ++it)
        {
            if (it->first.substr(0, transition.size()) == transition)
            {
                output[ck_it->first] += it->second;
            }
        }
    }
    for (ck_it = this->costerKronigRatios.begin();\
         ck_it != this->costerKronigRatios.end(); ++ck_it)
    {
        // Coster-Kronig pattern Origin-AnyDestination
        for (it = ck_it->second.begin(); it != ck_it->second.end(); ++it)
        {
            if (std::equal(destination.rbegin(), destination.rend(), it->first.rbegin()))
            {
                output[ck_it->first] += it->second;
            }
        }
    }

    // multiply by the respective shell constants

    tmpDouble = 0.0;
    for (it = this->shellConstants.begin(); it != this->shellConstants.end(); ++it)
    {
        tmpDouble += it->second;
    }
    if (tmpDouble > 1.0001)
    {
        for (it = this->shellConstants.begin(); it != this->shellConstants.end(); ++it)
        {
            std::cout << it->first << " = " << it->second << std::endl;
        }
        throw std::domain_error("Sum of CosterKronig and Fluorescence yields greater than 1.0");
    }

    std::map<std::string, double>::const_iterator const_it;
    for (it = output.begin(); it != output.end(); ++it)
    {
        if (it->first == "auger")
        {
            output[it->first] *= (1.0 - tmpDouble);
        }
        else
        {
            if (it->first == "fluorescence")
            {
                const_it = this->shellConstants.find("omega");
            }
            else
            {
                const_it = this->shellConstants.find(it->first);
            }
            output[it->first] *= const_it->second;
        }
    }
    return output;
}

// WARNING: The original constants are not cleared
//          Only those constants supplied will be overwritten!!!
void Shell::setShellConstants(std::map<std::string, double> shellConstants)
{
    std::string msg;
    std::map<std::string, double>::const_iterator c_it;

    for (c_it = shellConstants.begin(); c_it != shellConstants.end(); ++c_it)
    {
        if(this->shellConstants.find(c_it->first) != this->shellConstants.end())
        {
            this->shellConstants[c_it->first] = c_it->second;
        }
        else
        {
            msg = "Invalid constant " + c_it->first + " for " + this->name +" shell";
            throw std::invalid_argument(msg);
        }
    }
}

const std::map<std::string, double> & Shell::getShellConstants() const
{
    return this->shellConstants;

}

bool Shell::StringToInteger(const std::string& str, int & number)
{

    std::istringstream i(str);

    if (!(i >> number))
    {
        // Number conversion failed
        return false;
    }
    return true;
}

} // namespace fisx
