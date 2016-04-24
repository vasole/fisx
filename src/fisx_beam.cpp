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
#include "fisx_beam.h"
#include <algorithm>

namespace fisx
{

Beam::Beam()
{
    this->normalized = false;
}

void Beam::setBeam(const std::vector<double> & energy, \
                   const std::vector<double> & weight,\
                   const std::vector<int> & characteristic,\
                   const std::vector<double> & divergency)
{
    std::vector<double>::size_type i;
    std::vector<double>::size_type j;
    double defaultWeight;
    int defaultCharacteristic;
    double defaultDivergency;

    this->normalized = false;

    if(energy.size() > 0)
    {
        this->rays.resize(energy.size());
    }
    else
    {
        this->rays.clear();
        return;
    }

    defaultWeight = 1.0;
    if (weight.size())
    {
        defaultWeight = weight[0];
    }

    defaultCharacteristic = 1;
    if (characteristic.size())
    {
        defaultCharacteristic = characteristic[0];
    }

    defaultDivergency = 0.0;
    if (divergency.size())
    {
        defaultDivergency = divergency[0];
    }

    for(i = 0; i < this->rays.size(); i++)
    {
        this->rays[i].energy = energy[i];

        // weight is optional
        j = weight.size();
        if (j > 1)
        {
            this->rays[i].weight = weight[i];
        }
        else
        {
            this->rays[i].weight = defaultWeight;
        }

        // characteristic is optional
        j = characteristic.size();
        if (j > 1)
        {
            this->rays[i].characteristic = characteristic[i];
        }
        else
        {
            this->rays[i].characteristic = defaultCharacteristic;
        }

        // divergency is optional
        j = divergency.size();
        if (j > 1)
        {
            this->rays[i].divergency = divergency[i];
        }
        else
        {
            this->rays[i].divergency = defaultDivergency;
        }

    }
    this->normalizeBeam();
}

void Beam::setBeam(const double & energy, const double divergency)
{
    this->normalized = false;
    this->rays.clear();
    this->rays.resize(1);
    this->rays[0].energy = energy;
    this->rays[0].weight = 1.0;
    this->rays[0].characteristic = 1;
    this->rays[0].divergency = divergency;
    // it is already normalized
    this->normalizeBeam();
}

void Beam::setBeam(const int & nValues, const double *energy, const double *weight,
                 const int *characteristic, const double *divergency)
{
    int i;
    double tmpDouble;

    this->normalized = false;
    this->rays.clear();
    this->rays.resize(nValues);

    tmpDouble = 1.0;

    for (i=0; i < nValues; ++i)
    {
        this->rays[i].energy = energy[i];
        if (weight != NULL)
        {
            tmpDouble = weight[i];
        }
        this->rays[i].weight = tmpDouble;
        if (characteristic == NULL)
        {
            this->rays[i].characteristic = 1;
        }
        else
        {
            this->rays[i].characteristic = characteristic[i];
        }
        if (divergency == NULL)
        {
            this->rays[i].divergency = 0.0;
        }
        else
        {
            this->rays[i].divergency = divergency[i];
        }
    }
    this->normalizeBeam();
}

void Beam::setBeam(const int & nValues, const double *energy, const double *weight,
                   const double *characteristic, const double *divergency)
{
    int i;
    std::vector<int> flags;

    flags.resize(nValues);
    for (i = 0; i < nValues; i++)
    {
        flags[i] = (int) characteristic[i];
    }

    this->setBeam(nValues, energy, weight, &flags[0], divergency);
}

void Beam::normalizeBeam()
{
    std::vector<Ray>::size_type nValues;
    std::vector<Ray>::size_type i;
    double totalWeight;

    nValues = this->rays.size();
    totalWeight = 0.0;

    for (i = 0; i < nValues; ++i)
    {
        totalWeight += this->rays[i].weight;
    }
    if (totalWeight > 0.0)
    {
        for (i = 0; i < nValues; ++i)
        {
            this->rays[i].weight /= totalWeight;
        }
    }
    this->normalized = true;
    std::sort(this->rays.begin(), this->rays.end());
}

std::vector<std::vector<double> > Beam::getBeamAsDoubleVectors() const
{
    std::vector<double>::size_type nItems;
    std::vector<Ray>::size_type c_it;
    std::vector<std::vector<double> >returnValue;
    const Ray *ray;

    //if (!this->normalized)
    //{
    //    this->normalizeBeam();
    //}
    nItems = this->rays.size();
    returnValue.resize(4);
    if (nItems > 0)
    {
        returnValue[0].resize(nItems);
        returnValue[1].resize(nItems);
        returnValue[2].resize(nItems);
        returnValue[3].resize(nItems);
        for(c_it = 0; c_it < nItems; c_it++)
        {
            ray = &(this->rays[c_it]);
            returnValue[0][c_it] = (*ray).energy;
            returnValue[1][c_it] = (*ray).weight;
            returnValue[2][c_it] = (*ray).characteristic;
            returnValue[3][c_it] = (*ray).divergency;
        }
    }
    return returnValue;
}


const std::vector<Ray> & Beam::getBeam()
{
    //if (!this->normalized)
    //    this->normalizeBeam();
    return this->rays;
}

std::ostream& operator<< (std::ostream& o, Beam const & beam)
{
    std::vector<Ray>::size_type i;
    for(i = 0; i < beam.rays.size(); i++)
    {
        o << "E (keV) = " << beam.rays[i].energy << " weight = " << beam.rays[i].weight;
        if (i != (beam.rays.size() - 1)) o << std::endl;
    }
    return o;
}

} // namespace fisx
