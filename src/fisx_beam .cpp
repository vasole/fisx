#include "fisx_beam.h"
#include <algorithm>
#include <iostream>

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
    std::vector<int>::size_type i;
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
    nItems = rays.size();
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
