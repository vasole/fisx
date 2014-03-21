#include "fisx_beam.h"
#include <algorithm>

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
    this->rays.resize(energy.size());

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



void Beam::setBeam(int nValues, double *energy, double *weight,
                   int *characteristic, double *divergency)
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

void Beam::setBeam(int nValues, double *energy, double *weight,
                   double *characteristic, double *divergency)
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
            this->rays[i].characteristic = (int) characteristic[i];
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

std::vector<double> Beam::getBeamAsDoubleVector()
{
    std::vector<double>::size_type nItems;
    std::vector<Ray>::const_iterator c_it;
    std::vector<double> returnValue;
    Ray ray;

    if (!this->normalized)
    {
        this->normalizeBeam();
    }
    nItems = 4 * this->rays.size();
    returnValue.resize(nItems);
    nItems = 0;
    for(c_it = this->rays.begin(); c_it != this->rays.end(); ++c_it)
    {
        ray = *c_it;
        returnValue[nItems] = ray.energy;
        nItems++;
        returnValue[nItems] = ray.weight;
        nItems++;
        returnValue[nItems] = ray.characteristic;
        nItems++;
        returnValue[nItems] = ray.divergency;
        nItems++;
    }
    return returnValue;
}


const std::vector<Ray> & Beam::getBeam()
{
    if (!this->normalized)
        this->normalizeBeam();
    return this->rays;
}
