#include "fisx_beam.h"
#include <algorithm>

Beam::Beam()
{
    this->normalized = false;
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
            this->rays[i].characteristic = 1.0;
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

std::vector<double> Beam::getBeam()
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
