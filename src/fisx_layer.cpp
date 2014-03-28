#include "fisx_layer.h"
#include <math.h>
#include <stdexcept>

Layer::Layer(const std::string & name, const double & density, \
             const double & thickness, const double & funnyFactor)
{
    this->name = name;
    this->materialName = name;
    this->density = density;
    this->thickness = thickness;
    this->funnyFactor = funnyFactor;
    this->hasMaterial = false;
}

void Layer::setMaterial(const std::string & materialName)
{
    this->materialName = materialName;
    this->hasMaterial = false;
}

void Layer::setMaterial(const Material & material)
{
    this->material = material;
    if (this->density < 0.0)
    {
        this->density = this->material.getDefaultDensity();
    }
    if (this->thickness <= 0.0)
    {
        this->thickness = this->material.getDefaultThickness();
    }
    this->hasMaterial = true;
}

bool Layer::hasMaterialComposition() const
{
    return this->hasMaterial;
}

std::map<std::string, double> Layer::getMassAttenuationCoefficients(const double & energy, \
                                                const Elements & elements) const
{
    if (this->hasMaterial)
    {
        return elements.getMassAttenuationCoefficients(this->material.getComposition(), energy);
    }
    else
    {
        return elements.getMassAttenuationCoefficients(this->materialName, energy);
    }
}

double Layer::getTransmission(const double & energy, const Elements & elements, const double & angle) const
{
    // The material might not have been defined in the  current Elements instance.
    // However, its composition might be fine.
    const double PI = acos(-1.0);
    double muTotal;
    double tmpDouble;
    if (this->hasMaterial)
    {
        muTotal = elements.getMassAttenuationCoefficients(this->material.getComposition(), energy)["total"];
    }
    else
    {
        muTotal = elements.getMassAttenuationCoefficients(this->materialName, energy)["total"];
    }
    if (angle == 90.0)
    {
        tmpDouble = this->density * this->thickness;
    }
    else
    {
        if (angle < 0)
            tmpDouble = sin(-angle * PI / 180.);
        else
            tmpDouble = sin(-angle * PI / 180.);
        tmpDouble /= this->density * this->thickness;
    }
    if(tmpDouble <= 0.0)
    {
        std::string msg;
        msg = "Layer " + this->name + " thickness is " + elements.toString(tmpDouble) + " g/cm2";
        throw std::runtime_error( msg );
    }

    return (1.0 - this->funnyFactor) + (this->funnyFactor * exp(-(tmpDouble * muTotal)));
}

std::vector<double> Layer::getTransmission(const std::vector<double> & energy, const Elements & elements, \
                                           const double & angle) const
{
    const double PI = acos(-1.0);
    std::vector<double>::size_type i;
    std::vector<double> tmpDoubleVector;
    double tmpDouble;

    if (angle == 90.0)
    {
        tmpDouble = this->density * this->thickness;
    }
    else
    {
        if (angle < 0)
            tmpDouble = sin(-angle * PI / 180.);
        else
            tmpDouble = sin(-angle * PI / 180.);
        tmpDouble /= this->density * this->thickness;
    }

    if(tmpDouble <= 0.0)
    {
        std::string msg;
        msg = "Layer " + this->name + " thickness is " + elements.toString(tmpDouble) + " g/cm2";
        throw std::runtime_error( msg );
    }

    if (this->hasMaterial)
    {
        tmpDoubleVector = elements.getMassAttenuationCoefficients(this->material.getComposition(), energy)["total"];
    }
    else
    {
        tmpDoubleVector = elements.getMassAttenuationCoefficients(this->materialName, energy)["total"];
    }
    for (i = 0; i < tmpDoubleVector.size(); i++)
    {
        tmpDoubleVector[i] = (1.0 - this->funnyFactor) + \
                              (this->funnyFactor * exp(-(tmpDouble * tmpDoubleVector[i])));
    }
    return tmpDoubleVector;
}

