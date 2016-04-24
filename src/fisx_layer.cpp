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
#include "fisx_layer.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace fisx
{

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

void Layer::setDensity(const double & value)
{
    if (value <= 0.0)
    {
        throw std::invalid_argument("Density must be possitive value.");
    }
    this->density = value;
}

void Layer::setThickness(const double & value)
{
    if (value <= 0.0)
    {
        throw std::invalid_argument("Thickness must be possitive value.");
    }
    this->thickness = value;
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

std::map<std::string, double> Layer::getComposition(const Elements & elements)
{
    if (this->hasMaterial)
    {
        return this->material.getComposition();
    }
    else
    {
        return elements.getComposition(this->materialName);
    }
}

std::map<std::string, std::vector<double> > Layer::getMassAttenuationCoefficients( \
                                                   const std::vector<double> & energy, \
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

std::map<std::string, double> Layer::getComposition(const Elements & elements) const
{
    if (this->hasMaterial)
    {
        return this->material.getComposition();
    }
    else
    {
        return elements.getComposition(this->materialName);
    }
}

double Layer::getTransmission(const double & energy, const Elements & elements, const double & angle) const
{
    // The material might not have been defined in the  current Elements instance.
    // However, its composition might be fine.
    const double PI = std::acos(-1.0);
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
        {
            tmpDouble = std::sin(-angle * PI / 180.);
        }
        else
        {
            tmpDouble = std::sin(angle * PI / 180.);
        }
        tmpDouble = (this->density * this->thickness) / tmpDouble;
    }
    if(tmpDouble <= 0.0)
    {
        std::string msg;
        msg = "Layer " + this->name + " thickness is " + elements.toString(tmpDouble) + " g/cm2";
        throw std::runtime_error( msg );
    }

    return (1.0 - this->funnyFactor) + (this->funnyFactor * std::exp(-(tmpDouble * muTotal)));
}

std::vector<double> Layer::getTransmission(const std::vector<double> & energy, const Elements & elements, \
                                           const double & angle) const
{
    const double PI = std::acos(-1.0);
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
            tmpDouble = std::sin(-angle * PI / 180.);
        else
            tmpDouble = std::sin(angle * PI / 180.);
        tmpDouble = this->density * this->thickness / tmpDouble;
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

std::vector<std::pair<std::string, double> > Layer::getPeakFamilies(const double & energy, \
                                                                 const Elements & elements) const
{
    if (this->hasMaterial)
    {
        const std::map<std::string, double> & composition = this->material.getComposition();
        std::map<std::string, double> actualComposition;
        std::vector<std::string> elementsList;
        std::map<std::string, double>::const_iterator c_it;
        std::map<std::string, double>::const_iterator c_it2;

        for(c_it = composition.begin(); c_it != composition.end(); ++c_it)
        {
            actualComposition = elements.getComposition(c_it->first);
            for(c_it2 = actualComposition.begin(); c_it2 != actualComposition.end(); ++c_it2)
            {
                if (std::find(elementsList.begin(), elementsList.end(), c_it2->first) == elementsList.end())
                {
                    elementsList.push_back(c_it2->first);
                }
            }
        }
        return elements.getPeakFamilies(elementsList, energy);
    }
    else
    {
        return elements.getPeakFamilies(this->materialName, energy);
    }
}

std::ostream& operator<< (std::ostream& o, Layer const & layer)
{
    o << "Layer: " << layer.getMaterialName();
    o << " density(g/cm3) " << layer.getDensity();
    o << " thickness(cm) " << layer.getThickness();
    o << " funny " << layer.getFunnyFactor();
    return o;
}

} // namespace fisx
