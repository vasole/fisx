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
#include <stdexcept>
#include "fisx_material.h"
#include <iostream>

namespace fisx
{

Material::Material()
{
    this->initialized = false;
    this->name = "Unset name";
    this->comment = "";
    this->defaultDensity = 1.0;
    this->defaultThickness = 1.0;
}

Material::Material(const std::string & materialName, const double & density,\
                    const double & thickness, const std::string & comment)
{
    this->initialize(materialName, density, thickness, comment);
}


void Material::setName(const std::string & name)
{
    std::string msg;

    if (this->initialized)
    {
        msg = "Material::setName. Material is already initialized with name " + this->name;
        throw std::invalid_argument(msg);
    }
    this->initialize(name, this->defaultDensity, this->defaultThickness, this->comment);
}

void Material::initialize(const std::string & materialName, const double & density, \
                          const double & thickness, const std::string & comment)
{
    std::string msg;

/*
    if (this->initialized)
    {
        msg = "Material::initialize. Material is already initialized with name " + this->name;
        throw std::invalid_argument(msg);
    }
*/
    if(materialName.size() < 1)
    {
        throw std::invalid_argument("Material name should have at least one letter");
    }
    if (density <= 0.0)
    {
        throw std::invalid_argument("Material density should be positive");
    }
    if (thickness <= 0.0)
    {
        throw std::invalid_argument("Material thickness should be positive");
    }

    this->name = materialName;
    this->defaultDensity = density;
    this->defaultThickness = thickness;
    this->comment = comment;
    this->initialized = true;
}

void Material::setComposition(const std::map<std::string, double> & composition)
{
    std::map<std::string, double>::const_iterator c_it;
    std::vector<std::string> names;
    std::vector<double> amounts;


    for(c_it = composition.begin(); c_it != composition.end(); ++c_it)
    {
        names.push_back(c_it->first);
        amounts.push_back(c_it->second);
    }
    this->setComposition(names, amounts);
}

void Material::setComposition(const std::vector<std::string> & names, const std::vector<double> & amounts)
{
    std::vector<double>::size_type i;
    double total;

    if (names.size() != amounts.size())
    {
        for(i = 0; i < names.size(); i++)
        {
            std::cout << i << " name " << names[i] << std::endl;
        }
        for(i = 0; i < amounts.size(); i++)
        {
            std::cout << i << " amount " << amounts[i] << std::endl;
        }
        throw std::invalid_argument("Number of substances does not match number of amounts");
    }

    total = 0.0;

    for (i = 0; i < amounts.size(); ++i)
    {
        if (amounts[i] <= 0.0)
        {
            throw std::invalid_argument("Mass fractions cannot be negative");
        }
        total += amounts[i];
    }

    this->composition.clear();

    for (i = 0; i < amounts.size(); ++i)
    {
        this->composition[names[i]] = amounts[i] / total;
    }

}

std::map<std::string, double> Material::getComposition() const
{
    return this->composition;
}

std::string Material::getName() const
{
    return this->name;
}

std::string Material::getComment() const
{
    return this->comment;
}

} // namespace fisx
