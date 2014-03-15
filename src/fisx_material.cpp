#include <stdexcept>
#include "fisx_material.h"

Material::Material()
{
    this->initialized = false;
    this->name = "Unset name";
    this->comment = "";
    this->defaultDensity = 1.0;
    this->defaultThickness = 1.0;
}

Material::Material(const std::string & materialName)
{
    this->initialize(materialName, 1.0, 1.0, "");
}

Material::Material(const std::string & materialName, const double & density, const double & thickness)
{
    this->initialize(materialName, density, thickness, "");
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

    if (this->initialized)
    {
        msg = "Material::initialize. Material is already initialized with name " + this->name;
        throw std::invalid_argument(msg);
    }

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

std::map<std::string, double> Material::getComposition()
{
    return this->composition;
}

std::string Material::getName()
{
    return this->name;
}

