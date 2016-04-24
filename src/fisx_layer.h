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
#ifndef FISX_LAYER_H
#define FISX_LAYER_H
#include <string>
#include <vector>
#include <iostream>
#include "fisx_elements.h"

namespace fisx
{

/*!
  \class Layer
  \brief Class containing the composition of a layer
*/

class Layer
{
public:

    Layer(const std::string & name="", const double & density = 0.0,
                                       const double & thickness = 0.0,
                                       const double & funnyFactor = 1.0);

    friend std::ostream& operator<< (std::ostream& o, Layer const & layer);

    void setMaterial(const std::string & materialName);

    void setMaterial(const Material & material);

    std::string name;

    /*!
    Get the layer composition either from internal material or from name.
    In the later case the elements instance is not used.
    */
    std::map<std::string, double> getComposition(const Elements & elements) const;

    /*!
    Eventually get a handle to underlying material
    */
    const Material & getMaterial() const {return this->material;};

    /*!
    Get the layer transmission at the given energy using the elements library
    supplied.
    If the material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    double getTransmission(const double & energy, const Elements & elements, const double & angle = 90.0) const;


    /*!
    Get the layer mass attenuation coefficients transmission at the given energy using the elements library
    supplied. If the material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    std::map<std::string, double> getMassAttenuationCoefficients(const double & energy,
                                                                 const Elements & elements) const;


    /*!
    Get the layer mass attenuation coefficients transmission at the given energies using the elements
    library supplied. If the material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients( \
                                                                const std::vector<double> & energies,
                                                                const Elements & elements) const;


    /*!
    Get the layer transmissions at the given energies using the elements library
    supplied.
    If the material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    std::vector<double> getTransmission(const std::vector<double> & energy,
                                const Elements & elements, const double & angle = 90.0) const;

    /*!
    Return true if material composition was specified.
    Returns false if only the name or formula of the material was given.
    */
    bool hasMaterialComposition() const;

    /*!
    Given an energy and a reference to an elements library return an ordered vector of pairs.
    The first element is the peak family ("Si K", "Pb L1", ...) and the second the binding energy.
    */
    std::vector<std::pair<std::string, double> > getPeakFamilies(const double & energy, \
                                                                 const Elements & elements) const;


    const std::string & getMaterialName()const {return this->materialName;};

    void setDensity(const double & density);
    void setThickness(const double & thickness);
    const double & getDensity()const {return this->density;};
    const double & getThickness()const {return this->thickness;};
    const double & getFunnyFactor()const {return this->funnyFactor;};

    std::map<std::string, double> getComposition(const Elements & elements);

private:
    std::string materialName;
    bool hasMaterial;
    Material material;
    double funnyFactor;
    double density;
    double thickness;
};

} // namespace fisx

#endif //FISX_LAYER_H
