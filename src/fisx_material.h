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
#ifndef FISX_MATERIAL_H
#define FISX_MATERIAL_H
#include <string>
#include <vector>
#include <map>

namespace fisx
{

/*!
  \class Material
  \brief Class containing the composition of a material

   A material is nothing else than a name and a map of elements and mass fractions.

   The default density, default thickness and comment can be supplied for convenience purposes.
*/
class Material
{
public:
    /*!
    Minimalist constructor.
    */
    Material();

    /*!
    Expected constructor.
    */
    Material(const std::string & materialName, const double & density = 1.0, \
             const double & thickness = 1.0, const std::string & comment = "");

    void setName(const std::string & name);
    void initialize(const std::string & materialName, const double & density = 1.0,\
                    const double & thickness = 1.0, const std::string & comment="");

    /*!
    Set the composition of the material.
    This method normalizes the supplied amounts to make sure the sum is one.
    */
    void setComposition(const std::map<std::string, double> &);

    /*!
    Alternative method to set the composition of the material
    This method normalizes the supplied amounts to make sure the sum is one.
    */
    void setComposition(const std::vector<std::string> &, const std::vector<double> &);

    /*!
    Return the material composition as normalized mass fractions
    */
    std::map<std::string, double> getComposition() const;
    std::string getName() const;
    std::string getComment() const;
    double getDefaultDensity(){return this->defaultDensity;};
    double getDefaultThickness(){return this->defaultThickness;};

private:
    std::string name;
    bool initialized;
    std::map<std::string, double> composition;
    double defaultDensity;
    double defaultThickness;
    std::string    comment;
};

} // namespace fisx

#endif //FISX_MATERIAL_H
