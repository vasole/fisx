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
#ifndef FISX_EPDL97_H
#define FISX_EPDL97_H
#include <string>
#include <ctype.h>
#include <vector>
#include <map>

namespace fisx
{

class EPDL97
{
public:
    EPDL97();
    EPDL97(std::string directoryName);

    void setDataDirectory(std::string directoryName);
    // possibility to change binding energies
    void loadBindingEnergies(std::string fileName);
    void setBindingEnergies(const int & z, const std::map<std::string, double> & bindingEnergies);
    const std::map<std::string, double> & getBindingEnergies(const int & z);

    // the actual mass attenuation related functions
    std::map<std::string, double> getMassAttenuationCoefficients(const int & z, const double & energy) const;
    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients(const int & z,\
                                                const std::vector<double> & energy) const;

    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients(const int & z) const;

    // the vacancy distribution related functions
    std::map<std::string, double> getPhotoelectricWeights(const int & z, \
                                                          const double & energy);

    std::map<std::string, std::vector<double> > getPhotoelectricWeights(const int & z, \
                                                const std::vector<double> & energy);

    // utility functions
    std::string toUpperCaseString(const std::string &) const;
    std::pair<long, long> getInterpolationIndices(const std::vector<double> &,  const double &) const;

private:
    // internal function to load the data
    bool initialized;
    void loadData(std::string directoryName);
    void loadCrossSections(std::string fileName);

    // The directory name
    std::string directoryName;

    // The used file for binding energies
    std::string bindingEnergiesFile;

    // The used file for cross sections
    std::string crossSectionsFile;

    // The table containing all binding energies for all elements and shells
    // bindingEnergy[Z - 1]["K"] gives the binding energy for the K shell of element
    //                       with atomic number Z.
    std::vector<std::map<std::string, double> > bindingEnergy;

    // Mass attenuation data as read from the files
    // We have a table for each element but all of them share the same
    // file header structure
    std::vector<std::string> muInputLabels;
    std::map<std::string, int> muLabelToIndex;
    std::vector<std::vector<std::vector <double> > > muInputValues;
    std::vector<std::vector<double> > muEnergy;

    // Partial photoelectric mass attenuation coefficients
    // For each shell (= key), there is a vector for the energies
    // and a vector for the value of the mass attenuation coefficients
    // Expected map key values are:
    // K, L1, L2, L3, M1, M2, M3, M4, M5, "all other"
    void initPartialPhotoelectricCoefficients();
};

} // namespace fisx

#endif // FISX_EPDL97_H
