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
#ifndef FISX_SHELL_H
#define FISX_SHELL_H
#include <string>
#include <ctype.h>
#include <vector>
#include <map>


namespace fisx
{

class Shell
{
public:
    Shell(); // Do not use, just to be able to fit into a map container

    /*!
    Shell constructor. Only shells K, L1, L2, L3, M1, M2, M3, M4 or M5 accepted
    */
    Shell(std::string name);
    // ~Shell();

    /*!
    Set the radiative transitions probabilities.
    The transition labels have to start by the name used to instantiate the class.
    Valid examples: KL3, L3M5, ...
    The class normalizes the provided intensities.
    */
    void setRadiativeTransitions(std::map<std::string, double>);

    /*!
    Convenience method
    */
    void setRadiativeTransitions(std::vector<std::string>, std::vector<double> values);

    /*!
    Convenience method
    */
    void setRadiativeTransitions(const char *strings[], const double *values, int nValues);

    /*!
    Set the non-radiative transitions probabilities originating on this shell
    The transition labels have to start by the name used to instantiate the class, followed
    by a a hyphen, the origin shell of the filling electron and the shell from which the additional
    electron is ejected. (Examples: L1-L2M2, ...)
    The class normalizes the provided intensities separating Auger and Coster-Kronig ratios.
    */
    void setNonradiativeTransitions(std::vector<std::string>, std::vector<double> values);
    void setNonradiativeTransitions(std::map<std::string, double> values);
    void setNonradiativeTransitions(const char *strings[], const double *values, int nValues);

    /*!
    Get Auger transition ratios
    */
    const std::map<std::string, double> & getAugerRatios() const;

    /*!
    Get Coster-Kronig transition ratios
    */
    const std::map<std::string, std::map<std::string, double> > & getCosterKronigRatios() const;

    /*!
    Get X-ray fluorescence transition ratios
    */
    const std::map<std::string, double> & getFluorescenceRatios() const;

    const std::map<std::string, double> & getRadiativeTransitions() const;
    const std::map<std::string, double> & getNonradiativeTransitions() const;

    /*! Return the probabilities of direct transfer of a vacancy to a higher shell following
    an X-ray emission, an Auger transition and Coster-Kronig transitions (if any).
    Since the different types of transitions are normalized, it multiplies by the
    fluorescence, Auger or Coster-Kronig yields to get the probabilities.
    */
    std::map<std::string, double> getDirectVacancyTransferRatios(const std::string& destination) const;

    /*!
    Set the shell constants.
    omega denotes the fluorescence yield of the shell
    fij represents the Coster-Kronig yield
    The Auger yield is calculated as 1 - omega - sum(fij)
    WARNING: The original constants are not cleared
              Only those constants supplied will be overwritten!!!
    */
    void setShellConstants(std::map<std::string, double>);
    const std::map<std::string, double> & getShellConstants() const;
    bool StringToInteger(const std::string& str, int & number);

    double getFluorescenceYield() const;

private:
    std::string  name;
    int shellMainIndex;
    int subshellIndex;
    void _updateNonradiativeRatios();
    void _updateFluorescenceRatios();
    std::string toUpperCaseString(const std::string & str);
    // double    omega;
    // double* ck;
    // double  bindingEnergy;
    // map of the form {"omega": fluorescence_yield,
    //                  "fij": Coster-Kronig yield fij}
    std::map<std::string, double> shellConstants;
    std::map<std::string, double> radiativeTransitions;
    std::map<std::string, double> nonradiativeTransitions;
    std::map<std::string, double> augerRatios;
    std::map<std::string, std::map<std::string, double> > costerKronigRatios;
    std::map<std::string, double> fluorescenceRatios;
};

} // namespace fisx

#endif // FISX_SHELL_H
