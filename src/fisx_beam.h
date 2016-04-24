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
#ifndef FISX_BEAM_H
#define FISX_BEAM_H
#include <cstddef> // needed for NULL definition!!!
#include <vector>
#include <iostream>

namespace fisx
{

struct Ray
{
    double energy;
    double weight;
    int characteristic;
    double divergency;

    Ray()
    {
        energy = 0.0;
        weight = 0.0;
        characteristic = 0;
        divergency = 0.0;
    }

    bool operator < (const Ray &b) const
    {
        return (energy < b.energy);
    }
};

/*!
  \class Beam
  \brief Class describing an X-ray beam

   At this point a beam is described by a set of energies and weights. The characteristic flag just indicates if
   it is an energy to be considered for calculation of scattering peaks.
*/
class Beam
{
public:
    Beam();
    /*!
    Minimalist constructor.
    */

    /*!
    Beam description given as vectors.
    The beam is always internally ordered.
    */
    void setBeam(const std::vector<double> & energy, \
                 const std::vector<double> & weight = std::vector<double>(),\
                 const std::vector<int> & characteristic = std::vector<int>(),\
                 const std::vector<double> & divergency = std::vector<double>());

    friend std::ostream& operator<< (std::ostream& o, Beam const & beam);

    /*!
    Easy to wrap interface functions
    Except for the energy, you can use NULL pointers to use default values.
    */
    void setBeam(const int & nValues, const double *energy, const double *weight = NULL,
                 const int *characteristic = NULL, const double *divergency = NULL);

    void setBeam(const int & nValues, const double *energy, const double *weight = NULL,
                 const double *characteristic = NULL, const double *divergency = NULL);

    void setBeam(const double & energy, const double divergency = 0.0);

    /*!
    Returns a constant reference to the internal beam.
    */
    const std::vector<Ray> & getBeam();

    /*!
    Currently it returns a vector of "elements" in which each element is a vector of
    doubles with length equal to the number of energies.
    The first four elements are warranteed to exist and they are:
    [energy0, energy1, energy2, ...]
    [weight0, weight1, weight2, ...]
    [characteristic0, characteristic1, charactersitic2, ...]
    [divergency0, divergency1, divergency2, ...]
    */
    std::vector<std::vector<double> > getBeamAsDoubleVectors() const;

private:
    bool normalized;
    void normalizeBeam(void);
    std::vector<Ray> rays;
};

} // namespace fisx

#endif //FISX_BEAM_H
