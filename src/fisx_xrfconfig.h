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
#ifndef FISX_XRFCONFIG_H
#define FISX_XRFCONFIG_H
// TODO #include "fisx_version.h"
#include "fisx_detector.h"
#include "fisx_beam.h"

namespace fisx
{

class XRFConfig
{
public:
    XRFConfig();

    friend std::ostream& operator<< (std::ostream& o, XRFConfig const & config);

    void readConfigurationFromFile(const std::string & fileName);
    void saveConfigurationToFile(const std::string & fileName);

    /*!
    Set the excitation beam
    */
    void setBeam(const double & energy, const double & divergency);
    void setBeam(const std::vector<double> & energies, \
                 const std::vector<double> & weight, \
                 const std::vector<int> & characteristic = std::vector<int>(), \
                 const std::vector<double> & divergency = std::vector<double>());
    void setBeam(const Beam & beam);

    /*!
    Set the beam filters to be applied to the beam
    */
    void setBeamFilters(const std::vector<std::string> & names, \
                        const std::vector<double> & densities, \
                        const std::vector<double> & thicknesses, \
                        const std::vector<std::string> & comments);
    void setBeamFilters(const std::vector<Layer> & filters);

    /*!
    Set the excitation geometry.
    For the time being, just the incident, outgoing angles and scattering angle to detector center.
    */
    void setGeometry(const double & alphaIn, const double & alphaOut, const double & scatteringAngle = 90.);

    /*!
    Set the sample description.
    It consists on a set of layers representing different materials, densities and thicknesses.
    The first ( = top) layer will be taken as reference layer. This can be changed calling setRefenceLayer
    */
    void setSample(const std::vector<Layer> & layers, const int & referenceLayer = 0);
    void setSample(const std::vector<std::string> & names, \
                   const std::vector<double> & densities, \
                   const std::vector<double> & thicknesses, \
                   const std::vector<std::string> & comments,
                   const int & referenceLayer = 0);

    /*!
    Set the reference layer. The detector distance is measured from the reference layer surface.
    If not specified, the first layer is the reference layer (closest to the detector).
    */
    void setReferenceLayer(int referenceLayer);


    /*!
    Set the list of attenuators. Attenuators are layers between sample and detector.
    */
    void setAttenuators(const std::vector<Layer> & attenuators);
    void setAttenuators(const std::vector<std::string> & names, \
                        const std::vector<double> & densities, \
                        const std::vector<double> & thicknesses, \
                        const std::vector<std::string> & comments);

    /*!
    Collimators are not implemented yet. The collimators are attenuators that take into account their distance to
    the sample, their diameter, thickness and density
    */
    void setCollimators();
    void addCollimator();

    /*!
    Set the detector. For the time being it is very simple.
    It has active area, material, density, thickness and distance.
    */
    void setDetector(const Detector & detector);

    /*!
    Methods coordinating all the calculation
    */
    /*
    void detectedEmission()
    void expectedEmission():
    void expectedFluorescence();
    void expectedScattering();
    void peakRatios();
    */
    /*!
    Returns a constant reference to the internal beam.
    */
   const Beam & getBeam() const;
   const std::vector<Layer> & getBeamFilters() const {return this->beamFilters;};
   const std::vector<Layer> & getSample() const {return this->sample;};
   const std::vector<Layer> & getAttenuators() const {return this->attenuators;};
   const Detector & getDetector() const {return this->detector;};
   const double & getAlphaIn() const {return this->alphaIn;};
   const double & getAlphaOut() const {return this->alphaOut;};
   const double & getScatteringAngle() const {return this->scatteringAngle;};
   const int & getReferenceLayer() const {return this->referenceLayer;};

private:
    Beam beam;
    std::vector<Material> materials;
    std::vector<Layer> beamFilters;
    std::vector<Layer> sample;          // just other layer with funny factor set to 1.0
    std::vector<Layer> attenuators;
    int referenceLayer;
    double  alphaIn;
    double  alphaOut;
    double  scatteringAngle;
    // for the time being the detector is just other layer
    Detector detector;
    //collimators Not implemented;


    /*
    WARNING: If materials are not defined in terms of formulas or elemental compositions but in terms of other materials,
    the methods using a constant reference to an Elements instance can fail if all the materials are not present in the
    library or have been redefined:

    layer.getTransmission(energy, elementsInstance)
    layer.getTransmission(energies, elementsInstance)
    layer.getMassAttenuationCoefficients(energy, elementsInstance)
    layer.getMassAttenuationCoefficients(energies, elementsInstance)
    layer.getComposition(elementsInstance)
    detector.getEscape(energy, elementsInstance, const std::string & label = "", const int & update = 1)

    TODO: Implement those methods as XRF methods taken layers or detector as first argument.
    */
};

} // namespace fisx

#endif // FISX_SRF_CONFIG_H
