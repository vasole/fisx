#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2023 European Synchrotron Radiation Facility
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
#ifndef FISX_XRF_H
#define FISX_XRF_H
#include "fisx_xrfconfig.h"
#include "fisx_elements.h"
#include <iostream>

namespace fisx
{

class XRF
{

typedef std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
        expectedLayerEmissionType;
typedef std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > >::const_iterator \
        iteratorExpectedLayerEmissionType;

public:
    /*!
    Default constructor
    */
    XRF();

    /*!
    Constructor with configuration file
    */
    XRF(const std::string & configurationFile);

    /*!
    Read the configuration from file
    */
    void readConfigurationFromFile(const std::string & fileName);

    /*!
    Set the excitation beam
    */
    void setBeam(const Beam & beam);

    /*!
    Easy to wrap funtion to set the excitation beam
    \param energies Set of double values corresponding to energies (in keV) describing
           the incoming beam
    \param weight Set of weights with the relative intensity of each energy
    \param characteristic Integer flag currently ignored by fisx library
    \param divergency Beam divergency in degrees (currently ignored by fisx library)
    */
    void setBeam(const std::vector<double> & energies, \
                 const std::vector<double> & weight, \
                 const std::vector<int> & characteristic = std::vector<int>(), \
                 const std::vector<double> & divergency = std::vector<double>());

    /*!
    Funtion to set a single energy excitation beam
    \param energy Energy of the incoming beam
    */
    void setSingleEnergyBeam(const double & energy, const double & divergency = 0.0);


    /*!
    Set the beam filters to be applied
    */
    void setBeamFilters(const std::vector<Layer> & filters);

    /*!
    Set the user beam filters to be applied
    */
    void setUserBeamFilters(const std::vector<TransmissionTable> & userFilters);


    /*!
    Set the sample description.
    It consists on a set of layers of different materials, densities and thicknesses.
    The top layer will be taken as reference layer. This can be changed calling setRefenceLayer
    */
    void setSample(const std::vector<Layer> & layers, const int & referenceLayer = 0);


    /*!
    Convenience method for single layer samples.
    */
    void setSample(const std::string & name, \
                   const double & density = 1.0, \
                   const double & thickness = 1.0);

    /*!
    Convenience method for single layer samples.
    */
    void setSample(const Layer & layer);



    /*!
    It consists on a set of layers of different materials, densities and thicknesses and
    "funny" factors.
    */
    void setAttenuators(const std::vector<Layer> & attenuators);


    /*!
    Set the user transmission tables to be applied as attenuators
    */
    void setUserAttenuators(const std::vector<TransmissionTable> & userAttenuators);


    /*!
    Set the detector. For the time being it is very simple.
    It has active area/diameter, material, density, thickness and distance.
    */
    void setDetector(const Detector & detector);

    /*!
    Set the excitation geometry.
    For the time being, just the incident, outgoing angles and scattering angle to detector
    center. A negative scattering angle of 90 degrees indicates the scattering angle is the
    sum of alphaIn and alphaOut.
    */
    void setGeometry(const double & alphaIn, const double & alphaOut,\
                      const double & scatteringAngle = -90.);

    /*!
    Set the reference layer. The detector distance is measured from the reference layer surface.
    If not specified, the layer closest to the detector
    */
    void setReferenceLayer(const int & index);
    void setRefenceLayer(const std::string & name);

    /*!
    Set the elements library to be used.
    */
    //void setElementsReference(const Elements & elements);

    /*!
    Collimators are not implemented yet. The collimators are attenuators that take into account their distance to
    the sample, their diameter, thickness and density
    */
    void setCollimators();
    void addCollimator();

    /*!
    For debugging purposes
    */
    void printConfiguration() const;

    /*!
    Get the current configuration
    */
    const XRFConfig & getConfiguration() const;

    /*!
    Set the configuration
    */
    void setConfiguration(const XRFConfig & configuration);

    /*!
    Get the expected fluorescence emission coming from primary excitation per unit photon.
    It needs to be multiplied by the mass fraction and the total number of photons to get
    the actual primary fluorescence.

    The output is a map:

    Element -> Family -> Line -> energy: double, ratio: double

    */
    std::map< std::string, std::map< std::string, std::map<std::string, std::map<std::string, double> > > >\
                getExpectedPrimaryEmission(const std::vector<std::string> & elementList,
                                           const Elements & elements);

    /*!
    Get geometric efficiency term. Solid angle subtended by the detector at the sample divided by 4 * pi
    */
    double getGeometricEfficiency(const int & layerIndex = 0) const;


    /*!
    Return a complete output of the form
    [Element Family][Layer][line]["energy"] - Energy in keV of the emission line
    [Element Family][Layer][line]["primary"] - Primary rate prior to correct for detection efficiency
    [Element Family][Layer][line]["secondary"] - Secondary rate prior to correct for detection efficiency
    [Element Family][Layer][line]["rate"] - Overall rate
    [Element Family][Layer][line]["efficiency"] - Detection efficiency
    [Element Family][Layer][line][element line layer] - Secondary rate (prior to correct for detection efficiency)
    due to the fluorescence from the given element, line and layer index composing the map key.
    */
    std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                getMultilayerFluorescence(const std::string & element, \
                const Elements & elementsLibrary, const int & sampleLayerIndex = 0, \
                const std::string & lineFamily = "", const int & secondary = 0, \
                const int & useGeometricEfficiency = 1, const int & useMassFractions = 0, \
                const double & secondaryCalculationLimit = 0.0, \
                const Beam & overwritingBeam = Beam());
    /*!
    \param elementFamilyLayer - Vector of strings. Each string represents the information we are interested on.\n
    "Cr"     - We want the information for Cr, for all line families and sample layers\n
    "Cr K"   - We want the information for Cr, for the family of K-shell emission lines, in all layers.\n
    "Cr K 0" - We want the information for Cr, for the family of K-shell emission lines, in layer 0.
    \param elementsLibrary - Instance of library to be used for all the Physical constants
    \param secondary - Flag to indicate different levels of secondary excitation to be considered.\n
                0 Means not considered\n
                1 Consider secondary excitation\n
                2 Consider tertiary excitation\n
    \param useGeometricEfficiency - Take into account solid angle or not. Default is 1 (yes)

    \param useMassFractions - If 0 (default) the output corresponds to the requested information if the mass
    fraction of the element would be one on each calculated sample layer. To get the actual signal, one
    has to multiply the rates by the actual mass fraction of the element on each sample layer.
                       If set to 1, the rate will be already corrected by the actual mass fraction.

    \return Return a complete output of the form:\n
    [Element Family][Layer][line]["energy"] - Energy in keV of the emission line\n
    [Element Family][Layer][line]["primary"] - Primary rate prior to correct for detection efficiency\n
    [Element Family][Layer][line]["secondary"] - Secondary rate prior to correct for detection efficiency\n
    [Element Family][Layer][line]["rate"] - Overall rate\n
    [Element Family][Layer][line]["efficiency"] - Detection efficiency\n
    [Element Family][Layer][line][element line layer] - Secondary rate (prior to correct for detection efficiency)
    due to the fluorescence from the given element, line and layer index composing the map key.\n
    [Element Family][Layer][line]["massFraction"] - Mass fraction of the element in the considered layer
    */
    std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                getMultilayerFluorescence(const std::vector<std::string> & elementFamilyLayer, \
                const Elements & elementsLibrary, \
                const int & secondary = 0, \
                const int & useGeometricEfficiency = 1, \
                const int & useMassFractions = 0, \
                const double & secondaryCalculationLimit = 0.0,
                const Beam & overwritingBeam = Beam());

    /*!
    Basis method called by all the other convenience methods.
    \param elementList - Vector of strings. Each string represents one element.\n
    \param layerList - Vector of integers representing the sample layers where the calculation for the elements in
    the elements list takes place. Either a length of 1 or of the same length as the element list.
    \param familyList - Vector of strings. Each string represents one peak family.\n
    The vector has to be of the same size as the elements list.
    \param elementsLibrary - Instance of library to be used for all the Physical constants
    \param secondary - Flag to indicate different levels of secondary excitation to be considered.\n
                0 Means not considered\n
                1 Consider secondary excitation\n
                2 Consider tertiary excitation\n
    \param useGeometricEfficiency - Take into account solid angle or not. Default is 1 (yes)

    \param useMassFractions - If 0 (default) the output corresponds to the requested information if the mass
    fraction of the element would be one on each calculated sample layer. To get the actual signal, one
    has to multiply the rates by the actual mass fraction of the element on each sample layer.
                       If set to 1, the rate will be already corrected by the actual mass fraction.

    \param overwritingBeam - Optional parameter to use a different beam than the one in the configuration.

    \return Return a complete output of the form:\n
    [Element Family][Layer][line]["energy"] - Energy in keV of the emission line\n
    [Element Family][Layer][line]["primary"] - Primary rate prior to correct for detection efficiency\n
    [Element Family][Layer][line]["secondary"] - Secondary rate prior to correct for detection efficiency\n
    [Element Family][Layer][line]["rate"] - Overall rate\n
    [Element Family][Layer][line]["efficiency"] - Detection efficiency\n
    [Element Family][Layer][line][element line layer] - Secondary rate (prior to correct for detection efficiency)
    due to the fluorescence from the given element, line and layer index composing the map key.\n
    [Element Family][Layer][line]["massFraction"] - Mass fraction of the element in the considered layer

    */
    std::map<std::string, std::map<int, std::map<std::string, std::map<std::string, double> > > > \
                getMultilayerFluorescence(const std::vector<std::string> & elementList,
                                          const Elements & elementsLibrary, \
                                          const std::vector<int> & layerList, \
                                          const std::vector<std::string> &  familyList, \
                                          const int & secondary = 0, \
                                          const int & useGeometricEfficiency = 1, \
                                          const int & useMassFractions = 0, \
                                          const double & secondaryCalculationLimit = 0.0,
                                          const Beam & overwritingBeam = Beam()) const;


    double getEnergyThreshold(const std::string & elementName, const std::string & family, \
                                const Elements & elementsLibrary) const;

    /*!
    Get the beam
    */
    const Beam & getBeam() const {return this->configuration.getBeam();};

    /*!
    Get the beam filters to be applied
    */
    const std::vector<Layer> & getBeamFilters() const {return this->configuration.getBeamFilters();};

    /*!
    Get the user beam filters to be applied
    */
    const std::vector<TransmissionTable> & getUserBeamFilters() const {return this->configuration.getUserBeamFilters();};

    /*!
    Get the sample description
    */
    const std::vector<Layer> & getSample() const {return this->configuration.getSample();};

    /*!
    Get the attenuator layers
    */
    const std::vector<Layer> & getAttenuators() const {return this->configuration.getAttenuators();};

    /*!
    Get the user attenuators to be applied
    */
    const std::vector<TransmissionTable> & getUserAttenuators() const {return this->configuration.getUserAttenuators();};

    /*!
    Get the detector
    */
    const Detector & getDetector() const {return this->configuration.getDetector();};


    const double & getAlphaIn() const {return this->configuration.getAlphaIn();};
    const double & getAlphaOut() const {return this->configuration.getAlphaOut();};
    const double & getScatteringAngle() const {return this->configuration.getScatteringAngle();};
    const int & getReferenceLayer() const {return this->configuration.getReferenceLayer();};
    const std::vector<Material> & getMaterials() const {return this->configuration.getMaterials();};

    /*
    Replacement of layer methods accounting for materials defined in the configuration.
    Only the XRF instance knows about all the materials defined. That is not the case for the elements library.
    */

    /*!
    Get the layer composition either from internal material or from name.
    In the later case the elements instance is not used.
    */
    std::map<std::string, double> getLayerComposition(const Layer & layer, const Elements & elements) const;

    /*!
    Get the layer mass attenuation coefficients transmission at the given energy using the elements library
    supplied but accounting for the materials defined in the configuration.
    If the layer material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    std::map<std::string, double> getLayerMassAttenuationCoefficients(const Layer & layer,
                                                                      const double & energy,
                                                                      const Elements & elements,
                                                const std::map<std::string, double> & layerComposition = std::map<std::string, double>()) const;


    /*!
    Get the layer mass attenuation coefficients transmission at the given energies using the elements
    library supplied but accounting for the materials defined in the configuration.
    If a material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    std::map<std::string, std::vector<double> > getLayerMassAttenuationCoefficients( \
                                                const Layer & layer,
                                                const std::vector<double> & energies,
                                                const Elements & elements,
                                                const std::map<std::string, double> & layerComposition = std::map<std::string, double>()) const;
    /*!
    Get the layer transmissions at the given energy using the elements library
    supplied but accounting for the materials defined in the configuration.
    If a material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    double getLayerTransmission(const Layer & layer, const double & energy, const Elements & elements, const double & angle = 90.0, \
                                                const std::map<std::string, double> & layerComposition = std::map<std::string, double>()) const;

    /*!
    Get the layer transmissions at the given energies using the elements library
    supplied but accounting for the materials defined in the configuration.
    If a material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    std::vector<double> getLayerTransmission(const Layer & layer,
                                             const std::vector<double> & energy,
                                             const Elements & elements,
                                             const double & angle = 90.0,
                                             const std::map<std::string, double> & layerComposition = std::map<std::string, double>()) const;

    /*!
    Get the emitted peak families from a layer given an energy and a reference to an elements library.
    It uses the supplied Elements library accounting for the materials defined in the configuration.
    It returns an ordered vector of pairs.
    The first element is the peak family ("Si K", "Pb L1", ...) and the second the binding energy.
    */
    std::vector<std::pair<std::string, double> > getLayerPeakFamilies(const Layer & layer,
                                             const double & energy,
                                             const Elements & elements,
                                             const std::map<std::string, double> & layerComposition = std::map<std::string, double>()) const;
private:
    /*!
    Reference to elements library to be used for calculations
    */
    /*!
    The internal configuration
    */
    XRFConfig configuration;

    /*!
    Some optimization flags
    */
    bool recentBeam;

    expectedLayerEmissionType lastMultilayerFluorescence;
};

} // namespace fisx

#endif // FISX_XRF_H
