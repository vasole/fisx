#ifndef FISX_XRF_H
#define FISX_XRF_H
#include "fisx_xrfconfig.h"
#include "fisx_elements.h"
#include <iostream>

class XRF
{
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
    */
    void setBeam(const std::vector<double> & energies, \
                 const std::vector<double> & weight, \
                 const std::vector<int> & characteristic = std::vector<int>(), \
                 const std::vector<double> & divergency = std::vector<double>());

    /*!
    Funtion to set a single energy excitation beam
    */
    void setBeam(const double & energy, const double & divergency = 0.0);


    /*!
    Set the beam filters to be applied
    */
    void setBeamFilters(const std::vector<Layer> & filters);

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
    If not specified, the lauer closest to the detector
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
    Get the current configuration
    */
    const XRFConfig & getConfiguration();

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
    Methods coordinating all the calculation
    void detectedEmission()
    void expectedEmission():
    void expectedFluorescence();
    void expectedScattering();
    void peakRatios();
    */
    double getGeometricEfficiency(const int & layerIndex = 0) const;

    std::map<std::string, std::map<std::string, double> > getFluorescence(const std::string element, \
                const Elements & elementsLibrary, const int & sampleLayerIndex = 0, \
                const std::string & lineFamily = "", const int & secondary = 0, \
                const int & useGeometricEfficiency = 1);

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
};
#endif
