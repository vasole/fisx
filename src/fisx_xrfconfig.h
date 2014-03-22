#ifndef FISX_XRFCONFIG_H
#define FISX_XRFCONFIG_H
// TODO #include "fisx_version.h"
#include "fisx_layer.h"
#include "fisx_beam.h"

class XRFConfig
{
public:
    XRFConfig();

    void readConfigurationFromFile(const std::string & fileName);
    void saveConfigurationToFile(const std::string & fileName);

    /*!
    Set the excitation beam
    */
    void setBeam(const double & energy, const double & weight = 1.0, \
                 const int & characteristic = 1, const double & divergency = 0.0);
    void setBeam(const std::vector<double> & energies, \
                 const std::vector<double> & weight, \
                 const std::vector<int> & characteristic = std::vector<int>(), \
                 const std::vector<double> & divergency = std::vector<double>());
    void setBeam(const Beam & beam);

    /*!
    Set the beam filters to be applied to the beam
    */
    void setBeamFilters(const std::vector<std::string> & names,\
                        const std::vector<double> & densities,\
                        const std::vector<double> & thicknesses,\
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
    void setSample(const std::vector<std::string> & names,\
                   const std::vector<double> & densities,\
                   const std::vector<double> & thicknesses,\
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
    void setAttenuators(const std::vector<std::string> & names,\
                        const std::vector<double> & densities,\
                        const std::vector<double> & thicknesses,\
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
    //void setDetector();

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

private:
    Beam beam;
    std::vector<Material> materials;
    std::vector<Layer> beamFilters;
    std::vector<Layer> sample;
    std::vector<Layer> attenuators;
    int referenceLayer;
    double  alphaIn;
    double  alphaOut;
    double  scatteringAngle;
    Layer detector;
    //collimators Not implemented;


    /*
    The Attenuators have methods getTransmission(double Energy) and getTransmissions(std::vector<double> Energies)
    std::map<std::string, Attenuators>  beamFilters;
    geometry
    std::map<std::string, Attenuators> Sample; #Funny set to 1
    std::map<std::string, Attenuators>  attenuators;
    detector returns a set of peaks (with labels!) for each incoming energy. Escape peaks handled. Sum peaks not yet.
    */
};
#endif
