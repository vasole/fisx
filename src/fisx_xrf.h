#ifndef FISX_XRF_H
#define FISX_XRF_H

class XRF
{
public:
    XRF();

    /*!
    Set the excitation beam
    */
    void setBeam();

    /*!
    Set the beam filters to be applied
    */
    void setBeamFilters();


    /*!
    Add a new beam filter to list of beam filters
    */
    void addBeamFilter();


    /*!
    Set the excitation geometry.
    For the time being, just the incident, outgoing angles and scattering angle to detector center.
    */
    void setGeometry();

    /*
    Set the sample description.
    It consists on a set of layers of different materials, densities and thicknesses.
    The top layer will be taken as reference layer. This can be changed calling setRefenceLayer
    */
    void setSample():

    /*
    Set the reference layer. The detector distance is measured from the reference layer surface.
    If not specified, the lauer closest to the detector
    */
    void setReferenceLayer(int);
    void setRefenceLayer(std::string name);

    /*!
    It consists on a set of layers of different materials, densities and thicknesses and a "funny" factor.
    */
    void setAttenuators();

    /*!
    Add a new attenuator to the list of attenuators
    */
    void addAttenuator();


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
    setDetector();

    /*!
    Methods coordinating all the calculation
    */
    void detectedEmission()
    void expectedEmission():
    void expectedFluorescence();
    void expectedScattering();
    void peakRatios();

private:
    /*
    The Attenuators have methods getTransmission(double Energy) and getTransmissions(std::vector<double> Energies)
    Beam beam;
    std::map<std::string, Attenuators>  beamFilters;
    geometry
    std::map<std::string, Attenuators> Sample; #Funny set to 1
    std::map<std::string, Attenuators>  attenuators;
    collimators Not implemented;
    detector returns a set of peaks (with labels!) for each incoming energy. Escape peaks handled. Sum peaks not yet.
    */
};
#endif
