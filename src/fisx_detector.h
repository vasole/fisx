#ifndef FISX_DETECTOR_H
#define FISX_DETECTOR_H
#include "fisx_layer.h"

/*!
  \class Detector
  \brief Class describing the detector.
*/

class Detector:public Layer
{
public:

    /*!
    layer like initialization.
    The detector is assumed to be cylindrical and the diameter is calculated.
    */
    Detector(const std::string & name="",  const double & density = 0.0,
                                           const double & thickness = 0.0,
                                           const double & funnyFactor = 1.0);


    /*!
    Active area in cm2.
    The detector is assumed to be cylindrical and the diameter is calculated.
    */
    void setActiveArea(const double & area);

    /*!
    Diameter in cm2.
    For the time being the detector is assumed to be cylindrical.
    */
    void setDiameter(const double & diameter);

    /*!
    Returns the active area in cm2
    */
    double getActiveArea() const;

    /*!
    Returns the diameter in cm2
    */
    const double & getDiameter() const;

    /*!
    Sets the distance to reference layer in cm2
    */
    void setDistance(const double & distance);

    /*!
    Returns the distance to reference layer in cm2
    */
    const double & getDistance() const;

    /*!
    Returns escape peak energy and rate associated to given energy.

    The optional arguments label and update serve for caching purposes.
    */
    std::map<std::string, std::map<std::string, double> > getEscape(const double & energy, \
                                                            const Elements & elementsLibrary, \
                                                            const std::string & label = "", \
                                                            const int & update = 1);

    void setMinimumEscapePeakEnergy(const double & energy);
    void setMinimumEscapePeakIntensity(const double & intensity);
    void setMaximumNumberOfEscapePeaks(const int & nPeaks);

private:
    double diameter ;
    double distance ;
    // Escape peak related parameters
    double escapePeakEnergyThreshold;
    double escapePeakIntensityThreshold;
    int escapePeakNThreshold;
    double escapePeakAlphaIn;
    std::map< std::string, std::map<std::string, std::map<std::string, double> > > escapePeakCache;
    // TODO: Calibration, fano, noise, and so on.
};

#endif //FISX_DETECTOR_H
