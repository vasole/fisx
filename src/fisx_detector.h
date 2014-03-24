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
    double getActiveArea();

    /*!
    Returns the diameter in cm2
    */
    const double & getDiameter();

    /*!
    Sets the distance to reference layer in cm2
    */
    void setDistance(const double & distance);

    /*!
    Returns the distance to reference layer in cm2
    */
    const double & getDistance();


private:
    double diameter ;
    double distance ;
    // TODO: Calibration, fano, noise, and so on.
};

#endif //FISX_DETECTOR_H
