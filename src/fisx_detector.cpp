#include "fisx_detector.h"
#include <math.h>
#include <stdexcept>

Detector::Detector(const std::string & name, const double & density, const double & thickness, \
                  const double & funnyFactor): Layer(name, density, thickness, funnyFactor)
{
    this->diameter =  0.0;
    this->distance = 10.0;
}

double Detector::getActiveArea() const
{
    double pi;
    pi = acos(-1.0);
    return (0.25 * pi) * (this->diameter * this->diameter);
}

void Detector::setActiveArea(const double & area)
{
    double pi;
    pi = acos(-1.0);
    if (area < 0)
    {
        throw std::invalid_argument("Negative detector area");
    }
    this->diameter = 2.0 * sqrt(area/pi);
}

void Detector::setDiameter(const double & diameter)
{
    if (diameter < 0)
    {
        throw std::invalid_argument("Negative detector diameter");
    }
    this->diameter = diameter;
}

void Detector::setDistance(const double & distance)
{
    if (distance <= 0)
    {
        throw std::invalid_argument("Negative detector distance");
    }
    this->distance = distance;
}

const double & Detector::getDiameter() const
{
    return this->diameter;
}

const double & Detector::getDistance() const
{
    return this->distance;
}
