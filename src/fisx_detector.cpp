#include "fisx_detector.h"
#include <math.h>
#include <stdexcept>

Detector::Detector(const std::string & name, const double & density, const double & thickness, \
                  const double & funnyFactor): Layer(name, density, thickness, funnyFactor)
{
    ;
}

double Detector::getActiveArea()
{
    double pi;
    pi = acos(-1.0);
    return (0.5 * pi) * (this->diameter * this->diameter);
}

void Detector::setActiveArea(const double & area)
{
    double pi;
    pi = acos(-1.0);
    if (area < 0)
    {
        throw std::invalid_argument("Negative detector area");
    }
    this->diameter = sqrt(area/(0.5 * pi));
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

const double & Detector::getDiameter()
{
    return this->diameter;
}

const double & Detector::getDistance()
{
    return this->distance;
}
