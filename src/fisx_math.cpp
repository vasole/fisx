#include "fisx_math.h"
#include <cmath>
#include <stdexcept>
#include <iostream>

double Math::E1(const double & x)
{
    if (x <= 0)
    {
		throw std::invalid_argument("E1(x) Invalid argument. x should be positive");
    }
    if(x < 1)
    {
        return Math::AS_5_1_53(x);
    }
    else
    {
        return Math::AS_5_1_56(x);
    }
}

double Math::AS_5_1_53(const double & x)
{
    // Euler's gamma = 0.577215664901532860
    double a[6] = {-0.57721566, 0.99999193, -0.24991055,\
                    0.05519968, -0.00976004, 0.00107857};
    double result;
    if (x > 1)
    {
        throw std::invalid_argument("AS_5_1_53(x) Invalid argument. 0 < x <= 1");
    }

    result = a[5] * x;

    for (int i = 4; i > 0; --i)
    {
        result = (result + a[i]) * x;
    }
    result = result + a[0] - std::log(x);
    return result;
}

double Math::AS_5_1_56(const double & x)
{
    double a[4] = {8.5733287401, 18.0590169730, 8.6347608925, 0.2677737343};
    double b[4] = {9.5733223454, 25.6329561486, 21.0996530827, 3.9584969228};
    double num, den, result;

    if (x < 1)
    {
		throw std::invalid_argument("AS_5_1_56(x) Invalid argument. 1 <= x ");
    }
    num = x;
    den = x;
    for (int i = 4; i > 0; --i)
    {
        num = (num + a[i]) * x;
        den = (den + b[i]) * x;
    }
    result = (num / den) * (std::exp(-x)/x);
    return result;
}

double Math::deBoerD(const double & x)
{
    if  (x > 1)
    {
        // follow AS_5_1_56
        double a[4] = {8.5733287401, 18.0590169730, 8.6347608925, 0.2677737343};
        double b[4] = {9.5733223454, 25.6329561486, 21.0996530827, 3.9584969228};
        double num, den;
        num = x;
        den = x;
        for (int i = 4; i > 0; --i)
        {
            num = (num + a[i]) * x;
            den = (den + b[i]) * x;
        }
        return (num/den)/x;
    }
    else
    {
        return std::exp(x) * Math::AS_5_1_53(x);
    }
}


double Math::deBoerL0(const double & mu1, const double & mu2, const double & muj, \
                               const double & density, const double & thickness)
{
    double d;
    double tmpDouble;

    // express the thickness in g/cm2
    d = thickness * density;

    // deal with the problematic term
    tmpDouble = (muj - mu2) * d;
    if (tmpDouble > 0.0)
    {
        tmpDouble = Math::deBoerD(tmpDouble / (mu2 * (mu1 + mu2)));
    }
    else
    {
        tmpDouble = 0.0;
    }
    tmpDouble = -(Math::deBoerD(muj * d) / (mu1 * mu2)) + \
                 (Math::deBoerD((muj + mu1) * d) / (mu1 * (mu1 + mu2)));
    tmpDouble *= std::exp(-(mu1 + muj) * d);

    tmpDouble += std::log(1.0 + (mu1/muj)) / (mu1 * (mu1 + mu2));

    if (mu2 < muj)
    {
        tmpDouble += (std::exp(-(mu1 + mu2) * d) / (mu2 * (mu1 + mu2))) * \
                      std::log(1.0 - (mu2 / muj));
    }
    return tmpDouble;
}
