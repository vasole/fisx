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
/*
    if (d <= 0.0)
    {
        // thick target
        return (muj/mu1) * std::log(1 + mu1/muj);
    }
*/
    // deal with the problematic term
    tmpDouble = (muj - mu2) * d;
    if (tmpDouble > 0.0)
    {
        tmpDouble = Math::deBoerD(tmpDouble) / (mu2 * (mu1 + mu2));
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

double Math::deBoerX(const double & p, const double & q, const double & d1, const double & d2)
{
    return Math::deBoerV(p, q, d1, d2) - Math::deBoerV(p, q, d1, 0.0) \
                                       - Math::deBoerV(p, q, 0.0, d2) + Math::deBoerV(p, q, 0.0, 0.0);
}

double Math::deBoerV(const double & p, const double & q, const double & d1, const double & d2)
{
    double mu1j, mu2j, mubj, db;
    double tmpDouble1;
    double tmpDouble2;

    // case V(0,0) with db equal to 0
    if ((db == 0) && (p == 0) && (q == 0))
    {
        tmpDouble2 = (mu2j / p) * std::log(1.0 + (p / mu2j));
        tmpDouble1 =  1.0 - (q / mu1j);
        if (tmpDouble1 > 0.0)
            tmpDouble1 = - tmpDouble2 - (mu1j / q) * std::log(tmpDouble1);
        else
            tmpDouble1 = - tmpDouble2;
        return tmpDouble1 / (p * mu1j + q * mu2j);
    }

    // X(p, q, d1, d2) = V(d1, d2) - V(d1, 0) - V(0, d2) + V(0, 0)
    // V(inf, d2) = 0.0;
    // V(d1, inf) = 0.0; provided that q - muij < 0
    // Therefore, for a layer on thick substrate
    // X (p, q, d1, inf) = - V(d1, 0) + V(0,0)
    // and for small values of d1 (thin layer on thick substrate) that gives
    // X (p, q, d1, inf)X (p, q, d1, inf) is about (d1/p) * std::log(1.0 + (p/mu2j))


    tmpDouble1 = Math::deBoerD((1.0 + (p / mu2j)) * (mu1j*d1 + mubj * db + mu2j*d2));
    tmpDouble1 *= (mu2j /(p * (p * mu1j + q * mu2j)));
    tmpDouble1 *= std::exp((q - mu1j) * d1 - (p + mu2j) * d2 -mubj * db);

    tmpDouble2 = mu1j * d1 + mubj * db + mu2j * d2;
    tmpDouble2 = ((mu1j / (q * (p * mu1j + q * mu2j))) * Math::deBoerD(( 1.0 - (q/mu1j)) * tmpDouble2)) - \
                 (Math::deBoerD(tmpDouble2)/(p * q));

    return tmpDouble1 + tmpDouble2;
}
