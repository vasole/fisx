#include "fisx_math.h"
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <iostream>

double Math::E1(const double & x)
{
    if (x == 0)
    {
        throw std::invalid_argument("E1(x) Invalid argument. x cannot be 0");
    }
    //if (x < -1)
    //{
    // throw std::invalid_argument("Not implemented. x cannot be less than -1");
    //}
    if (x < 0)
    {
        // AS 5.1.11
        // Euler's gamma = 0.577215664901532860;
        // I decide to evaluate just 10 terms of the series
        double result;
        double factorial[11] = {1.0, 1.0, 2.0, 6.0, 24., 120.0, 720., 5040., 40320., 362880., 3628800.};
        result = -0.577215664901532860;
        for(int n = 10; n > 0; --n)
        {
            result -= pow(-x, n) /(n * factorial[n]);
        }
        return result - std::log(-x);
    }
    if(x < 1)
    {
        return Math::AS_5_1_53(x) - log(x);
    }
    else
    {
        if (0)
        {
            // rational approximation is less accurate
            return Math::AS_5_1_56(x) / (x * std::exp(x));
        }
        else
        {
            return std::exp(-x) * Math::_deBoerD(x);
        }
    }
}

double Math::AS_5_1_53(const double & x)
{
    // Returns E1(x) + log(x)
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
    result = result + a[0];
    return result;
}

double Math::AS_5_1_56(const double & x)
{
    // Returns x * exp(x) * E1(x)
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
    result = (num / den);
    return result;
}

double Math::deBoerD(const double & x)
{

#ifndef NDEBUG
    // AS 5.1.19
    // 1 / (x + 1) < exp(x) * E1(x) < (1.0 /x)
    // AS 5.1.20
    // 0.5 * log(1 + 2.0/x) < exp(x) * E1(x) < log(1 + 1.0 /x)
    double tmpResult;
    double limit0, limit1;
    if (x < 0)
    {
        return std::exp(x) * E1(x);
    }
    if (x > 1)
    {
        tmpResult = Math::_deBoerD(x);
        //tmpResult = Math::AS_5_1_56(x) / x;
    }
    else
        tmpResult = std::exp(x) * (Math::AS_5_1_53(x) - log(x));
    limit0 = 0.5 * log(1 + 2.0/x);
    limit1 = log(1 + 1.0 /x);
    if ((tmpResult < limit0) || (tmpResult > limit1))
    {
        std::cout << "deBoerD error with x = " << x << std::endl;
        std::cout << "old result = " << Math::AS_5_1_56(x) / x << std::endl;
        std::cout << "new result = " << Math::_deBoerD(x, 1.0e-5) << std::endl;
        std::cout << "limit0 = " << limit0 << std::endl;
        std::cout << "limit1 = " << limit1 << std::endl;
    }
    return tmpResult;
#else
    if (x < 0)
    {
        return std::exp(x) * E1(x);
    }
    if (x > 1)
    {
        return Math::_deBoerD(x);
    }
    else
        return std::exp(x) * (Math::AS_5_1_53(x) - log(x));
#endif
}


double Math::deBoerL0(const double & mu1, const double & mu2, const double & muj, \
                               const double & density, const double & thickness)
{
    double d;
    double tmpDouble;

    /*
    if (!Math::isFiniteNumber(mu1))
    {
        std::cout << "mu1 = " << mu1 << std::endl;
        throw std::runtime_error("Math::deBoerL0. Received not finite mu1 < 0");
    }
    if (!Math::isFiniteNumber(mu2))
    {
        std::cout << "mu2 = " << mu2 << std::endl;
        throw std::runtime_error("Math::deBoerL0. Received not finite mu2 < 0");
    }
    if (!Math::isFiniteNumber(muj))
    {
        std::cout << "muj = " << muj << std::endl;
        throw std::runtime_error("Math::deBoerL0. Received non finite muj < 0");
    }
    */

    // express the thickness in g/cm2
    d = thickness * density;
    if (((mu1 + mu2) * d) > 10.)
    {
        // thick target
        tmpDouble = (muj/mu1) * std::log(1 + mu1/muj) / ((mu1 + mu2) * muj);
        //std::cout << "THICK TARGET = " << tmpDouble << std::endl;
        return tmpDouble;
    }
    // std::cout << " (mu1 + mu2) * d = " << (mu1 + mu2) * d << std::endl;
    if (((mu1 + mu2) * d) < 0.01)
    {
        // very thin target, neglect enhancement
        //std::cout << "Very thin target, not considered = " << 0.0 << std::endl;
        return 0.0;
    }

    /*
    if ((mu1*d < 0.1) && (muj * d < 1) )
    {
        // thin target (this approximation only gives the order of magnitude.
        // it is not as good as the thick target one
        // std::cout << " d = " << d << " muj * d " << muj * d << " ";
        tmpDouble = -0.5 * (muj * d) * std::log(muj * d) / ((mu1 + mu2) * muj);
        tmpDouble *= (1.0 - exp(-(mu1 + mu2) * d));
        std::cout << "THIN TARGET = " << tmpDouble << std::endl;
        return tmpDouble;
    }
    */

    tmpDouble = Math::deBoerD((muj - mu2) * d) / (mu2 * (mu1 + mu2));
    tmpDouble = tmpDouble -(Math::deBoerD(muj * d) / (mu1 * mu2)) + \
                 (Math::deBoerD((muj + mu1) * d) / (mu1 * (mu1 + mu2)));
    tmpDouble *= std::exp(-(mu1 + muj) * d);

    tmpDouble += std::log(1.0 + (mu1/muj)) / (mu1 * (mu1 + mu2));

    if (mu2 < muj)
    {
        tmpDouble += (std::exp(-(mu1 + mu2) * d) / (mu2 * (mu1 + mu2))) * \
                      std::log(1.0 - (mu2 / muj));
    }
    else
    {
        tmpDouble += (std::exp(-(mu1 + mu2) * d) / (mu2 * (mu1 + mu2))) * \
                      std::log((mu2 / muj) - 1.0);
    }
    if (tmpDouble < 0)
    {
        std::cout << "CALCULATED = " << tmpDouble << std::endl;
        std::cout << " mu1 = " << mu1 << std::endl;
        std::cout << " mu2 = " << mu2 << std::endl;
        std::cout << " muj = " << muj << std::endl;
        std::cout << " d = " << d << std::endl;
    }
    return tmpDouble;
}

double Math::deBoerX(const double & p, const double & q, const double & d1, const double & d2, \
                     const double & mu_1_j, const double & mu_2_j, const double & mu_b_j_d_t)
{
    return Math::deBoerV(p, q, d1, d2, mu_1_j, mu_2_j, mu_b_j_d_t) - \
           Math::deBoerV(p, q, d1, 0.0, mu_1_j, mu_2_j, mu_b_j_d_t) - \
           Math::deBoerV(p, q, 0.0, d2, mu_1_j, mu_2_j, mu_b_j_d_t) + \
           Math::deBoerV(p, q, 0.0, 0.0, mu_1_j, mu_2_j, mu_b_j_d_t);
}

double Math::deBoerV(const double & p, const double & q, const double & d1, const double & d2, \
                     const double & mu1j, const double & mu2j, const double & mubjdt)
{
    double tmpDouble1;
    double tmpDouble2;

    // case V(0,0) with db equal to 0
    if ((mubjdt == 0) && (p == 0) && (q == 0))
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


    tmpDouble1 = Math::deBoerD((1.0 + (p / mu2j)) * (mu1j*d1 + mubjdt + mu2j*d2));
    tmpDouble1 *= (mu2j /(p * (p * mu1j + q * mu2j)));
    tmpDouble1 *= std::exp((q - mu1j) * d1 - (p + mu2j) * d2 -mubjdt);

    tmpDouble2 = mu1j * d1 + mubjdt + mu2j * d2;
    tmpDouble2 = ((mu1j / (q * (p * mu1j + q * mu2j))) * Math::deBoerD(( 1.0 - (q/mu1j)) * tmpDouble2)) - \
                 (Math::deBoerD(tmpDouble2)/(p * q));

    return tmpDouble1 + tmpDouble2;
}

bool Math::isNumber(const double & x)
{
    return (x == x);
}

bool Math::isFiniteNumber(const double & x)
{
    return (x <= DBL_MAX && x >= -DBL_MAX);
}

double Math::_deBoerD(const double &x, const double & epsilon, const int & maxIter)
{
    // Evaluate exp(x) * E1(x) for x > 1
    //
    // Adapted from continued fraction expression of En(x) from Mathematica wb site
    //
    // Modified Lentz algorithm following Numerical Recipes description
    //
    double f, D, C;
    // double tiny = 1.0e-30; not needed, we never get 0 denominator.
    double a, b, delta;

    if (x <= 1)
    {
        std::cout << "x = " << x << std::endl;
        throw std::runtime_error("_deBoerD algorithm converges for x > 1");
    }

    // In the Lentz algorithm, we have to provide b0, b(i) and a(i) for i = 1, ...
    // The rest of the algorithm is the same for each function
    //b = 1 + x;  // b0
    b = 1 + x;
    f = b;          // f = b0
    C = f;          // C = f0
    D = 0.0;        // D = 0.0
    for (int i = 1; i < maxIter; i++)
    {
        b = b + 2;      // b(i) = x + 2 * i + 1;
        a = - i * i;    // a(i) = - (i * i)
        C = b + a / C;              // C(i) = b(i) + a(i) / C(i-1)
        //if (C == 0)   // This check is not needed
        //    C = tiny;
        D = b + a * D;  // D(i) = b(i) - a(i) * D(i-1)
        //if (D == 0)   // This check is not needed
        //    D = tiny;
        D = 1.0 / D;
        delta = C * D;
        f *= delta;
        if (std::abs(delta - 1) < epsilon)
        {
            // The continued fraction is already summed in f
            // adapt to what we want to calculate
            return 1.0 / f;
        }
    }

    std::cout << " Continued fraction failed to converge for x = " << x << std::endl;
    // return average of quoted values
    double limit0, limit1;
    limit0 = 0.5 * log(1 + 2.0/x);
    limit1 = log(1 + 1.0 /x);
    return 0.5 * (limit0 + limit1);
}
