#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2018 European Synchrotron Radiation Facility
#
# This file is part of the fisx X-ray developed by V.A. Sole
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
#############################################################################*/
#include "fisx_math.h"
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <iostream>

namespace fisx
{

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
        //
        // E1(z) = - gamma - log(z) - sum ( (-1)^n * x^n / (n * n! ))
        //
        //
        // Euler's gamma = 0.577215664901532860;
        // I decide to evaluate just 10 terms of the series
        //
        // TODO: Unroll the loop with n * factorial(n) already computed.
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
        return Math::AS_5_1_53(x) - std::log(x);
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

double Math::En(const int & n, const double & x)
{
    if (n < 1)
    {
        throw std::runtime_error("Math::En(n, x). n Must be greater or equal to 1");
    }
    if (n == 1)
    {
        return Math::E1(x);
    }
    else
    {
        if (x == 0)
        {
            // special value
            return 1.0 / (n - 1);
        }
        else
        {
            // Use recurrence relation
            //              1
            // E(n+1, z) = ---- ( exp (-z) - z * E(n, z))
            //              n
            //
            return (std::exp(-x) - x * Math::En(n - 1, x)) / (n - 1);
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
    num = a[3] * x;
    den = b[3] * x;
    for (int i = 2; i >= 0; --i)
    {
        num = (num * x )+ a[i];
        den = (den * x )+ b[i];
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
        tmpResult = std::exp(x) * (Math::AS_5_1_53(x) - std::log(x));
    limit0 = 0.5 * std::log(1 + 2.0/x);
    limit1 = std::log(1 + 1.0 /x);
    if ((tmpResult < limit0) || (tmpResult > limit1))
    {
        std::cout << "deBoerD error with x = " << x << std::endl;
        std::cout << "old result = " << Math::AS_5_1_56(x) / x << std::endl;
        std::cout << "new result = " << Math::_deBoerD(x, 1.0e-5) << std::endl;
        std::cout << "limit0 = " << limit0 << std::endl;
        std::cout << "limit1 = " << limit1 << std::endl;
        return Math::_deBoerD(x, 1.0e-5);
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
    if ((mu1 <= 0.0) || (mu2 <= 0.0) || (muj <= 0.0))
    {
        std::cout << "mu1 = " << mu1 << std::endl;
        std::cout << "mu2 = " << mu2 << std::endl;
        std::cout << "muj = " << muj << std::endl;
        throw std::runtime_error("Math::deBoerL0 received negative input");
    }

    // express the thickness in g/cm2
    d = thickness * density;
    if (((mu1 + mu2) * d) > 10.)
    {
        // thick target
        tmpDouble = (muj/mu1) * std::log(1 + mu1/muj) / ((mu1 + mu2) * muj);
        //std::cout << "THICK TARGET = " << tmpDouble << std::endl;
        if (!Math::isFiniteNumber(tmpDouble))
        {
            std::cout << "Math::deBoerL0. Thick target. Not a finite result" << std::endl;
            std::cout << "Received parameters " << std::endl;
            std::cout << "mu1 = " << mu1 << std::endl;
            std::cout << "mu2 = " << mu2 << std::endl;
            std::cout << "muj = " << muj << std::endl;
            std::cout << "thickness = " << thickness << std::endl;
            std::cout << "density = " << density << std::endl;
            throw std::runtime_error("Math::deBoerL0. Thick target. Non-finite result");
        }
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
        std::cout << " Math::deBoerL0 CALCULATED = " << tmpDouble << std::endl;
        std::cout << " mu1 = " << mu1 << std::endl;
        std::cout << " mu2 = " << mu2 << std::endl;
        std::cout << " muj = " << muj << std::endl;
        std::cout << " d = " << d << std::endl;
        throw std::runtime_error("Math::deBoerL0. Negative result");
    }
    if (!Math::isFiniteNumber(tmpDouble))
    {
        std::cout << " Math::deBoerL0 CALCULATED = " << tmpDouble << std::endl;
        std::cout << " mu1 = " << mu1 << std::endl;
        std::cout << " mu2 = " << mu2 << std::endl;
        std::cout << " muj = " << muj << std::endl;
        std::cout << " d = " << d << std::endl;
        throw std::runtime_error("Math::deBoerL0. Non-finite result");
    }
    return tmpDouble;
}

double Math::deBoerX(const double & p, const double & q, const double & d1, const double & d2, \
                     const double & mu1j, const double & mu2j, const double & mubj_dt)
{
     if ((false))
     {
        double result;
        result =  Math::deBoerV(p, q, d1, d2, mu1j, mu2j, mubj_dt) - \
               Math::deBoerV(p, q, d1, 0.0, mu1j, mu2j, mubj_dt) - \
               Math::deBoerV(p, q, 0.0, d2, mu1j, mu2j, mubj_dt) + \
               Math::deBoerV(p, q, 0.0, 0.0, mu1j, mu2j, mubj_dt);
        if (d1 < 0.01)
        {
            // VERY THIN CASE:
            std::cout << " THIN CASE " << std::endl;
            std::cout << " EXPECTED = " << (d1 / p) * std::log(1 + p / mu2j) << std::endl;
            std::cout << " MEASURED = " << result << std::endl;
            std::cout << " V(d1,inf) = " << Math::deBoerV(p, q, d1, d2, mu1j, mu2j, mubj_dt) << std::endl;
            std::cout << " V(d1, 0) = " << Math::deBoerV(p, q, d1, 0.0, mu1j, mu2j, mubj_dt) << std::endl;
            std::cout << " V(0, d2) = " << Math::deBoerV(p, q, 0.0, d2, mu1j, mu2j, mubj_dt) << std::endl;
            std::cout << " V(0.0, 0) = " << Math::deBoerV(p, q, 0.0, 0.0, mu1j, mu2j, mubj_dt) << std::endl;
        }
        else
        {
            std::cout << "NOT THIN CASE " << std::endl;
        }
        if (result < 0)
        {
            std::cout << "p    " << p << std::endl;
            std::cout << "q    " << q << std::endl;
            std::cout << "d1   " << d1 << std::endl;
            std::cout << "d2   " << d2 << std::endl;
            std::cout << "mu1j " << mu1j << std::endl;
            std::cout << "mu2j " << mu2j << std::endl;
            std::cout << "mubjdt " << mubj_dt << std::endl;
            std::cout << " error  " << std::endl;
            throw std::runtime_error("negative contribution");
        }
        if (!Math::isFiniteNumber(result))
        {
            std::cout << "p    " << p << std::endl;
            std::cout << "q    " << q << std::endl;
            std::cout << "d1   " << d1 << std::endl;
            std::cout << "d2   " << d2 << std::endl;
            std::cout << "mu1j " << mu1j << std::endl;
            std::cout << "mu2j " << mu2j << std::endl;
            std::cout << "mubjdt " << mubj_dt << std::endl;
            std::cout << " error  " << std::endl;
            std::cout << "Math::deBoerV(p, q, d1, d2, mu1j, mu2j, mubj_dt) = ";
            std::cout << Math::deBoerV(p, q, d1, d2, mu1j, mu2j, mubj_dt) << std::endl;
            std::cout << "Math::deBoerV(p, q, d1, 0.0, mu1j, mu2j, mubj_dt) = ";
            std::cout << Math::deBoerV(p, q, d1, 0.0, mu1j, mu2j, mubj_dt) << std::endl;
            std::cout << "Math::deBoerV(p, q, 0.0, d2, mu1j, mu2j, mubj_dt) = ";
            std::cout << Math::deBoerV(p, q, 0.0, d2, mu1j, mu2j, mubj_dt) << std::endl;
            std::cout << "Math::deBoerV(p, q, 0.0, 0.0, mu1j, mu2j, mubj_dt) = ";
            std::cout << Math::deBoerV(p, q, 0.0, 0.0, mu1j, mu2j, mubj_dt)<< std::endl;
            throw std::runtime_error("Not finite contribution");
        }
        return result;
     }
     else
     {
        return Math::deBoerV(p, q, d1, d2, mu1j, mu2j, mubj_dt) - \
               Math::deBoerV(p, q, d1, 0.0, mu1j, mu2j, mubj_dt) - \
               Math::deBoerV(p, q, 0.0, d2, mu1j, mu2j, mubj_dt) + \
               Math::deBoerV(p, q, 0.0, 0.0, mu1j, mu2j, mubj_dt);

    }
}

double Math::deBoerV(const double & p, const double & q, const double & d1, const double & d2, \
                     const double & mu1j, const double & mu2j, const double & mubjdt)
{
    double tmpDouble1;
    double tmpDouble2;
    double tmpHelp;


    // case V(0,0) with db equal to 0
    if ((mubjdt == 0) && (d1 == 0) && (d2 == 0))
    {
        // V(0, 0) with db = 0;
        tmpDouble1 = 1.0 - (q / mu1j);
        tmpDouble2 = 1.0 + (p / mu2j);
        if (tmpDouble1 < 0)
            tmpDouble1 = - tmpDouble1;
        if (tmpDouble2 < 0)
            tmpDouble2 = - tmpDouble2;
        tmpHelp = (mu2j / p) * std::log(tmpDouble2) + (mu1j/q) * std::log(tmpDouble1);
        tmpHelp = -tmpHelp / (p * mu1j + q * mu2j);
        if (!Math::isFiniteNumber(tmpHelp))
        {
            std::cout << "p    " << p << std::endl;
            std::cout << "q    " << q << std::endl;
            std::cout << "d1   " << d1 << std::endl;
            std::cout << "d2   " << d2 << std::endl;
            std::cout << "mu1j " << mu1j << std::endl;
            std::cout << "mu2j " << mu2j << std::endl;
            std::cout << "mubjdt " << mubjdt << std::endl;
            std::cout << "1.0 + (p / mu2j) = " << 1.0 + (p / mu2j) << std::endl;
            std::cout << "tmpDouble1 = " << tmpDouble1 << std::endl;
            std::cout << "tmpDouble2 = " << tmpDouble2 << std::endl;
            std::cout << "p * mu1j + q * mu2j = " << p * mu1j + q * mu2j << std::endl;
            std::cout << "Error 0" << std::endl;
            throw std::runtime_error("Error 0: Error on V(0,0) with no intermediate layer");
        }
        return tmpHelp;
    }

    // X(p, q, d1, d2) = V(d1, d2) - V(d1, 0) - V(0, d2) + V(0, 0)
    // V(inf, d2) = 0.0;
    // V(d1, inf) = 0.0; provided that q - muij < 0
    // Therefore, for a layer on thick substrate
    // X (p, q, d1, inf) = - V(d1, 0) + V(0,0)
    // and for small values of d1 (thin layer on thick substrate) that gives
    // X (p, q, d1, inf)X (p, q, d1, inf) is about (d1/p) * std::log(1.0 + (p/mu2j))
    tmpDouble1 = (mu2j /(p * (p * mu1j + q * mu2j))) * \
                 Math::deBoerD((1.0 + (p / mu2j)) * (mu1j*d1 + mubjdt + mu2j*d2));
    if (!Math::isFiniteNumber(tmpDouble1))
    {
        std::cout << "p    " << p << std::endl;
        std::cout << "q    " << q << std::endl;
        std::cout << "d1   " << d1 << std::endl;
        std::cout << "d2   " << d2 << std::endl;
        std::cout << "mu1j " << mu1j << std::endl;
        std::cout << "mu2j " << mu2j << std::endl;
        std::cout << "mubjdt " << mubjdt << std::endl;
        std::cout << " error 1 " << std::endl;
        throw std::runtime_error("error1");
    }

    tmpHelp = mu1j * d1 + mubjdt + mu2j * d2;
    tmpDouble2 = (mu1j / (q * (p * mu1j + q * mu2j))) * \
                  Math::deBoerD(( 1.0 - (q / mu1j)) * tmpHelp);
    if (!Math::isFiniteNumber(tmpDouble2))
    {
        std::cout << "p    " << p << std::endl;
        std::cout << "q    " << q << std::endl;
        std::cout << "d1   " << d1 << std::endl;
        std::cout << "d2   " << d2 << std::endl;
        std::cout << "mu1j " << mu1j << std::endl;
        std::cout << "mu2j " << mu2j << std::endl;
        std::cout << "mubjdt " << mubjdt << std::endl;
        std::cout << " error 3 " << std::endl;
        throw std::runtime_error("error3");
    }

    tmpDouble2 -= Math::deBoerD(tmpHelp)/(p * q);
    if (!Math::isFiniteNumber(tmpDouble2))
    {
        std::cout << "p    " << p << std::endl;
        std::cout << "q    " << q << std::endl;
        std::cout << "d1   " << d1 << std::endl;
        std::cout << "d2   " << d2 << std::endl;
        std::cout << "mu1j " << mu1j << std::endl;
        std::cout << "mu2j " << mu2j << std::endl;
        std::cout << "mubjdt " << mubjdt << std::endl;
        std::cout << " error 4 " << std::endl;
        throw std::runtime_error("error4");
    }
    tmpHelp = std::exp((q - mu1j) * d1 - (p + mu2j) * d2 - mubjdt) * (tmpDouble1 + tmpDouble2);
    if (!Math::isFiniteNumber(tmpHelp))
    {
        std::cout << "p    " << p << std::endl;
        std::cout << "q    " << q << std::endl;
        std::cout << "d1   " << d1 << std::endl;
        std::cout << "d2   " << d2 << std::endl;
        std::cout << "mu1j " << mu1j << std::endl;
        std::cout << "mu2j " << mu2j << std::endl;
        std::cout << "mubjdt " << mubjdt << std::endl;
        std::cout << "(q - mu1j) * d1 - (p + mu2j) * d2 - mubjdt = " ;
        std::cout << (q - mu1j) * d1 - (p + mu2j) * d2 - mubjdt << std::endl;
        std::cout << "exp((q - mu1j) * d1 - (p + mu2j) * d2 - mubjdt) = " ;
        std::cout << exp((q - mu1j) * d1 - (p + mu2j) * d2 - mubjdt) << std::endl;
        std::cout << " error 5 " << std::endl;
        throw std::runtime_error("error5");
    }
    return tmpHelp;
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

double Math::getFWHM(const double & energy, const double & noise, \
                     const double & fano, const double & quantumEnergy)
{
    return sqrt(noise * noise + energy * fano * 2.3548 * 2.3548 * quantumEnergy);
}


double Math::erf(const double & x)
{
    double z;
    double t;
    double r;

    z = std::fabs(x);
    t = 1.0 / (1.0 + 0.5 * z);
    r = t * std::exp(- z * z - 1.26551223 + t * (1.00002368 + t * (0.3740916 + \
        t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + \
        t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    if (x < 0)
    {
       r = 2.0 - r;
    }
    return (1.0 - r);
}

double Math::erfc(const double & x)
{
    double z;
    double t;
    double r;

    z = std::fabs(x);
    t = 1.0 / (1.0 + 0.5 * z);
    r = t * std::exp(- z * z - 1.26551223 + t * (1.00002368 + t * (0.3740916 + \
        t * (0.09678418 + t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + \
        t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    if (x < 0)
    {
       r = 2.0 - r;
    }
    return (r);
}

} // namespace fisx
