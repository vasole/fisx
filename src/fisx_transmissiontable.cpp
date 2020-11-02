#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2020 European Synchrotron Radiation Facility
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
#include <stdexcept>
#include "fisx_transmissiontable.h"
#include <iostream>

namespace fisx
{

TransmissionTable::TransmissionTable()
{
    this->name = "";
    this->comment = "";
    // no attenuation at all
    this->transmissionTable[0.0] = 1.0;
}

TransmissionTable::TransmissionTable(const std::string & name, const std::string & comment)
{
    this->name = name;
    this->comment = comment;
    // no attenuation at all
    this->transmissionTable[0.0] = 1.0;
}


void TransmissionTable::setName(const std::string & name)
{
    this->name = name;
}

void TransmissionTable::setComment(const std::string & comment)
{
    this->comment = comment;
}


void TransmissionTable::setTransmissionTable(const std::map<double, double> & table, \
                          const std::string & name, const std::string & comment)
{
    std::string msg;
        throw std::invalid_argument(msg);
    std::map<double, double>::const_iterator c_it;

    // verify that no value is negative
    for (c_it = table.begin(); c_it != table.end(); ++c_it)
    {
        if (c_it->first < 0.0)
        {
            msg = "TransmissionTable::setTransmissionTable. Invalid energy";
            throw std::invalid_argument(msg);
        }
        if (c_it->second < 0.0)
        {
            msg = "TransmissionTable::setTransmissionTable. Invalid transmission";
            throw std::invalid_argument(msg);
        }
    }

    // verify the table was not already named
    if ((this->name.size() > 0) and (name.size() > 0))
    {
        msg = "TransmissionTable::setTransmissionTable. Table cannot be renamed";
        throw std::invalid_argument(msg);
    }
    if ((this->comment.size() > 0) and (comment.size() > 0))
    {
        msg = "TransmissionTable::setTransmissionTable. Table comment cannot be changed";
        throw std::invalid_argument(msg);
    }

    // assign the values
    this->transmissionTable = table;
    this->name = name;
    this->comment = comment;
}

void TransmissionTable::setTransmissionTable(const std::vector<double> & energy, \
                                             const std::vector<double> & transmission, \
                                             const std::string & name, \
                                             const std::string & comment)
{
    std::vector<double>::size_type i;
    std::string msg;

    // check input
    // verify the table was not already named
    if ((this->name.size() > 0) and (name.size() > 0))
    {
        msg = "TransmissionTable::setTransmissionTable. Table cannot be renamed";
        throw std::invalid_argument(msg);
    }
    if ((this->comment.size() > 0) and (comment.size() > 0))
    {
        msg = "TransmissionTable::setTransmissionTable. Table comment cannot be changed";
        throw std::invalid_argument(msg);
    }

    // verify the data are consistent
    if (energy.size() != transmission.size())
    {
        msg = "TransmissionTable::setTransmissionTable. Number of energies and attenuations not matching";
        throw std::invalid_argument(msg);
    }
    for (i = 0; i < energy.size(); ++i)
    {
        if (energy[i] < 0.0)
        {
            msg = "TransmissionTable::setTransmissionTable. Negative energy";
            throw std::invalid_argument(msg);
        }
        if (transmission[i] < 0.0)
        {
            msg = "TransmissionTable::setTransmissionTable. Negative transmission";
            throw std::invalid_argument(msg);
        }
    }

    // fill the table
    this->transmissionTable.clear();
    for (i = 0; i < energy.size(); ++i)
    {
        this->transmissionTable[energy[i]] = transmission[i];
    }
    this->name = name;
    this->comment = comment;
}


std::string TransmissionTable::getComment() const
{
    return this->comment;
}

std::string TransmissionTable::getName() const
{
    return this->name;
}

double TransmissionTable::getTransmission(const double & energy) const
{
    //typedef std::map<double, double>::const_iterator i_t;
    std::map<double, double>::const_iterator i1, i2;
    double x1, x2, factor;

    if (this->transmissionTable.size() < 1)
    {
        return 1.0;
    }

    if (energy <= 0.0)
    {
        return 1.0;
    }

    i2 = this->transmissionTable.upper_bound(energy);
    if (i2 == this->transmissionTable.end())
    {
        // return last value of the table
        // return (--i)->second;
        // or return 1 because we are dealing with a transmission?
        return 1.0;
    }
    i1 = this->transmissionTable.begin();
    if (energy < (i1->second))
    {
        // return first value of the table
        // return i->second;
        // or return 0 because we are dealing with a transmission?
        return 0.0;
    }
    i1 = i2;
    i1--;

    x1 = i1->first;
    x2 = i2->first;

    // perform a linear interpolation
    factor = (energy - x1) / (x2 - x1);
    return factor * (i2->second) + (1.0 - factor) * i1->second;
}

std::vector<double> TransmissionTable::getTransmission(const std::vector<double> & energy) const
{
    std::vector<double>::size_type n, i;
    std::vector<double> result;

    n = energy.size();
    result.resize(n);
    for (i = 0; i < n; ++i)
    {
        result[i] = this->getTransmission(energy[i]);
    }
    return result;
}

std::map<double, double> TransmissionTable::getTransmissionTable() const
{
    return this->transmissionTable;
}



} // namespace fisx
