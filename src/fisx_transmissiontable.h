#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2020 European Synchrotron Radiation Facility
#
# This file is part of the fisx X-ray library developed by V.A. Sole
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
#ifndef FISX_TRANSMISSIONTABLE_H
#define FISX_TRANSMISSIONTABLE_H
#include <string>
#include <vector>
#include <map>

namespace fisx
{

/*!
  \class TrasmissionTable
  \brief Class containing a transmission table

   A transmission table is a tabulated list of energies (in keV) and values between 0 and 1 to represent a transmissio.

   It can be used to take into account the transmission of optical elements or windows as provided by manufacturers.

*/
class TransmissionTable
{
public:
    /*!
    Minimalist constructor.
    */
    TransmissionTable();

    /*!
    Expected constructor.

    Just a name and a description.
    */
    TransmissionTable(const std::string & name, const std::string & comment = "");

    void setName(const std::string & name);

    void setComment(const std::string & comment);

    /*!
    Provide the table as a map of doubles.
    */
    void setTransmissionTable(const std::map<double, double> & table, \
                              const std::string & name="", const std::string & comment = "");

    /*!
    Alternative method to provide the table via two lists of doubles
    */
    void setTransmissionTable(const std::vector<double> & energy, const std::vector<double> & transmission, \
                              const std::string & name="", const std::string & comment = "");

    /*!
    Retrieves the internal table of energies and associated transmission
    */
    std::map<double, double> getTransmissionTable() const;

    /*!
    Return the transmission at a given energy by log-log interpolation into the internal table
    */
    double getTransmission(const double &energy) const;

    /*!
    Return the transmission at a set of energy by linear interpolation into the internal table
    */
    std::vector<double> getTransmission(const std::vector<double> & energy) const;

    std::string getName() const;
    std::string getComment() const;

private:
    std::string name;
    std::string comment;
    std::map<double, double> transmissionTable;
};

} // namespace fisx

#endif //FISX_TRANSMISSIONTABLE_H
