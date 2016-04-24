#/*##########################################################################
#
# The fisx library for X-Ray Fluorescence
#
# Copyright (c) 2014-2016 European Synchrotron Radiation Facility
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
#ifndef FISX_SIMPLE_SPECFILE_H
#define FISX_SIMPLE_SPECFILE_H
#include <string>
#include <vector>
#include <map>

namespace fisx
{

#ifdef SCAN_CLASS
class SpecfileScan
{
public:
    SpecfileScan(int number);
    std::vector<std::string> getHeader();
    std::vector<std::string> getHeader(std::string);
    std::vector<std::string> getAllLabels();
    std::map<std::string, std::vector<double>> getData();
    std::vector<double> getDataColumn(int column);
    std::vector<double> getDataRow(int row);

private:
    std::vector<std::string> scanBuffer;
};
#endif // SCAN_CLASS


class SimpleSpecfile
{
public:
    SimpleSpecfile();
    SimpleSpecfile(std::string fileName);
    void setFileName(std::string fileName);
    int getNumberOfScans();
    std::vector<std::string> getScanHeader(int scanIndex);
    std::vector<std::string> getScanLabels(int scanIndex);
    // std::map<std::string, std::vector<double>> getScanData(int scanIndex);
    std::vector<std::vector<double> > getScanData(int scanIndex);
    //std::vector<double> getScanDataColumn(int scanIndex, std::string label);
    //std::vector<double> getScanDataColumn(int scanIndex, int column);
    //std::vector<double> getScanDataRow(int scanIndex, int row);
    // it is the responsibility of the caller to delete the scan
    // SpecfileScan* getScan(int scanNumber);
private:
    std::string fileName;
    // the starting and ending points of each scan
    std::vector<std::pair<long, long> > scanFilePositions;
    std::vector<std::streampos> scanPosition;
};

} // namespace fisx

#endif // FISX_SIMPLE_SPECFILE_H
