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
#include "fisx_simplespecfile.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stdlib.h>
#include <cctype>
#define isNumber(x) ( isdigit(x) || x == '-' || x == '+' || x == '.' || x == 'E' || x == 'e')

namespace fisx
{

SimpleSpecfile::SimpleSpecfile()
{
    this->fileName = "";
    this->scanFilePositions.clear();
    this->scanPosition.clear();
}

SimpleSpecfile::SimpleSpecfile(std::string fileName)
{
    this->setFileName(fileName);
}


void SimpleSpecfile::setFileName(std::string fileName)
{
    std::string line;
    long startLine = -1;
    long endLine = -1;
    long numberOfLines = -1;
    std::ifstream fileInstance(fileName.c_str(), std::ios::in | std::ios::binary);
    std::streampos position;

    this->scanFilePositions.clear();
    this->scanPosition.clear();
    position = 0;
    while (std::getline(fileInstance, line))
    {
        //std::cout << line << std::endl;
        ++numberOfLines;
        if(line.size() > 1)
        {
            if (line.substr(0, 2) == "#S")
            {
                startLine = numberOfLines;
                this->scanFilePositions.push_back(std::make_pair(startLine, startLine));
                endLine = -1;
                this->scanPosition.push_back(position);
            }
        }
        else
        {
            if(startLine >= 0)
            {
                startLine = -1;
                endLine = numberOfLines;
                this->scanFilePositions.back().second = endLine;
            }
        }
        position = fileInstance.tellg();
    }

    if ((endLine == -1) && startLine >= 0)
    {
        this->scanFilePositions.back().second = numberOfLines + 1;
    }
    fileInstance.clear();
    if (fileInstance.is_open())
    {
        fileInstance.close();
    }
    // std::cout << "Number of scans: " << this->scanFilePositions.size();
    // std::cout << std::endl;
    this->fileName = fileName;
    //std::cout << "PASSED" << std::endl;
}


int SimpleSpecfile::getNumberOfScans()
{
    return (int) this->scanFilePositions.size();
}


std::vector<std::string> SimpleSpecfile::getScanLabels(int scanIndex)
{
    std::ifstream fileInstance(this->fileName.c_str(),  std::ios::in | std::ios::binary);
    std::string line;
    std::string::size_type iStart, iEnd;
    long i;
    long nLines;
    std::vector<std::string> result;

    // fileInstance.seekg(0, std::ifstream::beg);
    if((scanIndex >= (long) this->getNumberOfScans()) || (scanIndex < 0))
    {
        throw std::invalid_argument("Not a valid scan index");
    }

    // This is very slow
    //for (i = 0; i < this->scanFilePositions[scanIndex].first; i++)
    //{
    //    std::getline(fileInstance, line);
    //    //std::cout << line << std::endl;
    //}
    // This is faster but needs the file open in binary mode
    fileInstance.seekg(this->scanPosition[scanIndex], std::ios::beg);
    nLines = 1 + this->scanFilePositions[scanIndex].second - \
                        this->scanFilePositions[scanIndex].first;

    if (nLines < 0)
    {
        throw std::runtime_error("Negative number of lines to be read !!!");
    }

    i = 0;
    while (i < nLines)
    {
        if (line.size() > 1)
        {
            if (line.substr(0, 2) != "#L")
            {
                std::getline(fileInstance, line);
            }
            else
            {
                // go out of the loop
                i = nLines;
            }
        }
        else
        {
            std::getline(fileInstance, line);
        }
        i++;
    }
    if (line.size() < 2)
    {
        throw std::runtime_error("Label line not found");
    }

    if (line.substr(0, 2) != "#L")
    {
        throw std::runtime_error("Label line not found");
    }

    // trim trailing CR if present    
    if (line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);

    // trim leading and trailing spaces
    iStart = line.find_first_of(" ") + 1;
    iEnd = line.find_last_not_of(" ");
    line = line.substr(iStart, 1 + iEnd - iStart);

    // split on two spaces
    i = 0;
    iStart = 0;
    while ( i < (long) (line.size() - 2))
    {
        if (line.substr(i, 2) == "  ")
        {
            result.push_back(line.substr(iStart, i - iStart));
            while ((line.substr(i, 1) == " ") && (i < (long) line.size()))
            {
                i++;
            }
            iStart = i;
        }
        else
        {
            i++;
        }
    }
    if (iStart < line.size())
    {
        result.push_back(line.substr(iStart, line.size() - iStart));
    }
    return result;
}

std::vector<std::vector<double> > SimpleSpecfile::getScanData(int scanIndex)
{
    std::ifstream fileInstance(this->fileName.c_str(),  std::ios::in | std::ios::binary);
    std::string line;
    std::string::size_type iString;
    std::string tmpString;
    long i;
    bool replaceDot;
    std::vector<std::vector<double> > result;
    std::vector<double> lineNumbers;
    std::vector<std::vector<double> >::size_type row = 0;

    if((strtod("4.5", NULL) - 4.0) < 0.4)
    {
        replaceDot = true;
    }
    else
    {
        replaceDot = false;
    }

    if((scanIndex >= (long) this->getNumberOfScans()) || (scanIndex < 0))
    {
        throw std::invalid_argument("Not a valid scan index");
    }

    //for (i = 0; i < this->scanFilePositions[scanIndex].first; i++)
    //{
    //    std::getline(fileInstance, line);
    //}

    // If instead of the loop I use this it does not work unless the file is opened in binary mode
    fileInstance.seekg(this->scanPosition[scanIndex], std::ios::beg);
    i = this->scanFilePositions[scanIndex].first;
    for ( ; i < this->scanFilePositions[scanIndex].second; i++)
    {
        std::getline(fileInstance, line);
        if (!line.empty())
        {
            if(line[line.size() - 1] == '\r')
                line.erase(line.size() - 1);
        }
        // std::cout << line << std::endl;
        if ((i < this->scanFilePositions[scanIndex].first) ||\
            (line[0] == '#'))
        {
            // either we have not reached the scan
            // or it is a header or comment line
            ;
        }
        else
        {
            // std::cout << line.size() << "Numeric Data line: " << line << std::endl;
            // parse the line
            iString = 0;
            lineNumbers.clear();
            while (iString < line.size())
            {
                tmpString.clear();
                while( (!isNumber(line[iString])) && (iString < line.size()))
                {
                    iString++;
                }
                while( (isNumber(line[iString])) && (iString < line.size()))
                {
                    if (replaceDot && (line[iString] == '.'))
                    {
                        tmpString += ",";
                    }
                    else
                    {
                        tmpString += line[iString];
                    }
                    iString++;
                }
                if (tmpString.size() > 0)
                {
                    lineNumbers.push_back(strtod(tmpString.c_str(), NULL));
                }
            }
            if (lineNumbers.size() > 0)
            {
                if (result.size() > 0)
                {
                    if (result[0].size() != lineNumbers.size())
                    {
                        throw std::runtime_error("Badly formatted line");
                    }
                }
                result.push_back(lineNumbers);
            }
            row++;
        }
    }
    return result;
}

} // namespace fisx
