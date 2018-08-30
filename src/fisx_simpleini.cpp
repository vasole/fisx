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
#include "fisx_simpleini.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stdlib.h>

namespace fisx
{

SimpleIni::SimpleIni()
{
    this->fileName = "";
    this->sections.clear();
    this->sectionPositions.clear();
}

SimpleIni::SimpleIni(std::string fileName)
{
    this->readFileName(fileName);
}


void SimpleIni::readFileName(std::string fileName)
{
    std::string mainKey;
    std::string key;
    std::string line;
    std::string trimmedLine;
    std::string tmpString;
    std::string content;
    std::string::size_type i;
    // I should use streampos instead of size_type
    std::string::size_type p0, p1, length;
    int nQuotes;
    int nItems;
    int equalPosition;
    long numberOfLines;
    std::streampos position;
    std::ifstream fileInstance(fileName.c_str(), std::ios::in | std::ios::binary);
    std::string msg;

    if (!fileInstance.good())
    {
        msg = "File <" + fileName + "> mot accessible";
        throw std::ios_base::failure(msg);;
    }
    this->sections.clear();
    this->sectionPositions.clear();
    p0 = 0;
    position = 0;
    numberOfLines = -1;
    while (std::getline(fileInstance, line))
    {
        // std::cout << " INPUT <" << line << ">" << std::endl;
        ++numberOfLines;

        trimmedLine = "";
        if((line.size() > 0))
        {
            // trim leading and trailing spaces
            p1 = line.find_last_not_of(" \n\r\t");
            if (p1 != line.npos)
            {
                trimmedLine = line.substr(0, p1 + 1);
                p0 = trimmedLine.find_first_not_of(" \n\r\t");
                trimmedLine = trimmedLine.substr(p0);
            }
        }
        if (trimmedLine.size() > 0)
        {
            if (trimmedLine[0] == '[')
            {
                if (p0 != 0)
                {
                    msg = "Main key not as first line character";
                    std::cout << "WARNING: " << msg;
                    std::cout << " line = <" << line << ">" << std::endl;
                }
                mainKey = "";
                p0 = 0;
                p1 = 0;
                nItems = 1;
                for (i = 1; i < trimmedLine.size(); i++)
                {
                    if (trimmedLine[i] == '[')
                    {
                        msg = "Invalid line: <" + line + ">";
                        throw std::invalid_argument(msg);
                    }
                    else
                    {
                        if (trimmedLine[i] == ']')
                        {
                            nItems++;
                            p1 = i;
                        }
                    }
                }
                length = p1 - p0 - 1;
                if ((nItems != 2) || (length < 1))
                {
                    msg = "Invalid line: <" + line + ">";
                    throw std::invalid_argument(msg);
                }
                mainKey = trimmedLine.substr(p0 + 1, length);
                this->sections.push_back(mainKey);
                this->sectionPositions[mainKey] = numberOfLines;
                continue;
            }

            tmpString = "";
            nQuotes = 0;
            equalPosition = 0;
            for(i=0; i < trimmedLine.size(); i++)
            {
                if (trimmedLine[i] == '"')
                {
                    nQuotes++;
                    tmpString += trimmedLine[i];
                    continue;
                }
                if (nQuotes && (nQuotes % 2))
                {
                    // in between quotes
                    tmpString += trimmedLine[i];
                    continue;
                }
                else
                {
                    if ((trimmedLine[i] == '#') || (trimmedLine[i] == ';'))
                    {
                        i = trimmedLine.size();
                    }
                    if (equalPosition == 0)
                    {
                        if (!isspace(trimmedLine[i]))
                        {
                            tmpString += trimmedLine[i];
                        }
                    }
                    else
                    {
                        tmpString += trimmedLine[i];
                    }
                    if(trimmedLine[i] == '=')
                    {
                        if (equalPosition == 0)
                        {
                            equalPosition = tmpString.size() - 1;
                        }
                    }
                }
            }
            if (nQuotes % 2)
            {
                msg = "Unmatched double quotes in line: <" + line + ">";
                throw std::invalid_argument(msg);
            }
            // std::cout << " tmpString <" << tmpString << ">" << std::endl;

            if (tmpString.size() < 1)
            {
                // empty line
                key = "";
                continue;
            }
            if(equalPosition > 0)
            {
                // we have a key
                key = "";
                p0 = 0;
                p1 = equalPosition;
                length = p1 - p0;
                if (length < 1)
                {
                    msg = "Invalid line: <" + line + ">";
                    throw std::invalid_argument(msg);
                }
                key = tmpString.substr(p0, length);
                content = tmpString.substr(p1 + 1, tmpString.size() - p1 - 1);
                this->sectionContents[mainKey][key] = content;
            }
            else
            {
                // continuation line
                if ((key.size() > 0) && (mainKey.size() > 0))
                {
                    this->sectionContents[mainKey][key] += tmpString;
                }
                else
                {
                    std::cout << "Ignored line: <" + line + ">";
                }
            }
        }
    }
    fileInstance.clear();
    if (fileInstance.is_open())
    {
        fileInstance.close();
    }
    this->fileName = fileName;
}


const std::vector<std::string> & SimpleIni::getSections()
{
    return this->sections;
}

void SimpleIni::getSubsections(const std::string & parent, \
                               std::vector<std::string> & destination, \
                               const bool & caseSensitive)
{
    std::string targetString;
    std::string instanceString;
    std::locale loc;
    std::vector<std::string>::size_type i;

    destination.clear();
    if (parent.size() == 0)
    {
        destination.resize(this->sections.size());
        for (i = 0; i < this->sections.size(); i++)
        {
            destination[i] = this->sections[i];
        }
        return;
    }
    if (caseSensitive)
    {
        targetString = parent + ".";
        for (i = 0; i < this->sections.size(); i++)
        {
            if (this->sections[i].size() == targetString.size())
            {
                if (this->sections[i].substr(0, targetString.size()) == targetString)
                {
                    destination.push_back(this->sections[i]);
                }
            }
        }
    }
    else
    {
        targetString = parent + ".";
        this->toUpper(targetString, loc);
        for (i = 0; i < this->sections.size(); i++)
        {
            instanceString = this->sections[i];
            if (instanceString.size() >= targetString.size())
            {
                this->toUpper(instanceString, loc);
                if (instanceString.substr(0, targetString.size()) == targetString)
                {
                    destination.push_back(this->sections[i]);
                }
            }
        }
    }
}

const std::map<std::string, std::string > & SimpleIni::readSection(const std::string & key, \
                                                                   const bool & caseSensitive)
{
    std::string inputKey;
    std::string tmpKey;
    std::locale loc;
    std::vector<std::string>::size_type i, j;

    if (this->sectionContents.find(key) == this->sectionContents.end())
    {
        if (!caseSensitive)
        {
            inputKey = key;
            for (i = 0; i < key.size(); i++)
            {
                inputKey[i] = std::toupper(inputKey[i], loc);
            }
            for (i = 0; i < this->sections.size(); i++)
            {
                tmpKey = this->sections[i];
                if (tmpKey.size() == inputKey.size())
                {
                    j = 0;
                    while (std::toupper(tmpKey[j], loc) == inputKey[j])
                    {
                        j++;
                        if (j == inputKey.size())
                        {
                            break;
                        }
                    }
                    if (j < inputKey.size())
                    {
                        return this->sectionContents[tmpKey];
                    }
                }
            }
        }
        this->defaultContent.clear();
        return this->defaultContent;
    }
    return this->sectionContents[key];
}

} // namespace fisx
