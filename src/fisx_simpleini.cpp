#include "fisx_simpleini.h"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stdlib.h>
#include <cctype>
#define isNumber(x) ( isdigit(x) || x == '-' || x == '+' || x == '.' || x == 'E' || x == 'e')

SimpleIni::SimpleIni()
{
    this->fileName = "";
    this->keys.clear();
    this->keyPositions.clear();
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

    this->keys.clear();
    this->keyPositions.clear();
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
                this->keys.push_back(mainKey);
                this->keyPositions[mainKey] = numberOfLines;
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
                this->keyContents[mainKey][key] = content;
            }
            else
            {
                // continuation line
                if ((key.size() > 0) && (mainKey.size() > 0))
                {
                    this->keyContents[mainKey][key] += tmpString;
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


const std::vector<std::string> & SimpleIni::getKeys()
{
    return this->keys;
}

const std::map<std::string, std::string > & SimpleIni::readKey(const std::string & key,
                                                               const std::string & defaultValue)
{
    if (this->keyContents.find(key) == this->keyContents.end())
    {
        this->defaultContent.clear();
        this->defaultContent[key] = defaultValue;
        return this->defaultContent;
    }
    return this->keyContents[key];
}

template<typename T>
void SimpleIni::parseStringAsSingleValue(const std::string & keyContent,
                                      T & destination,
                                      const T & defaultValue)
{
    std::stringstream stream(keyContent);
    stream >> destination;
    if (stream.fail())
    {
        destination = defaultValue;
    }
}

/*
template<typename T>
void SimpleIni::parseStringAsMultipleValues(const std::string & keyContent,
                                            std::vector<T> & destination, const T & defaultValue)
{
    std::stringstream ss(keyContent);
    T result;
    std::string item;
    char delimiter;

    delimiter = ',';
    destination.clear();
    while (std::getline(ss, item, delimiter)) {
        if (this->stringConverter(item, result))
            destination.push_back(result);
        else
            destination.push_back(defaultValue);
    }
}
template<typename T>
bool SimpleIni::stringConverter(const std::string& str, T & number)
{
    std::istringstream i(str);
    if (!(i >> number))
    {
        // Number conversion failed
        return false;
    }
    return true;
}
*/
