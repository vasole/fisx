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
    std::string tmpString;
    std::string content;
    std::string::size_type i;
    // I should use streampos instead of size_type
    std::string::size_type p0, p1, length;
    int nQuotes;
    int nItems;
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
        //std::cout << line << std::endl;
        ++numberOfLines;

        if (line.size() > 1)
        {
            if (line[0] == '[')
            {
                mainKey = "";
                p0 = 0;
                p1 = 0;
                nItems = 1;
                for (i = 1; i < line.size(); i++)
                {
                    if (line[i] == '[')
                    {
                        msg = "Invalid line: <" + line + ">";
                        throw std::invalid_argument(msg);
                    }
                    else
                    {
                        if (line[i] == ']')
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
                mainKey = line.substr(p0 + 1, length);
                this->keys.push_back(mainKey);
                this->keyPositions[mainKey] = numberOfLines;
                continue;
            }
            tmpString = "";
            nQuotes = 0;
            for(i=0; i < line.size(); i++)
            {
                if (line[i] == '#')
                {
                    i = line.size();
                }
                if (line[i] == '"')
                {
                    nQuotes++;
                    continue;
                }
                if (nQuotes == 0)
                {
                    if ((line[i] != ' ') && (line[i] != '\r'))
                    {
                        tmpString += line[i];
                    }
                }
                else
                {
                    if(line[i] != '\r')
                    {
                        tmpString += line[i];
                    }
                }
            }
            if ((nQuotes != 0) && (nQuotes != 2))
            {
                msg = "Unmatched double quotes in line: <" + line + ">";
                throw std::invalid_argument(msg);
            }
            if (tmpString.size() < 1)
            {
                // empty line
                key = "";
                continue;
            }
            if(tmpString.find('=') != tmpString.npos)
            {
                // we have a key
                key = "";
                p0 = 0;
                p1 = tmpString.find('=');
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
