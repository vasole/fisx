#ifndef FISX_SIMPLE_INI_H
#define FISX_SIMPLE_INI_H
#include <string>
#include <vector>
#include <map>
#include <sstream>

class SimpleIni
{
public:
    SimpleIni();
    SimpleIni(std::string fileName);
    void readFileName(std::string fileName);
    /*!
    Get all the section names in the file.
    */
    const std::vector<std:: string> & getSections();

    /*!
    Read a particular section with the option to be case sensitive or not.
    It returns a map<string, string> with the key and the key content.
    Attention: subsections are not considered keys.
    If the section is not present, it returns an empty map.
    */
    const std::map<std::string, std::string> & readSection(const std::string & section,
                                                       const bool & caseSensitive = true);

    /*!
    Static method to parse a string
    */
    template<typename T>
    static void parseStringAsSingleValue(const std::string & keyContent,\
                                  T & destination,
                                  const T & defaultValue)
    {
        std::stringstream stream(keyContent);
        stream >> destination;
        if (stream.fail())
        {
            destination = defaultValue;
        }
    };

    template<typename T>
    static void parseStringAsMultipleValues(const std::string & keyContent,
                                            std::vector<T> & destination,
                                            const T & defaultValue,
                                            const char & separator = ',')
    {
        std::stringstream ss(keyContent);
        T result;
        std::string item;
        destination.clear();
        while (std::getline(ss, item, separator))
        {
            if (SimpleIni::stringConverter(item, result))
                destination.push_back(result);
            else
                destination.push_back(defaultValue);
        }
    };

    template<typename T>
    static bool stringConverter(const std::string& str, T & number)
    {
        std::istringstream i(str);
        if (!(i >> number))
        {
            // Number conversion failed
            return false;
        }
        return true;
    };

private:
    std::string fileName;
    std::map<std::string, std::map<std::string, std::string> > sectionContents;
    std::vector<std::string> sections;
    std::map<std::string, long> sectionPositions;
    std::map<std::string, std::string>  defaultContent;
};
#endif
