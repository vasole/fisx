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
    const std::vector<std:: string> & getKeys();
    const std::map<std::string, std::string> & readKey(const std::string & key,
                                                       const std::string & defaultValue="");
    template<typename T>
    void parseStringAsSingleValue(const std::string & keyContent,\
                                  T & destination,
                                  const T & defaultValue);

    template<typename T>
    static void parseStringAsMultipleValues(const std::string & keyContent,
                                            std::vector<T> & destination, const T & defaultValue)
    {
        std::stringstream ss(keyContent);
        T result;
        std::string item;
        char delimiter;

        delimiter = ',';
        destination.clear();
        while (std::getline(ss, item, delimiter))
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
    std::map<std::string, std::map<std::string, std::string> > keyContents;
    std::vector<std::string> keys;
    std::map<std::string, long> keyPositions;
    std::map<std::string, std::string>  defaultContent;
};
#endif
