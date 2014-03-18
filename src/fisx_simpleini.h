#ifndef FISX_SIMPLE_INI_H
#define FISX_SIMPLE_INI_H
#include <string>
#include <vector>
#include <map>

class SimpleIni
{
public:
    SimpleIni();
    SimpleIni(std::string fileName);
    void readFileName(std::string fileName);
    const std::vector<std:: string> & getKeys();
    const std::map<std::string, std::string> & readKey(const std::string & key,
                                                       const std::string & defaultValue="");
private:
    std::string fileName;
    std::map<std::string, std::map<std::string, std::string> > keyContents;
    std::vector<std::string> keys;
    std::map<std::string, long> keyPositions;
    std::map<std::string, std::string>  defaultContent;
};
#endif
