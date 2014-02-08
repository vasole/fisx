#ifndef FISX_SIMPLE_SPECFILE_H
#define FISX_SIMPLE_SPECFILE_H
#include <string>
#include <vector>
#include <map>

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
#endif


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
#endif
