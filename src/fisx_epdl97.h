#ifndef FISX_EPDL97_H
#define FISX_EPDL97_H
#include <string>
#include <ctype.h>
#include <vector>
#include <map>

class EPDL97
{
public:
    EPDL97();
    EPDL97(std::string directoryName);

    void setDataDirectory(std::string directoryName);
    // possibility to change binding energies
    void loadBindingEnergies(std::string fileName);
    void setBindingEnergies(const int & z, const std::map<std::string, double> & bindingEnergies);
    const std::map<std::string, double> & getBindingEnergies(const int & z);

    // the actual mass attenuation related functions
    std::map<std::string, double> getMassAttenuationCoefficients(const int & z, const double & energy);
    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients(const int & z,\
                                                const std::vector<double> & energy);

    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients(const int & z);

    // the vacancy distribution related functions
    std::map<std::string, double> getPhotoelectricWeights(const int & z, \
                                                          const double & energy);

    std::map<std::string, std::vector<double> > getPhotoelectricWeights(const int & z, \
                                                const std::vector<double> & energy);

    // utility functions
    std::string toUpperCaseString(const std::string &);
    std::pair<long, long> getInterpolationIndices(const std::vector<double> &,  const double &);

private:
    // internal function to load the data
    bool initialized;
    void loadData(std::string directoryName);
    void loadCrossSections(std::string fileName);

    // The directory name
    std::string directoryName;

    // The used file for binding energies
    std::string bindingEnergiesFile;

    // The used file for cross sections
    std::string crossSectionsFile;

    // The table containing all binding energies for all elements and shells
    // bindingEnergy[Z - 1]["K"] gives the binding energy for the K shell of element
    //                       with atomic number Z.
    std::vector<std::map<std::string, double> > bindingEnergy;

    // Mass attenuation data as read from the files
    // We have a table for each element but all of them share the same
    // file header structure
    std::vector<std::string> muInputLabels;
    std::map<std::string, int> muLabelToIndex;
    std::vector<std::vector<std::vector <double> > > muInputValues;
    std::vector<std::vector<double> > muEnergy;

    // Partial photoelectric mass attenuation coefficients
    // For each shell (= key), there is a vector for the energies
    // and a vector for the value of the mass attenuation coefficients
    // Expected map key values are:
    // K, L1, L2, L3, M1, M2, M3, M4, M5, "all other"
    void initPartialPhotoelectricCoefficients();
};
#endif
// FISX_EPDL97_H
