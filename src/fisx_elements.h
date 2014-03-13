#ifndef FISX_ELEMENTS_H
#define FISX_ELEMENTS_H
/*!
  \file fisx_elements.h
  \brief Elements properties
  \author V.A. Sole
  \version 0.1
 */

#include <string>
#include <vector>
#include <map>
#include "fisx_simplespecfile.h"
#include "fisx_element.h"
#include "fisx_epdl97.h"
#include "fisx_material.h"

/*!
  \class Elements
  \brief Class handling the physical properties

   This class initializes a default library of physical properties and allows the user
   to modify and to access those properties.
 */
class Elements
{

public:
    /*!
    Initialize the library from the EPDL97 data files found in the provided directory.
    */
    Elements(std::string dataDirectory);

    /*!
    Initialize the library from the EPDL97 data files found in the provided directory.
    It forces the EPDL97 photon photoelectric cross sections to follow the shell binding
    energies found in the provided binding energies file (full path needed because it is
    not necessarily found in the the same directory as the EPDL97 data).
    */
    Elements(std::string dataDirectory, std::string bindingEnergiesFile);

    /*!
    Convenience constructor. It is equivalent to use the previous constructor and
    - to call the setMassAttenuationCoefficientsFile with that file name.
    - to call the setMassAttenuationCoefficientes method with the values extracted from the crossSectionFile.
    */
    Elements(std::string dataDirectory, std::string bindingEnergiesFile, std::string crossSectionsFile);

    // Direct element handling
    /*!
    Returns true if the element with name elementName is already defined in the library.
    */
    bool isElementNameDefined(const std::string & elementName);
    /*!
    Returns a reference to the element with name elementName if defined in the library.
    */
    const Element & getElement(const std::string & elementName);
    /*!
    Get a copy of the element with name elementName.
    */
    Element getElementCopy(const std::string & elementName);

    // function to ADD or REPLACE if already existing an element
    /*!
    Add an element instance to the library.
    */
    void addElement(Element & elementInstance);

    /*!
    Retrieve the names of the elements already defined in the library.
    */
    std::vector<std::string> getElementNames();

    // mass attenuation related functions
    /*!
    Update the total mass attenuation coefficients of the default elements with those found
    in the given file.
    */
    void setMassAttenuationCoefficientsFile(std::string fileName);

    /*!
    Update the total mass attenuation coefficients of the supplied element.
    The partial mass attenuation photoelectric coefficients are updated by the library in order
    to be consistent with the supplied mass attenuation coefficients
    */
    void setMassAttenuationCoefficients(std::string elementName, \
                                        std::vector<double> energy, \
                                        std::vector<double> photoelectric, \
                                        std::vector<double> coherent, \
                                        std::vector<double> compton,\
                                        std::vector<double> pair);

    /*!
    Retrieve the internal table of photon mass attenuation coefficients of the requested element.
    */
    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients(std::string elementName);

    /*!
    Given an element and a set of energies, give back the mass attenuation coefficients at the given
    energies as a map where the keys are the different physical processes and the values are the vectors
    of the calculated values via log-log interpolation in the internal table.
    */
    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients(std::string elementName,
                                                                         std::vector<double> energies);

    /*!
    Given a map of elements and mass fractions element and a set of energies, give back the mass attenuation
    coefficients at the given energies as a map where the keys are the different physical processes and the
    values are the vectors of the calculated values via log-log interpolation in the internal table.
    */
    std::map<std::string, std::vector<double> > getMassAttenuationCoefficients(\
                                                std::map<std::string, double> elementMassFractions,\
                                                std::vector<double> energies);


    /*!
    Convenience method.
    Given an element and an energy, give back the mass attenuation coefficients at the given
    energy as a map where the keys are the different physical processes and the values are calculated
    values via log-log interpolation in the internal table.
    */
    std::map<std::string, double> getMassAttenuationCoefficients(std::string elementName, double energy);

    /*!
    Convenience method.
    Given a map of elements and mass fractions element and one energy, give back the mass attenuation
    coefficients at the given energy as a map where the keys are the different physical processes and the
    values are the calculated values via log-log interpolation in the internal table.    */
    std::map<std::string, double> getMassAttenuationCoefficients(\
                                                std::map<std::string, double> elementMassFractions,\
                                                double energies);

    // Material handling
    /*!
    Greate a new Material with the given name.
    */
    void createMaterial(std::string name);

    /*!
    Greate a new Material with the given name, density, thickness and comment.
    */
    void createMaterial(const std::string & name, const double & density,
                        const double & thickness, const std::string  & comment);

    /*!
    Set the material composition of the material with name materialName.
    A composition is a map where the keys are elements/materials already defined in the library
    and the values are mass amounts.
    The library is supposed to normalize to make sure the mass fractions sum unity.
    */
    void setMaterialComposition(const std::string & materialName, \
                                const std::map<std::string, double> & composition);

    /*!
    Set the material composition of the material with name materialName.
    It is composed of the elements/materials with the given names and mass amounts.
    The library is supposed to normalize to make sure the mass fractions sum unity.
    */
    void setMaterialComposition(const std::string & materialName, \
                                const std::vector<std::string> & names,\
                                const std::vector<double> & amounts);

    /*!
    Retrieve a reference to instance of the material identified by materialName
    */
    const Material & getMaterial(const std::string & materialName);

    /*!
    Retrieve a copy of the instance of the material identified by materialName
    */
    Material getMaterialCopy(const std::string & materialName);

    /*!
    Add a Material instance to the set of defined materials.
    */
    void addMaterial(Material & materialInstance);

    /*!
    Try to interpret a given string as a formula, returning the associated mass fractions
    as a map of elements and mass fractions. In case of failure, it returns an empty map.
    */
    std::map<std::string, double> getCompositionFromFormula(const std::string & formula);

    /*!
    Try to parse a given string as a formula, returning the associated number of "atoms"
    per single molecule. In case of failure, it returns an empty map.
    */
    std::map<std::string, double> parseFormula(const std::string & formula);

    /*!
    Utility to convert from string to double.
    */
    bool StringToDouble(const std::string & str, double& number);

    /*!
    Utility to convert from double to string.
    */
    std::string toString(const double& number);

private:
    void initialize(std::string, std::string);

    // The EPDL97 library
    EPDL97 epdl97;

    // The map of defined elements
    // The map has the form Element Name, Index
    std::map<std::string , int> elementDict;

    // The vector of defined elements
    std::vector<Element> elementList;

    // The vector of defined Materials
    std::vector<Material> materialList;

    // The map Material name, index
    std::map<std::string , int> materialDict;

};
#endif
