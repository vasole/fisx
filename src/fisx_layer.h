#ifndef FISX_LAYER_H
#define FISX_LAYER_H
#include <string>
#include <vector>
#include "fisx_elements.h"

/*!
  \class Layer
  \brief Class containing the composition of a layer
*/

class Layer
{
public:

    Layer(const std::string & name="", const double & density = 0.0,
                                       const double & thickness = 0.0,
                                       const double & funnyFactor = 1.0);

    void setMaterial(const std::string & materialName);

    void setMaterial(const Material & material);

    std::string name;
    double density;
    double thickness;

    /*!
    Eventually get a handle to underlying material
    */
    const Material & getMaterial() {return this->material;};

    /*!
    Get the layer transmission at the given energy using the elements library
    supplied.
    If the material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    double getTransmission(const double & energy, Elements & elements) const;

    /*!
    Get the layer transmissions at the given energies using the elements library
    supplied.
    If the material is not defined or cannot be handled, it will throw the
    relevant error.
    */
    std::vector<double> getTransmission(const std::vector<double> & energy, Elements & elements) const;

    /*!
    Return true if material composition was specified.
    Returns false if only the name or formula of the material was given.
    */
    bool hasMaterialComposition() const;

    const std::string & getMaterialName(){return this->materialName;};

private:
    std::string materialName;
    bool hasMaterial;
    Material material;
    double funnyFactor;
};

#endif //FISX_LAYER_H
