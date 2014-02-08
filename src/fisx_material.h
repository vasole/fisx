#ifndef FISX_MATERIAL_H
#define FISX_MATERIAL_H
#include <string>
#include <vector>
#include <map>

/*!
  \class Material
  \brief Class containing the composition of a material

   A material is nothing else than a name and a map of elements and mass fractions.

   The default density, default thickness and comment can be supplied for convenience purposes.
*/
class Material
{
public:
    Material();
    /*!
    Minimalist constructor.
    */
    Material(const std::string &  materialName);
    Material(const std::string & materialName, const double & density, const double & thickness);
    Material(const std::string & materialName, const double & density, \
             const double & thickness, const std::string & comment);

    void setName(const std::string & name);
    void initialize(const std::string & materialName, const double & density,\
                    const double & thickness, const std::string & comment);

    /*!
    Set the composition of the material.
    This method normalizes the supplied amounts to make sure the sum is one.
    */
    void setComposition(const std::map<std::string, double> &);

    /*!
    Alternative method to set the composition of the material
    This method normalizes the supplied amounts to make sure the sum is one.
    */
    void setComposition(const std::vector<std::string> &, const std::vector<double> &);

    /*!
    Return the material composition as normalized mass fractions
    */
    std::map<std::string, double> getComposition();
    std::string getName();
    std::string getComment();
    double getDefaultDensity();
    double getDefaultThickness();

private:
    std::string name;
    bool initialized;
    std::map<std::string, double> composition;
    double defaultDensity;
    double defaultThickness;
    std::string    comment;
};

#endif //FISX_MATERIAL_H
