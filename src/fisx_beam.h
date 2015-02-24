#ifndef FISX_BEAM_H
#define FISX_BEAM_H
#include <cstddef> // needed for NULL definition!!!
#include <vector>
#include <iostream>

namespace fisx
{

struct Ray
{
    double energy;
    double weight;
    int characteristic;
    double divergency;

    Ray()
    {
        energy = 0.0;
        weight = 0.0;
        characteristic = 0;
        divergency = 0.0;
    }

    bool operator < (const Ray &b) const
    {
        return (energy < b.energy);
    }
};

/*!
  \class Beam
  \brief Class describing an X-ray beam

   At this point a beam is described by a set of energies and weights. The characteristic flag just indicates if
   it is an energy to be considered for calculation of scattering peaks.
*/
class Beam
{
public:
    Beam();
    /*!
    Minimalist constructor.
    */

    /*!
    Beam description given as vectors.
    The beam is always internally ordered.
    */
    void setBeam(const std::vector<double> & energy, \
                 const std::vector<double> & weight = std::vector<double>(),\
                 const std::vector<int> & characteristic = std::vector<int>(),\
                 const std::vector<double> & divergency = std::vector<double>());

    friend std::ostream& operator<< (std::ostream& o, Beam const & beam);

    /*!
    Easy to wrap interface functions
    Except for the energy, you can use NULL pointers to use default values.
    */
    void setBeam(const int & nValues, const double *energy, const double *weight = NULL,
                 const int *characteristic = NULL, const double *divergency = NULL);

    void setBeam(const int & nValues, const double *energy, const double *weight = NULL,
                 const double *characteristic = NULL, const double *divergency = NULL);

    void setBeam(const double & energy, const double divergency = 0.0);

    /*!
    Returns a constant reference to the internal beam.
    */
    const std::vector<Ray> & getBeam();

    /*!
    Currently it returns a vector of "elements" in which each element is a vector of
    doubles with length equal to the number of energies.
    The first four elements are warranteed to exist and they are:
    [energy0, energy1, energy2, ...]
    [weight0, weight1, weight2, ...]
    [characteristic0, characteristic1, charactersitic2, ...]
    [divergency0, divergency1, divergency2, ...]
    */
    std::vector<std::vector<double> > getBeamAsDoubleVectors() const;

private:
    bool normalized;
    void normalizeBeam(void);
    std::vector<Ray> rays;
};

} // namespace fisx

#endif //FISX_BEAM_H
