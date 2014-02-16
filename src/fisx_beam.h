#ifndef FISX_BEAM_H
#define FISX_BEAM_H
#include <cstddef> // needed for NULL definition!!!
#include <vector>


class Ray
{
public:
    Ray()
    {
        this->energy = 0.0;
        this->weight = 0.0;
        this->characteristic = 0.0;
        this->divergency = 0.0;
    }
    double energy;
    double weight;
    double characteristic;
    double divergency;
    bool operator < (const Ray &b) const
    {
        if (energy < b.energy)
            return true;
        if (weight < b.weight)
            return true;
        return false;
    };
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

    /*! Beam description given as vectors
    */
    void setBeam(const std::vector<double> & energy, \
                 const std::vector<double> & weight,\
                 const std::vector<double> & characteristic,\
                 const std::vector<double> & divergency);

    /*! Easy to wrap interface functions
    Except for the energy, you can use NULL pointers to use default values.
    */
    void setBeam(int nValues, double *energy, double *weight = NULL,
                 double *characteristic = NULL, double *divergency = NULL);

    /*! Returns a vector of length 4 * the number of rays describing the beam as:
    [energy0, weight0, characteristic0, divergency0, energy1, weight1, ...]
    */
    std::vector<double> getBeam();

private:
    bool normalized;
    void normalizeBeam(void);
    std::vector<Ray> rays;
};

#endif //FISX_BEAM_H
