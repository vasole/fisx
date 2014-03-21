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
        this->characteristic = 0;
        this->divergency = 0.0;
    }
    double energy;
    double weight;
    int characteristic;
    double divergency;
    bool operator < (const Ray &b) const
    {
        if (characteristic < b.characteristic)
            return true;
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
                 const std::vector<double> & weight = std::vector<double>(),\
                 const std::vector<int> & characteristic = std::vector<int>(),\
                 const std::vector<double> & divergency = std::vector<double>());

    /*!
    Easy to wrap interface functions
    Except for the energy, you can use NULL pointers to use default values.
    */
    void setBeam(int nValues, double *energy, double *weight = NULL,
                 int *characteristic = NULL, double *divergency = NULL);

    void setBeam(int nValues, double *energy, double *weight = NULL,
                 double *characteristic = NULL, double *divergency = NULL);

    /*!
    Returns a constant reference to the internal beam.
    */
    const std::vector<Ray> & getBeam();



    /*! Returns a vector of doubles of length 4 * the number of rays describing the beam as:
    [energy0, weight0, characteristic0, divergency0, energy1, weight1, ...]
    */
    std::vector<double> getBeamAsDoubleVector();

private:
    bool normalized;
    void normalizeBeam(void);
    std::vector<Ray> rays;
};

#endif //FISX_BEAM_H
