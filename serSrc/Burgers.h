#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include "Model.h"

/**
 * @class Burgers
 * @brief Creates a Burgers instance that does computations on Burger's equation
 * */
class Burgers {
public:
    Burgers(Model &m);
    ~Burgers();

    void SetInitialVelocity();
    void SetIntegratedVelocity();
    void WriteVelocityFile();
    void SetEnergy();
    double GetE()     const { return E; }
private:
    void ComputeNextVelocityState();
    void wrap(double* A, int Nyr, int Nxr, double** res);

    /// Burger parameters
    Model* model;
    double* U;
    double* V;
    double* NextU;
    double* NextV;
    double E;
};
#endif //CLASS_BURGERS
