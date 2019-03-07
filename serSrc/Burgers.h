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
    double* NextVelocityState(bool SELECT_U);
    void SetLinearTerms(double* Vel, double* NextVel);
    void SetNonLinearTerms(double* Vel, double* Other, double* NextVel, bool SELECT_U);
    /// Burger parameters
    Model* model;
    double* U;
    double* V;
    double E;
};
#endif //CLASS_BURGERS
