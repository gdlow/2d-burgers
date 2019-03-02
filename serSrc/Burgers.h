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
    double ComputeR(double x, double y);
    double* NextVelocityState(double* Ui, double* Vi, bool U_OR_V);
    void SetMatrixCoefficients();

    /// Burger parameters
    Model* model;
    double* U;
    double* V;
    double* dVel_dx_2_coeffs;
    double* dVel_dy_2_coeffs;
    double* dVel_dx_coeffs;
    double* dVel_dy_coeffs;
    double E;

    /// Term arrays
    double* dVel_2;
    double* dVel;
};
#endif //CLASS_BURGERS
