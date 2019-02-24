#ifndef CLASS_BURGERS
#define CLASS_BURGERS

#include "Model.h"

class Burgers {
public:
    Burgers(Model &m);
    ~Burgers();

    void SetInitialVelocity();
    void SetIntegratedVelocity();
    void WriteVelocityFile();
    void SetEnergy();
private:
    double ComputeR(double x, double y);
    double NextEnergyState(double* Ui, double* Vi);
    double* NextVelocityState(double* Ui, double* Vi, bool U_OR_V);
    void SetMatrixCoefficients();
    void SetCache(double* Vel, double* Cache);

    // Burger parameters
    Model* model;
    double** U;
    double** V;
    double* U0;
    double* E;
    double* dVel_dx_2_coeffs;
    double* dVel_dy_2_coeffs;
    double* dVel_dx_coeffs;
    double* dVel_dy_coeffs;

    // Cache for dual partition
    double* Vel_c;
};
#endif //CLASS_BURGERS
